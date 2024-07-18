import numpy as np
import scipy.sparse
import tqdm
from sklearn.cluster import KMeans
from collections import Counter
import pandas as pd
import warnings
import scipy
import os
import time
from setga import utils, select_subset
from Bio import SeqIO
from functools import partial
from scipy.sparse import csr_matrix
import pickle


def comp_vars_sampled(expression_data,random_number_generator,rounds,phylostrata):
    rows_A = rounds
    # Generate the predefined matrix B
    matrix_B = expression_data.expressions_n
    # Define a generator to generate rows of matrix A
    def generate_rows_A():
        for _ in range(rows_A):
            yield phylostrata[random_number_generator(expression_data.full.shape[0])]

    # Initialize the result matrix
    result = np.zeros((rows_A, expression_data.expressions_n.shape[1]))

    # Generate rows of matrix A and perform matrix multiplication with matrix B
    for i, row_A in tqdm.tqdm(enumerate(generate_rows_A())):
        result[i] = np.dot(row_A, matrix_B)
    return np.var(result/expression_data.expressions_n.sum(axis=0),axis=1)
        

def comp_vars(expression_data,rounds):
    """Computes the min-variances of a TAI patterns for permuted phylostrata

    :param expression_data: expression data
    :type expression_data: pd.DataFrame
    :param rounds: number of permutations of phylostrata
    :type rounds: int
    :return: variances for the TAI patterns, used to determine the empirical p-value
    :rtype: np.array
    """
    avgs = []
    phil = expression_data.full["Phylostratum"]
    print("Running permuations")
    for _ in tqdm.trange(rounds):
        perm = np.random.permutation(phil)
        weighted = expression_data.expressions.mul(perm, axis=0)
        avg = weighted.sum(axis=0)/expression_data.expressions_n_sc.sum(axis=0)
        avgs.append(avg)
    return np.var(avgs, axis=1)

def compute_permutation_variance_sc(expression_data, rounds):
    # Precompute the sum of expressions_n along axis 0
    expressions_n_sum = expression_data.expressions_n.sum(axis=0)
    # Ensure it's a 1D array for division later
    expressions_n_sum = np.array(expressions_n_sum).flatten()

    avgs = []
    phil = expression_data.full["Phylostratum"]
    
    print("Running permutations")
    for _ in tqdm.trange(rounds):
        # Generate a permutation of Phylostratum
        perm = np.random.permutation(phil)
        
        # Create a diagonal matrix from the permutation
        perm_diag = scipy.sparse.diags(perm)
        
        # Multiply the expressions_n matrix by the permuted diagonal matrix
        weighted = perm_diag.dot(expression_data.expressions_n)
        
        # Compute the weighted sum along axis 0
        weighted_sum = weighted.sum(axis=0)
        
        # Compute the average
        avg = np.array(weighted_sum / expressions_n_sum).flatten()
        avgs.append(avg)
    
    # Convert avgs to a numpy array for variance computation
    avgs = np.array(avgs)
    print(avgs.shape)
    
    # Compute and return the variance along axis 1
    return np.var(avgs, axis=1)

def extract_similar(args):
    """identified genes, that have similar expression patterns to the extracted ones and generates a file that includes those and the genes identified by the minimizer

    :param args: gets args from the main cli script
    :type args: argparse.Namespace
    :return: saves the similar genes and the originally identified genes into a file
    :rtype: None
    """
    genes = np.array([line.strip() for line in open(args.genes, 'r')])
    arr = pd.read_csv(args.input, delimiter="\t")

    def remove_one_type_clusters(clusters):
        """removes clusters with just one type of genes (extracted/not extracted) 

        :param clusters: set of co-clustered genes
        :type clusters: list
        """
        def same_type_clusters(clusters):
            """tests if there is just one type of genes (extracted/not extracted) in the set

            :param clusters: set of co-clustered genes
            :type clusters: list
            :return: True if there is just one type of genes
            :rtype: bool
            """
            types = set(["ext" if x in genes else "edg" for x in clusters])
            return len(types) == 1

        valid_clusts = []      
        for clust in clusters:
            if not same_type_clusters(clust):
                valid_clusts.append(clust)
        return valid_clusts
    

    df_sorted = arr 
    df_sorted= df_sorted.reindex(columns=["GeneID","Phylostratum"] + list(df_sorted.columns[2:]))
    similars = []
    runs = 5
    for _ in tqdm.trange(runs):
        kmeans = KMeans(n_clusters=round(arr.shape[0]/100),n_init = 5).fit_predict(df_sorted.iloc[:,1:].to_numpy())
        clusters = df_sorted.GeneID.groupby(kmeans).apply(list)

        valid_clusts = remove_one_type_clusters(clusters)
        similar = []
        for cluster in valid_clusts: 

            clust = arr[arr.GeneID.isin(cluster)]
            clust.set_index('GeneID', inplace=True)
            corr = clust.iloc[:,2:].T.corr()

            ex_genes = list(set(cluster).intersection(set(genes)))

            phylostratum_threshold = 1
            correlation_threshold = 0.95

            def is_close(value, target_value, threshold):
                """returns if a value is close to a trashold value by a given threshold

                :param value: given value
                :type value: float
                :param target_value: target value
                :type target_value: float
                :param threshold: tolerance threshold
                :type threshold: fload
                :return: true if a value is close to a trashold value by a given threshold
                :rtype: bool
                """
                return abs(value - target_value) <= threshold
            
            for id_to_check in cluster:
                target_phylostratum = clust.loc[clust.index == id_to_check, 'Phylostratum'].iloc[0]
                close_phylostratum_rows = clust[clust.index.isin(ex_genes) & clust['Phylostratum'].apply(lambda x: is_close(x, target_phylostratum, phylostratum_threshold))]
                
                if not close_phylostratum_rows.empty:
                    max_corr_id = corr.loc[id_to_check, close_phylostratum_rows.index].idxmax()
                    correlation_value = corr.loc[id_to_check, max_corr_id]
                    if correlation_value > correlation_threshold:
                        if id_to_check not in genes:
                            similar.append(id_to_check)
        similars.append(similar)
    similars = dict(Counter([item for similar in similars for item in similar]))
    add_genes = np.array([key for key, value in similars.items() if value >= runs * 0.7])
    np.savetxt(os.path.join(args.output,"extracted_genes_added.txt"),np.concatenate([genes, add_genes]), fmt="%s")


def extract_coexpressed(args):
    """Finds all genes, that are co-expressed with the identified set and saves them in a file

    :param args: gets args from the main cli script
    :type args: argparse.Namespace
    """
    genes = np.array([line.strip() for line in open(args.genes, 'r')])
    arr = pd.read_csv(args.input, delimiter="\t")
    pearson_threshold = 30
    if arr.shape[1] < pearson_threshold + 2:
        warnings.warn(f"Cannot analyze coexpression for less than {pearson_threshold} stages")
        return
    exps = arr.iloc[:, 2:]
    exps = exps[exps.apply(lambda row: np.nanmax(row.values) >= 100, axis=1)]
    pg = arr.loc[exps.index, ['Phylostratum',"GeneID"]]
    arr = pd.concat([pg, exps], axis=1)

    arr['GeneID'] = pd.Categorical(arr['GeneID'], categories=list(set(genes)) + list(set(arr.GeneID).difference(set(genes))), ordered=True)

    # Sort the DataFrame based on column 'B'
    df_sorted = arr.sort_values(by='GeneID')
    df_sorted=df_sorted.reindex(columns=["GeneID","Phylostratum"] + list(df_sorted.columns[2:]))
    df_sorted.set_index('GeneID', inplace=True)
    corr = df_sorted.iloc[:,2:].T.corr(method='pearson')
    cross_cor = corr.iloc[len(genes) :,:len(genes)]
    matching_pairs = cross_cor.stack()[cross_cor.stack() > 0.95].index.tolist()
    ex_genes =  {ex_gene: [v for k, v in matching_pairs if k == ex_gene] for ex_gene, _ in matching_pairs}
    arrays = [(key, np.array(ex_genes[key])) for key in ex_genes]
    coexpressed = np.concatenate([np.column_stack((np.full_like(arr[1], arr[0]), arr[1])) for arr in arrays])
    df = pd.DataFrame(coexpressed,columns=["extracted_genes", "coexpressed"])
    df.to_csv(os.path.join(args.output,"coexpressed.tsv"),sep="\t")



def get_extracted_genes(args):
    """extracts genes, that are significantly influencing the TAI pattern

    :param args: gets args from the main cli script
    :type args: argparse.Namespace
    :return: Saves the identified genes, the (best) solution and run summary into files 
    :rtype: None
    """
    class Expression_data:
        """class to store the expression dataset with some precomputations
        """

        def quantilerank(xs):
            """computes the quantile rank for the phylostrata

            :param xs: numpy array of values
            :type xs: np.array
            :return: quantile ranked values
            :rtype: np.array
            """
            ranks = scipy.stats.rankdata(xs, method='average')
            quantile_ranks = [scipy.stats.percentileofscore(ranks, rank, kind='weak') for rank in ranks]
            return np.array(quantile_ranks)/100

        def __init__(self,expression_data) -> None:
            """
            :param expression_data: expression dataset
            :type expression_data: pd.DataFrame
            """
            expression_data["Phylostratum"] = Expression_data.quantilerank(expression_data["Phylostratum"])
            self.full = expression_data
            exps = expression_data.iloc[:, 2:]
            #exps = exps.applymap(lambda x: np.sqrt(x))
            #exps = exps.applymap(lambda x: np.log(x + 1))
            age_weighted = exps.mul(expression_data["Phylostratum"], axis=0).to_numpy()
            self.age_weighted = age_weighted
            self.expressions_n = exps.to_numpy()
            self.expressions = exps
            self.weighted_sum = np.sum(exps.mul(expression_data["Phylostratum"], axis=0).to_numpy(),axis=0)
            self.exp_sum = np.sum(exps.to_numpy(),axis=0)
            self.expressions_n_sc = exps.to_numpy()
            if args.single_cell:
                self.expressions_n = csr_matrix(exps.to_numpy())  # Define your sparse matrix 'a'
                self.age_weighted = csr_matrix(age_weighted)  # Define your sparse matrix 'a_w'
            print(exps.shape)


    arr = pd.read_csv(args.input,
                    delimiter="\t")
    expression_data = Expression_data(arr)

    def compute_distance(expression_data,n_sol):
    # Ensure n_sol is a dense array
        n_sol = np.asarray(n_sol)

        # Stage-wise sum of expression of selected genes
        a_w_result = n_sol @ expression_data.age_weighted
        a_result = n_sol @ expression_data.expressions_n

        # "Removing" the selected genes from the dataset
        numerator = expression_data.weighted_sum - a_w_result
        denominator = expression_data.exp_sum - a_result
        return np.var(np.divide(numerator, denominator))


    if args.variances:
        permuts = np.loadtxt(args.variances)
    else:
        #permuts = comp_vars_sampled(expression_data,phylostrata_sampler,10000,phylostrata)
        perm_start = time.perf_counter()
        if args.single_cell:
            
            permuts = compute_permutation_variance_sc(expression_data,1000)
            
        else:
            permuts = comp_vars(expression_data,100000)
        perm_stop = time.perf_counter()
    

    ind_length = expression_data.full.shape[0]



    def get_distance(solution):
        """computes variance of the TAI for the particular solution

        :param solution: binary encoded, which genes belong in the solution
        :type solution: array
        :return: variance
        :rtype: float
        """
        
        sol = np.array(solution)
        sol = np.logical_not(sol).astype(int)
        up = sol.dot(expression_data.age_weighted)
        down = sol.dot(expression_data.expressions_n)
        avgs = np.divide(up,down)
        return np.var(avgs)

    if args.single_cell:
        max_value = compute_distance(expression_data,np.zeros(ind_length))
    else:
        max_value = get_distance(np.zeros(ind_length))



    def end_evaluate_individual(individual):
        """individual fitness without the cutoff, just pure p-value

        :param individual: binary encoded, which genes belong in the solution
        :type individual: array
        :return: fitness
        :rtype: float
        """
        individual = np.array(individual)
        if args.single_cell:
            distance = compute_distance(expression_data,individual)
        else:
            distance = get_distance(individual)
        fit =  np.count_nonzero(permuts < distance)/len(permuts)
        # Return the fitness values as a tuple
        return np.sum(individual), fit

        
    def evaluate_individual(individual,permuts,expression_data):
        """computes the overall fitness of an individual

        :param individual: binary encoded, which genes belong in the solution
        :type individual: array
        :param permuts: precomputed variances from flat-line test
        :type permuts: np.array
        :param expression_data: dataset of expression of the genes
        :type expression_data: pd.DataFrame
        """
        def get_fit(res):
            """computes empirical p-value of an individual

            :param res: variance of an individual
            :type res: np.array
            :return: empirical p-value 
            :rtype: float
            """
            p = np.count_nonzero(permuts < res)/len(permuts)
            r = (res) / (max_value)
            r = r + p
            return r if p > 0.2 else 0
        sol = np.array(individual)
        sol = np.logical_not(sol).astype(int)
        distance = np.var(np.divide(sol.dot(expression_data.age_weighted),sol.dot(expression_data.expressions_n)))
        fit = get_fit(distance)
        # Return the fitness values as a tuple
        return fit
    

    def evaluate_individual_sc(individual,permuts,expression_data):
        """computes the overall fitness of an individual

        :param individual: binary encoded, which genes belong in the solution
        :type individual: array
        :param permuts: precomputed variances from flat-line test
        :type permuts: np.array
        :param expression_data: dataset of expression of the genes
        :type expression_data: pd.DataFrame
        """
        def get_fit(res):
            """computes empirical p-value of an individual

            :param res: variance of an individual
            :type res: np.array
            :return: empirical p-value 
            :rtype: float
            """
            p = np.count_nonzero(permuts < res)/len(permuts)
            r = (res) / (max_value)
            r = r + p
            return r if p > 0.2 else 0
        sol = np.array(individual)
        if args.single_cell:
            distance = compute_distance(expression_data,sol)
        else:
            distance = get_distance(sol)
        fit = get_fit(distance)
        # Return the fitness values as a tuple
        return [fit]
    

    def get_skewed_reference(num_points, skew):
        y_values = np.linspace(skew, 1, num_points+1)
        # Calculate corresponding y values such that the sum of x and y is 1
        x_values = 1 - y_values

        # Create the numpy array with two columns
        return np.column_stack((x_values, y_values))[:-1]
    
    def get_uniform_reference(num_points):
        y_values = np.linspace(0, 1, num_points)
        # Calculate corresponding y values such that the sum of x and y is 1
        x_values = 1 - y_values

        # Create the numpy array with two columns
        return np.column_stack((x_values, y_values))

    ref_points = get_uniform_reference(10)
    ref_points = np.append(ref_points,get_skewed_reference(4,0.75)[:-1],axis=0)

    population_size = 150
    num_generations = 15000
    num_islands = 4
    if args.single_cell:
        population_size = 80
        num_islands = 3
    mut  = 0.005
    cross = 0.02
    stop_after = 200
    if args.single_cell:
        stop_after = 40

    tic = time.perf_counter()
    if args.single_cell:
        evaluation_function = evaluate_individual_sc
    else:
        evaluation_function = evaluate_individual
    eval_part = partial(evaluation_function, permuts = permuts, expression_data = expression_data)
    pop,_,gens,logbook, best_sols = select_subset.run_minimizer(expression_data.full.shape[0],eval_part,1,"Variance",
                    mutation_rate = mut,crossover_rate = cross, 
                    pop_size = population_size, num_gen = num_generations, num_islands = num_islands, mutation = ["weighted","weighted","bit-flip","bit-flip"], 
                    crossover =  "uniform",
                    selection = "NSGA3",frac_init_not_removed = 0.2,ref_points = ref_points, stop_after = stop_after,weights = np.sqrt(np.var(expression_data.expressions_n,axis=1)))

    toc = time.perf_counter()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    ress = np.array([end_evaluate_individual(x) for x in pop])
    min_ex = min(pop, key=lambda ind: ind.fitness.values[1]).fitness.values[0]
    pop = np.array(pop)
    genes = utils.get_results(pop,ress,expression_data.full.GeneID)
    np.savetxt(os.path.join(args.output,"extracted_genes.txt"),genes, fmt="%s")

    if args.save_stats:
        np.savetxt(os.path.join(args.output,"permuts.txt"),permuts)
        with open(os.path.join(args.output, "summary.txt"), 'w') as file:
            # Write the first line
            file.write(f'Time: {toc - tic:0.4f} seconds\n')
            file.write(f'Permutation Time: {perm_stop - perm_start:0.4f} seconds\n')
            # Write the second line
            file.write(f'Min. genes in sol: {min_ex}\n')
            file.write(f'Number of genes: {len(genes)}\n')
            file.write(f'Number of generations: {gens}\n')
            np.savetxt(os.path.join(args.output,"best_sols.csv"), best_sols, delimiter="\t")
            np.savetxt(os.path.join(args.output,"complete.csv"), np.array(pop), delimiter="\t")

            with open(os.path.join(args.output, "logbook.pickle"), 'wb') as file:
                pickle.dump(logbook, file)





    

def get_fastas(args):
    """Makes a fasta file with all the extracted genes

    :param args: gets args from the main cli script
    :type args: argparse.Namespace
    """
    genes = np.array([line.strip() for line in open(args.genes, 'r')])
    filtered_records = []
    with open(args.fastas, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if any(record.id.startswith(gene) for gene in genes):
                filtered_records.append(record)

    with open(os.path.join(args.output,"extracted_fastas.fasta"), "w") as output_file:
        SeqIO.write(filtered_records, output_file, "fasta")
