import numpy as np
import tqdm
import matplotlib.pyplot as plt
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
        avg = weighted.sum(axis=0)/expression_data.expressions_n.sum(axis=0)
        avgs.append(avg)
    return np.var(avgs, axis=1)

def comp_min_max(expression_data,rounds):
    """Computes the min-max value of a TAI patterns for permuted phylostrata

    :param expression_data: expression data
    :type expression_data: pd.DataFrame
    :param rounds: number of permutations of phylostrata
    :type rounds: int
    :return: min-max values for the TAI patterns, used to determine the empirical p-value
    :rtype: np.array
    """
    avgs = []
    phil = expression_data.full["Phylostratum"]
    print("Running permuations")
    for _ in tqdm.trange(rounds):
        perm = np.random.permutation(phil)
        weighted = expression_data.expressions.mul(perm, axis=0)
        avg = weighted.sum(axis=0)/expression_data.expressions_n.sum(axis=0)
        avgs.append(avg)
    return np.max(avgs, axis=1) - np.min(avgs, axis=1)


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

    # Concatenate the arrays


def get_extracted_genes(args):
    """extracts genes, that are significantly influencing the TAI pattern

    :param args: gets args from the main cli script
    :type args: argparse.Namespace
    :return: Saves the identified genes, the (best) solution and run summary into files 
    :rtype: None
    """
    def mutWeightedFlipBit(individual, weights, mutation_rate):
        """Mutate the input individual by flipping the value of its attributes based on weighted probabilities for each index.

        The `individual` is expected to be a sequence, and the values of the attributes shall stay valid after the `not` operator is called on them. The overall mutation rate is preserved.

        :param individual: list
            Individual to be mutated.
        :type individual: list
        :param weights: list
            List of weights for each index of the individual. The higher the weight, the more probable the mutation.
        :type weights: list
        :param mutation_rate: float
            Overall mutation rate for the individual.
        :type mutation_rate: float

        :returns:
            tuple
                A tuple containing the mutated individual.

        :notes:
            This function uses the `numpy.random.choice` function from the numpy library to select indices based on their weights and mutates them with the specified mutation rate.
        """
        # Calculate the number of mutations based on the mutation rate
        num_mutations = int(len(individual) * mutation_rate)

        # Select indices to mutate based on weights
        indices_to_mutate = np.random.choice(len(individual), size=num_mutations, replace=False, p=weights)

        # Mutate selected indices
        for i in indices_to_mutate:
            individual[i] = type(individual[i])(not individual[i])

        return individual,

    def cxWeightedUniform(ind1, ind2, weights,cx_rate):
        """Executes a weighted uniform crossover that modifies in place the two
        sequence individuals.

        The attributes are swapped according to the *weights* probability.

        :param ind1: The first individual participating in the crossover.
        :type ind1: list
        :param ind2: The second individual participating in the crossover.
        :type ind2: list
        :param weights: List of weights for each attribute to be exchanged.
                        The higher the weight, the more probable the exchange.
        :type weights: list
        :returns: A tuple of two individuals.
        :rtype: tuple

        This function uses the :func:`numpy.random.choice` function from the numpy
        library.
        """
        num_cross = int(len(ind1) * cx_rate)
        crossover_indices = np.random.choice(len(ind2), size=num_cross, p=weights)
        for i in crossover_indices:
            ind1[i], ind2[i] = ind2[i], ind1[i]

        return ind1, ind2

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
            self.age_weighted = exps.mul(expression_data["Phylostratum"], axis=0).to_numpy()
            self.expressions_n = exps.to_numpy()
            self.expressions = exps


    arr = pd.read_csv(args.input,
                    delimiter="\t")
    expression_data = Expression_data(arr)
    if args.variances:
        permuts = np.loadtxt(args.variances)
    else:
        permuts = comp_vars(expression_data,100000)

    ind_length = expression_data.full.shape[0]

    population_size = 150
    num_generations = 8000
    num_islands = 7


    def get_distance(solution):
        """computes variance of the TAI for the particular solution

        :param solution: binary encoded, which genes belong in the solution
        :type solution: array
        :return: variance
        :rtype: float
        """
        sol = np.array(solution)
        up = sol.dot(expression_data.age_weighted)
        down = sol.dot(expression_data.expressions_n)
        avgs = np.divide(up,down)
        return np.var(avgs)


    max_value = get_distance(np.ones(ind_length))



    def end_evaluate_individual(individual):
        """individual fitness without the cutoff, just pure p-value

        :param individual: binary encoded, which genes belong in the solution
        :type individual: array
        :return: fitness
        :rtype: float
        """
        individual = np.array(individual)
        num_not_removed = np.sum(individual)
        len_removed = ind_length - num_not_removed
        distance = get_distance(individual)
        fit =  np.count_nonzero(permuts < distance)/len(permuts)
        # Return the fitness values as a tuple
        return len_removed, fit

        
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
            return r if p > 0.1 else 0
        sol = np.array(individual)
        distance = np.var(np.divide(sol.dot(expression_data.age_weighted),sol.dot(expression_data.expressions_n)))
        fit = get_fit(distance)
        # Return the fitness values as a tuple
        return [fit]

    mut  = 0.001
    cross = 0.02
    weights = np.log(np.var(expression_data.expressions_n,axis=1) + 1)
    mutation_part = partial(mutWeightedFlipBit,weights=weights/np.sum(weights),mutation_rate = mut)
    crossover_part = partial(cxWeightedUniform,weights=weights/np.sum(weights),cx_rate = cross)
    tic = time.perf_counter()
    pop,pareto_front = select_subset.run_minimizer(expression_data.full.shape[0],evaluate_individual,1,["Variance"], 
                    eval_func_kwargs={"permuts": permuts, "expression_data": expression_data},
                    mutation_rate = mut,crossover_rate = cross, 
                    pop_size = population_size, num_gen = num_generations, num_islands = num_islands, mutation = mutation_part, 
                    crossover =  crossover_part,
                    selection = "SPEA2",frac_init_not_removed = 0.005)

    toc = time.perf_counter()
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    np.savetxt(os.path.join(args.output,"complete.csv"), np.array(pop), delimiter="\t")
    ress = np.array([end_evaluate_individual(x) for x in pop])

    pop = np.array(pop)
    par = np.array([list(x) for x in pareto_front[0]])
    parr = np.array([end_evaluate_individual(x) for x in par])

    np.savetxt(os.path.join(args.output,"pareto.csv"), par, delimiter="\t")


    if args.save_plot:
        plot = utils.plot_pareto(ress,parr)
        plot.savefig(os.path.join(args.output, "pareto_front.png")) 
    genes = utils.get_results(pop,ress,expression_data.full.GeneID)
    np.savetxt(os.path.join(args.output,"extracted_genes.txt"),genes, fmt="%s")

    with open(os.path.join(args.output, "summary.txt"), 'w') as file:
        # Write the first line
        file.write(f'Time: {toc - tic:0.4f} seconds\n')
        
        # Write the second line
        file.write(f'Number of genes: {len(genes)}\n')

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