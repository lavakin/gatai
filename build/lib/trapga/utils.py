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
from setminga import utils, select
from Bio import SeqIO

def comp_vars(expression_data,rounds):
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
    genes = np.array([line.strip() for line in open(args.genes, 'r')])
    arr = pd.read_csv(args.input, delimiter="\t")

    def remove_one_type_clusters(clusters):
        def same_type_community(community):
            types = set(["ext" if x in genes else "edg" for x in community])
            return len(types) == 1

        valid_clusts = []      
        for clust in clusters:
            if not same_type_community(clust):
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
    class Expression_data:

        def quantilerank(xs):
            ranks = scipy.stats.rankdata(xs, method='average')
            quantile_ranks = [scipy.stats.percentileofscore(ranks, rank, kind='weak') for rank in ranks]
            return np.array(quantile_ranks)/100

        def __init__(self,expression_data) -> None:
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
    #parents_ratio = 0.2
    num_generations = 8000
    init_num_removed = 150
    num_islands = 6


    def get_distance(solution):
        sol = np.array(solution)
        up = sol.dot(expression_data.age_weighted)
        down = sol.dot(expression_data.expressions_n)
        avgs = np.divide(up,down)
        return np.var(avgs)


    max_value = get_distance(np.ones(ind_length))



    def end_evaluate_individual(individual):
        individual = np.array(individual)
        num_not_removed = np.sum(individual)
        len_removed = ind_length - num_not_removed
        distance = get_distance(individual)
        fit =  np.count_nonzero(permuts < distance)/len(permuts)
        # Return the fitness values as a tuple
        return len_removed, fit

        
    def evaluate_individual(individual,permuts,expression_data):
        def get_fit(res):
            p = np.count_nonzero(permuts < res)/len(permuts)
            r = (res) / (max_value)
            r = r + p
            return r if p > 0.2 else 0
        sol = np.array(individual)
        distance = np.var(np.divide(sol.dot(expression_data.age_weighted),sol.dot(expression_data.expressions_n)))
        fit = get_fit(distance)
        # Return the fitness values as a tuple
        return [fit]

    mut  = 0.001
    cross = 0.02
    tic = time.perf_counter()
    pop,pareto_front = select.run_minimizer(expression_data.full.shape[0],evaluate_individual,1,["Variance"], 
                    eval_func_kwargs={"permuts": permuts, "expression_data": expression_data},
                    mutation_rate = mut,crossover_rate = cross, 
                    pop_size = 150, num_gen = num_generations, num_islands = 8, mutation = "bit_flip" , 
                    crossover =  "uniform_partialy_matched",
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
        plot = utils.plot_pareto(ress,parr,args.output)
        plot.savefig(os.path.join(args.output, "pareto_front.png")) 
    genes = utils.get_results(pop,ress,args.output,expression_data.full.GeneID)
    np.savetxt(os.path.join(args.output,"extracted_genes.txt"),genes, fmt="%s")

    with open(os.path.join(args.output, "summary.txt"), 'w') as file:
        # Write the first line
        file.write(f'Time: {toc - tic:0.4f} seconds\n')
        
        # Write the second line
        file.write(f'Number of genes: {len(genes)}\n')

def get_fastas(args):
    genes = np.array([line.strip() for line in open(args.genes, 'r')])
    filtered_records = []
    with open(args.fastas, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in genes:
                filtered_records.append(record)

    with open(os.path.join(args.output,"extracted_fastas.fasta"), "w") as output_file:
        SeqIO.write(filtered_records, output_file, "fasta")