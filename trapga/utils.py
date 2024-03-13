import numpy as np
import tqdm
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from collections import Counter
import pandas as pd
import warnings

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


def extract_similar(genes, arr):
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
    for _ in range(10):
        kmeans = KMeans(n_clusters=round(arr.shape[0]/100)).fit_predict(df_sorted.iloc[:,1:].to_numpy())
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
    return [key for key, value in similars.items() if value >= 7]

def extract_coexpressed(genes, arr):
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
    return {ex_gene: [v for k, v in matching_pairs if k == ex_gene] for ex_gene, _ in matching_pairs}