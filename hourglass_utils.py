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

def plot(GA):
        fig11 = plt.figure(figsize=(13, 8))
        x = np.linspace (0.0000005, 0.000001, 300) 

        #calculate pdf of Gamma distribution for each x-value
        y = GA.gamma.pdf(x)
        best_inds = np.argsort(GA.fitness)[-15:]
        best = best_inds[-1]
        best_len = GA.ind_length - np.sum(GA.population[best]) 
        fitness = GA.fitness[best_inds]
        best_sols = GA.population[best_inds]
        col = [150 if x == max(fitness) else 40 for x in fitness]

        plt.style.use('seaborn-v0_8-pastel')
        fig =plt.figure(GA.curr_gen,figsize=(13, 8))
        grid = fig.add_gridspec(3, 10, wspace=1.5, hspace=0.5)
        ax1 = plt.subplot(grid[:2, :-1])
        ax2 = plt.subplot(grid[:, -1])
        ax3 = plt.subplot(grid[2, :-1])
        varr,_ = GA.get_var_and_p_single(GA.population[best])

        for sol,f in zip(GA.get_avgs(best_sols),col):
            ax1.plot(["Zygote", "Quadrant","Globular","Heart","Torpedo","Bent","Mature"], sol, lw=3,c=plt.cm.Greens(f))
            ax1.plot(["Zygote", "Quadrant","Globular","Heart","Torpedo","Bent","Mature"], GA.get_avgs(np.ones(GA.ind_length)), lw=3,c="Grey")
            ax1.set_ylim([3, 3.42])
            ax1.set_xlabel("Stage")
            ax1.set_ylabel("TAI")
        
        ax2.bar(["Removed"], best_len,color="Green")
        ax2.set_ylim([0, 900])
        ax2.yaxis.tick_right()
        ax3.plot(x, y)
        ax3.set_xlabel("Variance")
        ax3.set_ylabel("pdf")
        if varr  < 0.00000194:
            plt.vlines(varr, 0, GA.gamma.pdf(varr),colors=["Green"])
        plt.savefig("./best_graphs/best" + str(GA.curr_gen) + ".png")
        plt.close()

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