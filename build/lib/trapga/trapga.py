#!usr/bin/env python3 
import numpy as np
import scipy
import pandas as pd
from setminga import select
from setminga import utils
from trapga import utils as hg_utils
import argparse
import os
import time

parser = argparse.ArgumentParser()
    
# Add arguments for input and output files
parser.add_argument("input", type=str, help="Input file path")
parser.add_argument("output", type=str, help="Output file path")
parser.add_argument("--variances", type=str, help="Precomputed variances file stored as numbers in one line of", nargs='?')
parser.add_argument("--save_plot", action="store_true", help="Save pareto plot")
args = parser.parse_args()

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
    permuts = np.loadtxt(args.permuts)
else:
    permuts = hg_utils.comp_vars(expression_data,100000)

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
        return r if p > 0.1 else 0
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
                  mutation_rate = 0.001,crossover_rate = 0.02, 
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