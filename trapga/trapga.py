#!usr/bin/env python3 
from trapga import utils as hg_utils
import argparse


parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(dest='command')

# Sub-parser for the add_genes command
add_genes_parser = subparsers.add_parser('add_genes', help='Add genes functionality')
add_genes_parser.add_argument("input", type=str, help="Input file path to the dataset")
add_genes_parser.add_argument("genes", type=str, help="Input file path to the extracted genes")
add_genes_parser.add_argument("output", type=str, help="Output file path")

# Sub-parser for the run_minimizer command
run_minimizer_parser = subparsers.add_parser('run_minimizer', help='Run minimizer functionality')
run_minimizer_parser.add_argument("input", type=str, help="Input file path")
run_minimizer_parser.add_argument("output", type=str, help="Output file path")
run_minimizer_parser.add_argument("--variances", type=str, help="Precomputed variances file stored as numbers in one line", nargs='?')
run_minimizer_parser.add_argument("--save_plot", action="store_true", help="Save pareto plot")

# Sub-parser for the find_coexpressed command
find_coexpressed_parser = subparsers.add_parser('find_coexpressed', help='Find coexpressed genes functionality')
find_coexpressed_parser.add_argument("input", type=str, help="Input file path to the dataset")
find_coexpressed_parser.add_argument("genes", type=str, help="Input file path to the extracted genes")
find_coexpressed_parser.add_argument("output", type=str, help="Output file path")

# Sub-parser for the get_fastas command
get_fastas_parser = subparsers.add_parser('get_fastas', help='Get fastas functionality')
get_fastas_parser.add_argument("genes", type=str, help="Input gene list")
get_fastas_parser.add_argument("fastas", type=str, help="Input fasta file path")
get_fastas_parser.add_argument("output", type=str, help="Output file path")

args = parser.parse_args()

if args.command == "run_minimizer":
    hg_utils.get_extracted_genes(args)
elif args.command == "add_genes":
    hg_utils.extract_similar(args)
elif args.command == "find_coexpressed":
    hg_utils.extract_coexpressed(args)
elif args.command == "get_fastas":
    hg_utils.get_fastas(args)

def cli():
    """just some mandatory thing for it to work
    """
    pass

if __name__ == '__main__':
    """main
    """
    cli()
