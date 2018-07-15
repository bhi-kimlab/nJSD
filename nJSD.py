#!/usr/bin/python

import argparse, sys
from tITH import simple_nJSD as run

def f_parser(option):
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument('-n', action='store', dest='network',              \
                        help='Location to network file: geneA    geneB')
    pre_parser.add_argument('-r', action='store', dest='r_GEP',                \
                        help='File name of Refernece gene-expression profile')
    pre_parser.add_argument('-i', action='store', dest='q_GEP',                \
                        help='File name of Query gene-expression profile')
    
    if option == "whole":
        whole_parser = argparse.ArgumentParser(parents=[pre_parser])
        return whole_parser

    elif option == "geneset":
        geneset_parser = argparse.ArgumentParser(parents=[pre_parser])
        geneset_parser.add_argument('-t', action='store', dest='geneset',      \
                        help='Add Gene set list')
        return geneset_parser

    else:
        ''' unreachable '''
        exit("Unexpected termination: arg parser error")

if __name__ == "__main__":
    run_type = set(["whole", "geneset"])
    assert sys.argv[1] in run_type,                                            \
            "\n Only 'whole' or 'geneset' is available."

    result = f_parser(sys.argv[1]).parse_args(sys.argv[2:]) 

    if sys.argv[1] == "whole":
        run.whole_cal(result.network, result.r_GEP, result.q_GEP)
    elif sys.argv[1] == "geneset":
        run.gene_set_cal(result.network, result.r_GEP, result.q_GEP,           \
                         result.geneset)
