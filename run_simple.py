#!/usr/bin/python

import argparse, sys
from tITH import simple_nJSD as run

def f_parser():
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-n', action='store', dest='network', help='Location to network file: geneA	geneB')
	parser.add_argument('-r', action='store', dest='r_GEP', help='File name of Refernece gene-expression profile')
	parser.add_argument('-i', action='store', dest='q_GEP', help='File name of Query gene-expression profile')


	return parser

if __name__ == "__main__":
	result = f_parser().parse_args(sys.argv[1:]) 
	run.simple_cal(result.network, result.r_GEP, result.q_GEP)
