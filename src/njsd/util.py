import networkx as nx
import numpy as np


def parse_network(network_fp):
    """Parses network file and returns a network instance and a gene set.
    Attribute:
        network_fp (str): File path to a network file.
    """
    graph = nx.Graph()
    gene_set = set()

    with open(network_fp) as inFile:
        inFile.readline()  # Skip header.
        for line in inFile.readlines():
            gene1, gene2 = line.strip().split()

            graph.add_edge(gene1, gene2)
            gene_set.add(gene1)
            gene_set.add(gene2)

    return graph, gene_set


def parse_gene_set(gene_set_fp):
    """Parses gene set file and returns a (group)-(gene set) dictionary.
    Attribute:
        gene_set_fp (str): File path to a gene set file.
    """
    group_gene_set_dict = dict()

    with open(gene_set_fp) as inFile:
        for line in inFile.readlines():
            tokens = line.strip().split('\t')
            group = tokens[0]
            gene_set = set(tokens[1:])

            group_gene_set_dict[group] = gene_set

    return group_gene_set_dict


def parse_gene_expression(gene_expression_fp, mean=False):
    """Parses gene expression file and returns a (gene identifier)-(expression) dictionary.
    Attribute:
        gene_expression_fp (str): File path to a gene expression file.
        mean (bool): When making a normal(reference) gene expression profile, you might use
            average values of gene expressions for each gene. In this case, pass mean=True.
    """
    gene_expression_dict = dict()
    with open(gene_expression_fp) as inFile:
        inFile.readline()  # Skip header.
        for line in inFile.readlines():
            tokens = line.strip().split('\t')
            gene_identifier = tokens[0]

            if mean:
                expression = np.log2(np.mean([float(t) for t in tokens[1:]]) + 1.0)
            else:
                expression = np.log2(float(tokens[1]) + 1.0)

            gene_expression_dict[gene_identifier] = expression

    return gene_expression_dict


def get_maximally_ambiguous_network(query_gene_expression_dict):
    """Return maximally ambiguous network by assigning expression value of 1.0 to all of the genes.
    """
    return {gene: 1.0 for gene in query_gene_expression_dict.keys()}
