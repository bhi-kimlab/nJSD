__version__ = '0.1.3'

import numpy as np

import njsd.util as util
import njsd.entropy as entropy


def njsd_all(network, ref, query, file):
    """Compute transcriptome-wide nJSD between reference and query expression profiles.
    Attribute:
        network (str): File path to a network file.
        ref (str): File path to a reference expression file.
        query (str): File path to a query expression file.
    """
    graph, gene_set_total = util.parse_network(network)
    ref_gene_expression_dict = util.parse_gene_expression(ref, mean=True)
    query_gene_expression_dict = util.parse_gene_expression(query, mean=False)
    maximally_ambiguous_gene_experession_dict = util.get_maximally_ambiguous_network(query_gene_expression_dict)
    gene_set_present = set(query_gene_expression_dict.keys())

    with open(file, 'w') as outFile:
        print('nJSD_NT', 'nJSD_TA', 'tITH', sep='\t', file=outFile)

    normal_to_tumor_njsd = entropy.njsd(network=graph,
                                     ref_gene_expression_dict=ref_gene_expression_dict,
                                     query_gene_expression_dict=query_gene_expression_dict,
                                     gene_set=gene_set_present)

    tumor_to_ambiguous_njsd = entropy.njsd(network=graph,
                                        ref_gene_expression_dict=maximally_ambiguous_gene_experession_dict,
                                        query_gene_expression_dict=query_gene_expression_dict,
                                        gene_set=gene_set_present)
    tITH = normal_to_tumor_njsd / (normal_to_tumor_njsd + tumor_to_ambiguous_njsd)

    with open(file, 'a') as outFile:
        print(normal_to_tumor_njsd, tumor_to_ambiguous_njsd, tITH, sep='\t', file=outFile)
    return normal_to_tumor_njsd / (normal_to_tumor_njsd + tumor_to_ambiguous_njsd)


def njsd_geneset(network, ref, query, gene_set, file):
    """Compute gene set-specified nJSD between reference and query expression profiles.
    Attribute;
        network (str): File path to a network file.
        ref (str): File path to a reference expression file.
        query (str): File path to a query expression file.
        geneset (str): File path to a gene set file.
    """
    graph, gene_set_total = util.parse_network(network)
    ref_gene_expression_dict = util.parse_gene_expression(ref, mean=True)
    query_gene_expression_dict = util.parse_gene_expression(query, mean=False)
    group_gene_set_dict = util.parse_gene_set(gene_set)
    maximally_ambiguous_gene_experession_dict = util.get_maximally_ambiguous_network(query_gene_expression_dict)

    with open(file, 'w') as outFile:
        print('Gene_set_ID', 'nJSD_NT', 'nJSD_TA', 'tITH', sep='\t', file=outFile)

    for group, gene_set in group_gene_set_dict.items():
        normal_to_tumor_njsd = entropy.njsd(network=graph,
                                         ref_gene_expression_dict=ref_gene_expression_dict,
                                         query_gene_expression_dict=query_gene_expression_dict,
                                         gene_set=gene_set)

        tumor_to_ambiguous_njsd = entropy.njsd(network=graph,
                                            ref_gene_expression_dict=maximally_ambiguous_gene_experession_dict,
                                            query_gene_expression_dict=query_gene_expression_dict,
                                            gene_set=gene_set)

        tITH = normal_to_tumor_njsd / (normal_to_tumor_njsd + tumor_to_ambiguous_njsd)
        with open(file, 'a') as outFile:
            print(group, normal_to_tumor_njsd, tumor_to_ambiguous_njsd, tITH, sep='\t', file=outFile)
