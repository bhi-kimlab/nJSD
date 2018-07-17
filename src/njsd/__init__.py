__version__ = '0.1.0'

import numpy as np

import util
import entropy

def njsd_all(network, ref, query):
    """Compute transcriptome-wide nJSD between reference and query expression profiles.
    Attribute:
        network (str): File path to a network file.
        ref (str): File path to a reference expression file.
        query (str): File path to a query expression file.
    """
    graph, gene_set_total = util.parse_network(network)
    ref_gene_expression_dict = util.parse_gene_expression(ref, mean=True)
    query_gene_expression_dict = util.parse_gene_expression(ref, mean=False)
    maximally_ambiguous_gene_experession_dict = util.get_maximally_ambiguous_network(query_gene_expression_dict)
    gene_set_present = set(query_gene_expression_dict.keys())

    normal_to_tumor_njsd = entropy.njsd(network=graph,
                                     ref_gene_expression_dict=ref_gene_expression_dict,
                                     query_gene_expression_dict=query_gene_expression_dict,
                                     gene_set_present=gene_set_present)

    tumor_to_ambiguous_njsd = entropy.njsd(network=g,
                                        ref_gene_expression_dict=maximally_ambiguous_gene_experession_dict,
                                        query_gene_expression_dict=query_gene_expression_dict,
                                        gene_set_present=gene_set_present)

    return normal_to_tumor_njsd / (normal_to_tumor_njsd + tumor_to_ambiguous_njsd)


if __name__ == '__main__':
    main()
