__version__ = '0.1.3'

import logging

import njsd.util as util
import njsd.entropy as entropy


logger = logging.getLogger('nJSD')
formatter = logging.Formatter('[%(name)s %(levelname)s %(asctime)s] %(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


def njsd_all(network, ref, query, file, verbose=True):
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


def njsd_geneset(network, ref, query, gene_set, file, verbose=True):
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
    gene_set_present = set(query_gene_expression_dict.keys())

    with open(file, 'w') as outFile:
        print('Gene_set_ID', 'nJSD_NT', 'nJSD_TA', 'tITH', sep='\t', file=outFile)

    for group, gene_set in group_gene_set_dict.items():
        gene_set_to_be_analyzed = gene_set.intersection(gene_set_present)
        # If no genes are available for the group, just ignore it.
        if len(gene_set_to_be_analyzed) == 0:
            logger.warning('%s has no genes available for analysis. Ignoring the group.' % group)
            continue

        # If every gene has a single neighbor, just ignore it.
        if all([graph.degree(gene) == 1 for gene in gene_set_to_be_analyzed]):
            logger.warning('%s has no genes with enough neighbors. Ignoring the group.' % group)
            continue

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
