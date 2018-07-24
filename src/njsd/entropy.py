import numpy as np


def find_neighbors(graph, node):
    """Find neighbors of a node in a graph.
    Attribute:
        graph (networkx.Graph): A graph instance.
        node (str): Name of a node.
    """
    return list(graph.neighbors(node))


def get_neighbor_expression_vector(neighbors, gene_expression_dict):
    """Get an expression vector of neighboring genes.
    Attribute:
        neighbors (list): List of gene identifiers of neighboring genes.
        gene_expression_dict (dict): (Gene identifier)-(gene expression) dictionary.
    """
    expressions = []  # Expression vector.
    for gene in neighbors:
        try:
            expression = gene_expression_dict[gene]
        except KeyError:
            continue

        expressions.append(expression)
    return expressions


def exp2prob(expression_vector):
    """Convert an expression vector into a probability vector.
    Attribute:
        expression_vector (list): List of expression values.
    """
    v = np.asarray(expression_vector)
    if np.sum(v) == 0:
        return np.zeros(len(expression_vector))
    else:
        return v / np.sum(v)


def kld(p1, p2):
    """Compute Kullback-Leibler divergence between p1 and p2.
    It assumes that p1 and p2 are already normalized that each of them sums to 1.
    """
    return np.sum(np.where(p1 != 0, p1 * np.log(p1 / p2), 0))


def jsd(p1, p2):
    """Compute Jensen-Shannon divergence between p1 and p2.
    It assumes that p1 and p2 are already normalized that each of them sums to 1.
    """
    m = (p1 + p2) / 2
    return (kld(p1, m) + kld(p2, m)) / 2


def njsd(network, ref_gene_expression_dict, query_gene_expression_dict, gene_set):
    """Calculate Jensen-Shannon divergence between query and reference gene expression profile.
    """
    gene_jsd_dict = dict()
    reference_genes = ref_gene_expression_dict.keys()
    assert len(reference_genes) != 'Reference gene expression profile should have > 0 genes.'

    for gene in gene_set:
        if gene not in network.nodes:
            continue

        neighbors = find_neighbors(network, gene)

        query_expression_vec = get_neighbor_expression_vector(neighbors, query_gene_expression_dict)
        ref_expression_vec = get_neighbor_expression_vector(neighbors, ref_gene_expression_dict)
        assert len(query_expression_vec) == len(ref_expression_vec), 'Topology of reference network and query network differs. Please check.'

        # A gene which has non-expressed neighbors is ignored.
        if np.sum(query_expression_vec) == 0 and np.sum(ref_expression_vec) == 0:
            continue

        query_p_vec = exp2prob(query_expression_vec)
        ref_p_vec = exp2prob(ref_expression_vec)

        gene_jsd_dict[gene] = jsd(query_p_vec, ref_p_vec)

    return np.mean(list(gene_jsd_dict.values()))
