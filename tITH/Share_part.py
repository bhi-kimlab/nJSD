import numpy

def load_string(WHICH_NETWORK):
    import networkx as Nx

    GRAPH = Nx.Graph()

    FD = open(WHICH_NETWORK,"r")
    FD.readline()

    gene_set = set()
    for line in FD.readlines():
        block = line.rstrip("\n").split()
        gene1 = block[0]
        gene2 = block[1]

        GRAPH.add_edge(gene1,gene2)

        gene_set.add(gene1)
        gene_set.add(gene2)

    FD.close()
    return GRAPH, gene_set 

def load_kegg(FILENAME):
    FD = open(FILENAME,"r")

    module_set_dic = {}
    for line in FD.readlines():
        block = line.rstrip("\n").split("\t")
        
        module_set_dic[block[0]] = set(block[2:])

    FD.close()
    return module_set_dic

def personal_gene_info_load(FILENAME):
    FD = open(FILENAME, "r")
    FD.readline()

    edge_dic = {}
    for line in FD.readlines(): 
        block = line.rstrip("\n").split()
        gn    = block[0]
        
        edge_dic[ gn ] = numpy.log2(float(block[1]) + 1.0)

    FD.close()
    return edge_dic

def normal_gene_info_set_load(FILENAME):
    FD = open(FILENAME,"r")
    FD.readline()

    ndic = {}
    for line in FD.readlines():
        block = line.rstrip("\n").split()
        nid   = block[0]
        
        ndic[nid] = numpy.log2(numpy.average(float(block[1])) + 1.0)

    FD.close()
    return ndic

def make_Istate(original_dic):
    new_dic = {}
    keys    = original_dic.keys()
    for k in keys:
        new_dic[k] = 1.0

    return new_dic
