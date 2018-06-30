import numpy
import networkx as Nx

def cal_JSD(GRAPH, QUERY_DIC, NORMAL_DIC, MODULE_LIST):
    entropy_dic = {}
    nkeys       = NORMAL_DIC.keys()
    assert nkeys != 0, "GEP has 0 nodes"

    for node in MODULE_LIST:
        '''Find neighbor in GRAPH.
        Based on netighbor information, calculate JSD of node
        
        Neighbor which has no expression value
            in both reference GEP and query GEP will be ignored.

        KL : Kullback Liebler Divergence
        '''

        try:
            assert Nx.info(GRAPH, node), "Unexpected termination:              \
                                            node is not exist in GRAPH"

            node_neighbors = [n for n in GRAPH.neighbors(node)]

            if len(node_neighbors) == 0:
                continue
        except:
            print "Unexpected termination: Graph_neighbor error"
            exit(0)

        query_score  = []
        normal_score = []
        for neighbor in node_neighbors:
            try:
                q = QUERY_DIC[neighbor]
                n = NORMAL_DIC[neighbor]
            except KeyError:
                continue
                
            query_score.append(q)
            normal_score.append(n)
        
        assert len(query_score) == len(normal_score)
        total_neighbors = len(query_score)

        total_q = numpy.sum(query_score)
        total_n = numpy.sum(normal_score)
        if (total_q == 0) and (total_n == 0) :
            continue

        ''' convert gene expression level into probability '''
        tp_l = []
        np_l = []
        mp_l = []
        for ei in range(total_neighbors):
            if(total_q == 0):
                tp = 0
            else:
                tp = float(query_score[ei]) / total_q

            if(total_n == 0):
                np = 0
            else:
                np = float(normal_score[ei]) / total_n

            if(total_n + total_q) == 0 :
                mp = 0
            else:
                mp = float((tp+np) / (2))
            
            tp_l.append(tp)
            np_l.append(np)
            mp_l.append(mp)

        KL_tm = 0.0
        KL_nm = 0.0
        for ni in range(total_neighbors):
            if(mp_l[ni] == 0 ):
                continue

            if(tp_l[ni] != 0 ):
                KL_tm = KL_tm + tp_l[ni]*numpy.log(tp_l[ni]/mp_l[ni])

            if(np_l[ni] != 0):
                KL_nm = KL_nm + np_l[ni]*numpy.log(np_l[ni]/mp_l[ni])
        
        jsd = KL_tm + KL_nm

        entropy_dic[node] = jsd/2
    return entropy_dic

def cal_whole_JSD(FIRST_DIC, SECOND_DIC, MODULE_LIST):
    '''It doesn't use network information
    just calculate JSD with whole transcriptome
    '''

    first_d  = []
    second_d = []
    for each_gene in MODULE_LIST:
        try:
            first_d.append(FIRST_DIC[each_gene]) 
            second_d.append(SECOND_DIC[each_gene]) 
        except KeyError:
            continue

    total_f = sum(first_d)
    total_s = sum(second_d)

    fp_l = []
    sp_l = []
    mp_l = []
    for each_gene in MODULE_LIST:
        try:
            fp = float(FIRST_DIC[each_gene]) / total_f
            sp = float(SECOND_DIC[each_gene]) / total_s
        except KeyError:
            continue

        if(total_f + total_s ) == 0 :
            mp = 0
        else:
            mp = float( (fp+sp) / (2) )
        
        fp_l.append(fp)
        sp_l.append(sp)
        mp_l.append(mp)

    KL_fm = 0.0
    KL_sm = 0.0
    for ni in range(len(mp_l)):
        if(mp_l[ni] == 0 ):
            continue

        if(fp_l[ni] != 0 ):
            KL_fm = KL_fm + fp_l[ni]*numpy.log(fp_l[ni]/mp_l[ni])

        if(sp_l[ni] != 0):
            KL_sm = KL_sm + sp_l[ni]*numpy.log(sp_l[ni]/mp_l[ni])
    
    jsd = KL_fm + KL_sm

    return (jsd/2)
