import numpy

def simple_cal(W_NETWORK, REF_EXP, QRY_EXP):
    import Calculate_entropy as CE
    import Share_part as SHARE

    ''' BUILD Network '''
    network_name      = W_NETWORK
    G, whole_gene_set = SHARE.load_string(network_name)

    ''' Gene expression level loading'''
    normal_gep_dic = SHARE.normal_gene_info_set_load(REF_EXP)    
    query_gep_dic  = SHARE.personal_gene_info_load(QRY_EXP)

    ''' find module in built network '''
    module_geneset_list = [query_gep_dic.keys()]

    ''' Cal. '''
    outstream = QRY_EXP
    for m in module_geneset_list:
        entropy_dic = CE.cal_JSD(G, query_gep_dic, normal_gep_dic, m )
        
        ''' nJSD from Ref to Query gene expression profile '''
        total_entropy = numpy.average(entropy_dic.values())

        i_dic           = SHARE.make_Istate(query_gep_dic)
        i_entropy_dic   = CE.cal_JSD(G, i_dic, normal_gep_dic, m )
        i_total_entropy = numpy.average(i_entropy_dic.values())
        
        ''' nJSD from <query> gene expression profile to state H 
                                 that all genes have equal expression level '''
        tITH = (total_entropy) / (i_total_entropy + total_entropy)

        res  = "[Ref. -> Query: %f]\t[Query -> stateH: %f]\t<tITH: %f>" \
                                       % (total_entropy, i_total_entropy, tITH)
        outstream += "\t" + res + "\n"
    print outstream
