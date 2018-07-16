import numpy
import datetime
import Calculate_entropy as CE
import Share_part as SHARE


def gene_set_cal(W_NETWORK, REF_EXP, QRY_EXP, GENESET_LIST):
    ''' BUILD Network '''
    network_name      = W_NETWORK
    G, whole_gene_set = SHARE.load_string(network_name)

    ''' Gene expression level loading'''
    normal_gep_dic = SHARE.normal_gene_info_set_load(REF_EXP)    
    query_gep_dic  = SHARE.personal_gene_info_load(QRY_EXP)

    ''' Gene set list loading '''
    geneset_dic = SHARE.load_geneset(GENESET_LIST)

    ''' find module in built network '''
    module_geneset_list = query_gep_dic.keys()

    dt = datetime.datetime.now()
    ''' Cal. '''
    geneset_names = geneset_dic.keys()
    for geneset_name in geneset_names:
        geneset_module = geneset_dic[geneset_name]
        entropy_dic, dump_gene_set = CE.cal_JSD(G, query_gep_dic, normal_gep_dic, geneset_module)
        
        ''' dump '''
        if len(dump_gene_set) != 0:
            try:
                DUMP_FD = open("dumpgene_"+dt.strftime("%S%M%H%d%m%Y"), "w")
                print >> DUMP_FD, geneset_name, " ".join(list(dump_gene_set))
            except:
                DUMP_FD = open("dumpgene_"+dt.strftime("%S%M%H%d%m%Y"), "w")
                
        ''' nJSD from Ref to Query gene expression profile '''
        total_entropy = numpy.average(entropy_dic.values())

        i_dic           = SHARE.make_Istate(query_gep_dic)
        i_entropy_dic   = CE.cal_JSD(G, i_dic, query_gep_dic,geneset_module)[0]
        i_total_entropy = numpy.average(i_entropy_dic.values())
        
        ''' nJSD from <query> gene expression profile to state H 
                                 that all genes have equal expression level '''
        tITH = (total_entropy) / (i_total_entropy + total_entropy)

        res  = "[Ref. -> Query: %f]\t[Query -> stateH: %f]\t<tITH: %f>" \
                                       % (total_entropy, i_total_entropy, tITH)
        print "%s \t %s \t %s " % (QRY_EXP, geneset_name, res)
    try:
        DUMP_FD.close()
    except:
        ''' padding '''


def whole_cal(W_NETWORK, REF_EXP, QRY_EXP):
    ''' BUILD Network '''
    network_name      = W_NETWORK
    G, whole_gene_set = SHARE.load_string(network_name)

    ''' Gene expression level loading'''
    normal_gep_dic = SHARE.normal_gene_info_set_load(REF_EXP)    
    query_gep_dic  = SHARE.personal_gene_info_load(QRY_EXP)

    ''' find module in built network '''
    module_geneset_list = query_gep_dic.keys()

    ''' Cal. '''
    outstream = QRY_EXP
    
    dt = datetime.datetime.now()

    entropy_dic, dump_gene_set = CE.cal_JSD(G, query_gep_dic, normal_gep_dic, module_geneset_list)

    ''' dump '''
    if len(dump_gene_set) != 0:
        DUMP_FD = open("dumpgene_"+dt.strftime("%S%M%H%d%m%Y"), "w")
        print >> DUMP_FD, " ".join(list(dump_gene_set))
        DUMP_FD.close() 

    ''' nJSD from Ref to Query gene expression profile '''
    total_entropy = numpy.average(entropy_dic.values())

    i_dic           = SHARE.make_Istate(query_gep_dic)
    i_entropy_dic   = CE.cal_JSD(G, i_dic, query_gep_dic, module_geneset_list)[0]
    i_total_entropy = numpy.average(i_entropy_dic.values())
    
    ''' nJSD from <query> gene expression profile to state H 
                             that all genes have equal expression level '''
    tITH = (total_entropy) / (i_total_entropy + total_entropy)

    res  = "[Ref. -> Query: %f]\t[Query -> stateH: %f]\t<tITH: %f>" \
                                   % (total_entropy, i_total_entropy, tITH)
    outstream += "\t" + res + "\n"

    print outstream
