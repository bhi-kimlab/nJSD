import numpy

def simple_cal(W_NETWORK, REF_EXP, Q_EXP):
	import Calculate_entropy as CE
	import Share_part as SHARE
	network_name = W_NETWORK
	G, whole_gene_set = SHARE.load_string(network_name)
	## Build Network

	q_dic = SHARE.normal_express(REF_EXP)	
	edgevalue_dic = SHARE.personal_gene_info_load(Q_EXP)
	## GEPs load


	module_geneset_list = [edgevalue_dic.keys()]

	
	outstream = Q_EXP
	for m in module_geneset_list:
		entropy_dic = CE.cal_JSD(G, edgevalue_dic, q_dic, m )
		total_entropy = numpy.average(entropy_dic.values())
		## nJSD from Ref to Query GEP.

		i_dic = SHARE.make_Istate(edgevalue_dic)
		i_entropy_dic = CE.cal_JSD(G, i_dic, q_dic, m )
		i_total_entropy = numpy.average(i_entropy_dic.values())
		## nJSD from Query GEP to state H that all genes have equal expression level.

		tITH = (total_entropy) / (i_total_entropy + total_entropy)

		res = "[Ref. -> Query: %f]\t[Query -> stateH: %f]\t<tITH: %f>" % (total_entropy, i_total_entropy, tITH)
		outstream += "\t" + res
		outstream += "\n"
	print outstream


