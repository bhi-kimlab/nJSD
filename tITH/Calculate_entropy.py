import numpy

def cal_JSD(GRAPH, EDGE_DIC, NORMAL_DIC, MODULE_LIST):
	entropy_dic = {}
	nkeys = NORMAL_DIC.keys()
	for node in MODULE_LIST:
		'''Find neighbor in GRAPH.
		Based on netighbor information, calculate JSD of node
		
		Neighbor which has no expression value in both reference GEP and query GEP
			will be ignored.
		KL : Kullback Liebler Divergence
		'''

		try:
			node_neighbors = GRAPH.neighbors(node)
			if(len(node_neighbors) ==0):
				continue
		except:
			continue

		edge_score = []
		normal_score = []
		for neighbor in node_neighbors:
			try:
				e = EDGE_DIC[ neighbor ]
				n = NORMAL_DIC[ neighbor ]
			except KeyError:
				continue
			edge_score.append((e))
			normal_score.append((n))

		total_e = numpy.sum(edge_score)
		total_n = numpy.sum(normal_score)
		if(total_e ==0) & (total_n ==0 ) :
			continue

		if( len(edge_score) != len(normal_score) ):
			print("diff len of neighbor")
			continue
		total_neighbors = len(edge_score)

		tp_l = []
		np_l = []
		mp_l = []
		for ei in range(total_neighbors):
			if(total_e == 0):
				tp = 0
			else:
				tp = float(edge_score[ei]) / total_e

			if(total_n == 0):
				np = 0
			else:
				np = float(normal_score[ei]) / total_n

			if(total_n + total_e ) == 0 :
				mp = 0
			else:
				mp = float( (tp+np) / (2) )
			
			tp_l.append(tp)
			np_l.append(np)
			mp_l.append(mp)
		## made probability distribution with gene-expression level


		KL_tm = 0.0
		KL_nm = 0.0
		for ni in range(total_neighbors):
			if(mp_l[ni] == 0 ):
				continue

			if(tp_l[ni] != 0 ):
				KL_tm = KL_tm + tp_l[ni]*numpy.log( tp_l[ni]/mp_l[ni] ) #natural log

			if(np_l[ni] != 0):
				KL_nm = KL_nm + np_l[ni]*numpy.log( np_l[ni]/mp_l[ni] ) #natural log
		
		jsd = KL_tm + KL_nm

		entropy_dic[node] = jsd/2
#		print(jsd)
	return entropy_dic
#--#


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
			KL_fm = KL_fm + fp_l[ni]*numpy.log( fp_l[ni]/mp_l[ni] ) #natural log

		if(sp_l[ni] != 0):
			KL_sm = KL_sm + sp_l[ni]*numpy.log( sp_l[ni]/mp_l[ni] ) #natural log
	
	jsd = KL_fm + KL_sm

	
#	print(jsd)
	return (jsd/2)

