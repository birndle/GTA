from GTAdb import Mlpredictions as ML, Ptttablesaccid as PTT, Gihomolog as GTAs
from peewee import SelectQuery, UpdateQuery, InsertQuery

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

from matplotlib import pyplot
import numpy as np


def insert_new_rows(data, homo):
	genes = []
	for row in data:
		d = {'homolog':homo, 'gi':row[0], 'thous_best_bns_selected_feats':row[1], 'confidence_1000_feats':row[2]}
		genes.append(d)
	iq = InsertQuery(ML, rows=genes)
	iq.execute()

if __name__ ==  '__main__':

	""" get cluster sizes and add info to column 'cluster_size' in MLpredictions"""

	# unk_ids = [gene.gi for gene in ML.select().where(ML.cluster_size >> None)] # gis for genes in predictions table without cluster sizes
	# queries = SelectQuery(PTT, PTT.gi, PTT.accid, PTT.position).where(PTT.gi << unk_ids)
	# d = {}

	
	# for query in queries:
	# 	neighbors = SelectQuery(PTT, PTT.gi).where(PTT.accid == query.accid).where(PTT.position.between(query.position-10, query.position+10)).where(PTT.gi != query.gi)
	# 	d[query] = neighbors
	
	# i = 0
	# for query in d:
	# 	gene_id = query.gi
	# 	neighbors = d[query]
	# 	neighbor_ids = [n.gi for n in neighbors]
	# 	gta_cluster_size = GTAs.select().where(GTAs.gi << neighbor_ids).count()
	# 	uq = ML.update(cluster_size=gta_cluster_size).where(ML.gi == gene_id)
	# 	i += uq.execute()
	# 	print i
	

	# data = [(gene.confidence_1000_feats, gene.cluster_size) for gene in ML.select(ML.confidence_1000_feats, ML.cluster_size).where(~(ML.cluster_size >> None))]
	# x = np.array([conf for (conf, clust) in data])
	# y = np.array([clust for (conf, clust) in data])
	# pyplot.scatter(x,y)
	# pyplot.xlabel('SVM Decision Function')
	# pyplot.ylabel('Surrounding GTA cluster size')
	# pyplot.savefig('confidence_cluster_corr')


	