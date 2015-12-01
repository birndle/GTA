from collections import defaultdict
from GTAdb import Ptttablesaccid as PTT

filename = '../PSSM/orfg%d/RAxML_distances.orfg%d_raxml'

def get_dist_matrix(gene_num):

	f = open(filename % (gene_num, gene_num), 'r')
	matrix = defaultdict(dict)
	acc_ids = set()
	for line in f:
		words = line.split(' ')
		id1, id2, dist =  (words[0], words[1], float(words[3]))
		acc_ids.add(id1)
		acc_ids.add(id2)
		matrix[id1][id2] = dist
		matrix[id2][id1] = dist
	
	return matrix
		

		