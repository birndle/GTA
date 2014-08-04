import Sequence
import mmh3 

from sklearn.feature_selection import SelectKBest, chi2
import numpy as np
from scipy.stats import norm as norm

class Vectorizer:

	def __init__(self, k, num_bins, var=True, hashed=False):
		self.k = k
		self.b = num_bins
		self.var = var
		self.bag = None	 # a list of all unique kmers appearing in the training data
		self.hashed = hashed

	def hash_features(self, bag, weight=1):
		bins = [0]*self.b
		for word in bag:
			hash_key = mmh3.hash(word) % self.b
			bins[hash_key] += bag[word]*weight
		return bins


	def vectorize(self, seq):
		seq.get_kmers(self.k, self.var)
		feat_vec = []
		corres_kmers = []
		if self.hashed:
			feat_vec = self.hash_features(bag)
		else:
			for word in self.bag:
				feat_vec.append(seq.counts[word])
				corres_kmers.append(word)
		return feat_vec, corres_kmers


	def set_bag(self, bag):
		self.bag = bag
		

	def get_sparse_features(self, train_set, thr):
		zipped = zip(*train_set) 
		to_prune = [i for i in range(len(train_set[0])) if sum(zipped[i]) < thr]
		to_prune.sort(reverse=True)
		return to_prune


	def select_features(self, X, y, n=500, num_classes=2):
		
		def spread(func):
			ranked_feats = [func(i) for i in range(num_classes)]
			best_feats_idx = []
			i = 0
			
			while len(best_feats_idx) < n:
				rankings = list(ranked_feats[i%(num_classes+1)])
				fx = rankings.pop(0)
				while fx in best_feats_idx:
					fx = rankings.pop(0)
				best_feats_idx.append(fx)

			return best_feats_idx

		def BNS(clas):
			zipped = zip(*X)
			Y = np.array(y)
			scores = []
			pos_scores = []
			neg_scores = []
			
			for col in zipped:
				col = np.array(list(col))
				thr = np.average(col)
				
				pos = col[Y==clas]
				tp = pos[pos>thr]				
				neg = col[Y!=clas]
				fp = neg[neg>thr]


				prev_pos = float(len(tp))/len(pos)
				prev_neg = float(len(fp))/len(neg)

				bns_score = abs(norm.ppf(prev_pos) - norm.ppf(prev_neg))
				scores.append(bns_score)

				if prev_pos > prev_neg:
					pos_scores.append(bns_score)
					neg_scores.append('nope')
				else:
					pos_scores.append(bns_score)
					neg_scores.append('nope')


			# scores contains list of BNS scores for every feature

			# get the indicies of the n best features, and sort them from most best to least best
			n_best_feats_idx = np.array(scores).argsort()[-n:]
			n_best_feats_idx = n_best_feats_idx[::-1]

			# extract the indicies corresponding to GTA-indicative features and Virus-indicative features
			pos_feats_idx = [i for i in n_best_feats_idx if type(pos_scores[i]) is not str]
			neg_feats_idx = [i for i in n_best_feats_idx if type(neg_scores[i]) is not str]
			
			return n_best_feats_idx

		if num_classes > 2:
			return spread(BNS)
		else:
			return BNS(0)

		# elif method == 'chi':
		# 	selector = SelectKBest(score_func=chi2, k=n)
		# 	new_X = selector.fit_transform(X,y)
		# 	return new_X





