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


	def select_features(self, X, y, n=500, method='BNS'):
		if method == 'unan':
			zipped = zip(*X)
			y = np.array(y)
			scores = []
			for col in zipped:
				col = np.array(list(col))
				pos = col[y==1]
				tp = pos[pos>0]
				neg = col[y==0]
				fp = neg[neg>0]
				if float(len(tp))/len(pos) > 0.2:
					print "go GTAs!"
				if float(len(fp))/len(neg) > 0.2:
					print "go Viruses"
				print ">"

				# print "Percentage of GTAs containing kmer: %f" % (float(len(tp))/len(pos))
				# print "Percentage of Viruses containing kmer: %f" % (float(len(fp))/len(neg))

		if method == 'BNS':
			zipped = zip(*X)
			y = np.array(y)
			scores = []
			for col in zipped:
				col = np.array(list(col))
				thr = np.average(col)
				pos = col[y==1]
				tp = pos[pos>thr]
				neg = col[y==0]
				fp = neg[neg>thr]

				prev_pos = float(len(tp))/len(pos)
				prev_neg = float(len(fp))/len(neg)

				bns_score = abs(norm.ppf(prev_pos)-norm.ppf(prev_neg))
				scores.append(bns_score)
			
			n_best_feats_idx = np.array(scores).argsort()[-n:]
			n_best_feats_idx = n_best_feats_idx[::-1]
			new_X = []
			for x in X:
				new_x = list(x[n_best_feats_idx])
				new_X.append(new_x)
			return np.array(new_X), n_best_feats_idx

		else:
			selector = SelectKBest(score_func=chi2, k=n)
			new_X = selector.fit_transform(X,y)
			return new_X





