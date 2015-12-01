from collections import defaultdict

class Sequence:

	def __init__(self, seq, w=1):
		self.seq = seq
		self.weight = w
		self.counts = None

	def get_kmers(self, k, var=False):
		bag = defaultdict(int)
		if var == True:
			K = range(1, k+1)
		else:
			K = [k]
		for k in K:
			for x in range(len(self.seq)-k+1):
				kgram = self.seq[x:x+k]
				bag[kgram] += self.weight
		self.counts = bag
		return bag