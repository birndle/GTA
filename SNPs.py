import pandas as pd
from pandas import *

listA = []
listB = []
listC = []
if __name__ == '__main__':
	df = pd.read_table('reordered_snps.txt')
	
	def regex_fcn(row):		
		import re

		# first 3 sequences conserved, second 3 conserved but different
		def critA(s): 
			regex = r'(\w)\1{2}(?!\1)(\w)\2{2}'
			if re.search(regex, s):
				return True
			else:
				return False

		# first 3 conserved, second 3 not conserved
		def critB(s):
			regex = r'(\w)\1{2}(\w)(\w)(?:(?!\2)|(?!\3))\w'
			if re.search(regex, s):
				return True
			else:
				return False

		# second 3 conserved, first 3 not conserved
		def critC(s):
			regex= r'(\w)(\w)(?:(?!\1)|(?!\2))\w(\w)\3{2}'

		seq = row['SNP pattern']
		pos = tuple([row[i] for i in range(1,7)])

		if critA(seq):
			listA.append(pos)
		elif critB(seq):
			listB.append(pos)
		elif critC(seq):
			listC.append(pos)

	df.apply(regex_fcn, axis=1)
