"""
INPUT: Name of file containing data table
RETURNS: tuple containing 3 lists, e.g. (list1, list2, list3)

Each list contains tuples of 7 values. Each of these tuples corresponds to a different conserved column.
The 7 values contained in the tuples are the positions within the respective, unaligned genomes 
corresponding to the conserved column of the alignment.

list1: first 4 conserved AND second 3 conserved, BUT different
list2: first 4 conserved, second 3 NOT conserved
list3: first 4 NOT conserved, second 3 conserved
"""

def get_conserved_positions(table_file):
	import pandas as pd
	from pandas import *

	listA = []
	listB = []
	listC = []

	df = pd.read_table(table_file)
	
	def regex_fcn(row):		
		import re

		# first 4 sequences conserved, second 3 conserved but different
		def critA(s): 
			regex = r'(\w)\1{3}(?!\1)(\w)\2{2}'
			if re.search(regex, s):
				return True
			else:
				return False

		# first 4 conserved, second 3 not conserved
		def critB(s):
			regex = r'(\w)\1{3}(\w)(\w)(?:(?!\2)|(?!\3))\w'
			if re.search(regex, s):
				return True
			else:
				return False

		# first 4 not conserved, second 3 conserved
		def critC(s):
			regex= r'(\w)(\w)(\w)(?:(?!\1)|(?!\2)|(?!\3))\w(\w)\4{2}'
			if re.search(regex, s):
				return True
			else:
				return False

		seq = row['SNP pattern']
		pos = tuple([row[i] for i in range(1,8)])

		if critA(seq):
			listA.append(pos)
		elif critB(seq):
			listB.append(pos)
		elif critC(seq):
			listC.append(pos)

	df.apply(regex_fcn, axis=1)
	
	return (listA, listB, listC)


