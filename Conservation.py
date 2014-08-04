import math

# Shannon Entropy of a string
# = minimum average number of bits per symbol
# required for encoding the string
def shannon(st, stereo=False):

	# st = '00010101011110' # Shannon entropy for this would be 1 bit/symbol

	stList = list(st)
	alphabet = list(set(stList)) # list of symbols in the string
	partitions = [[a] for a in alphabet]
	
	if stereo:
		partitions = []
		partitions.append(['A', 'V', 'L', 'I', 'M', 'C']) # aliphatic
		partitions.append(['F', 'W', 'Y', 'H']) # aromatic
		partitions.append(['S', 'T', 'N', 'Q']) # polar
		partitions.append(['K', 'R']) # positive
		partitions.append(['D', 'E']) # negative
		partitions.append(['G', 'P']) # special conformations

	# calculate the frequency of each symbol in the string
	freqList = []
	for part in partitions:
	    ctr = 0
	    for sym in stList:
	        if sym in part:
	            ctr += 1
	    freqList.append(float(ctr) / len(stList))
	# Shannon entropy
	ent = 0.0
	for freq in freqList:
		if freq > 0:
			ent = ent + freq * math.log(freq, 2)
	ent = -1*ent
	return ent
