from Bio import SeqIO
from GTAdb import Accgitable as AccGi, Acctitletable as AccTitle
import subprocess
import random
import argparse

def cmdline_parser():
        """
        creates an argparse instance
        """
        parser = argparse.ArgumentParser(description=""" """)
        parser.add_argument("-g", "--gta",
                            help="""gta sequences""",
                            dest="gta",
                            required=True)
        return parser
"""
Randomly selects other sequences from genomes of provided genes and writes them to
a new, user-named output file in FASTA format.

INPUT:
gta_files: FASTA file containing genes in your training set
output: name of file where you want to write the random other sequences
	

"""
def pull_others(filename, output):

	def get_accnums(homolog_file):
		homologs = SeqIO.parse(homolog_file, "fasta")
		# for h in homologs:
		# 	print h.name
		gis = [homolog.name.split('|')[1] for homolog in homologs]
		sq = AccGi.select(AccGi.accid).where(AccGi.gi << gis)
		accs = [x.accid for x in sq]
		return accs, gis

	def find_file(accid):
		path = '/data2/mshakya/GTA/data/rsync_complete_genomes'
		proc = subprocess.Popen(['find', path, '-name', '%s.faa' % (accid)], stdout=subprocess.PIPE)
		return proc.stdout.read().strip('\n')

	others = [] # list of sequence record objects

	# get genomes corresponding to genes contained in this portion of your training set
	acc_ids, homolog_gis = get_accnums(filename)
	acc_ids = [a.split('.')[0] for a in acc_ids]
	genome_files = [find_file(a) for a in acc_ids]
	genomes = [SeqIO.parse(genome_file, 'fasta') for genome_file in genome_files]

	for genome in genomes:
		chosen = []
		genome = list(genome)
		choice = random.choice(genome)
		while len(chosen) < 2:
			if choice.name.split('|')[1] not in homolog_gis:
				chosen.append(choice)
			choice = random.choice(genome)
		others.extend(chosen)
	
	output_handle = open(output, 'w')
	SeqIO.write(others, output_handle, 'fasta')
	output_handle.close()

if __name__ == '__main__':
	parser = cmdline_parser()
	args = parser.parse_args()
	pull_others(args.gta, args.gta.split('.')[0]+'_others.fasta')



