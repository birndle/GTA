from Bio.Align.Applications import MuscleCommandline as musc
from Bio import AlignIO
from Bio import SeqIO
from StringIO import StringIO
from Bio.Align import MultipleSeqAlignment as msa
import argparse

"""
Input: filename for file containing unaligned sequences to be aligned via MUSCLE
Returns: newly-created multiple sequence alignment object 
"""

def muscle_aln(seqs_filename, output):
	cmd = musc(input=seqs_filename, out=output)
	cmd()
	aln = AlignIO.read(output, "fasta")
	return aln


def label_seqs(aln, label, output):
	for record in aln:
		record.id = label
	SeqIO.write(aln, output, 'fasta')


"""
Input: filenames for two alignments to be profile aligned
Returns: profile aligned MSA object 
"""
def profile_aln(aln1_filename, aln2_filename, output):
	cmd = musc(in1=aln1_filename, in2=aln2_filename, out=output, profile=True)
	cmd()
	prof_aln = AlignIO.read(output, "fasta")
	return prof_aln


"""
Input: 

prof_aln = MSA object, contains profile alignment 
label1 = contains id name of records to extract from profile alignment, i.e. 'gta'
label2 = contains id name of other records to extract from profile alignment, i.e. 'viral'

Returns: tuple containing the two extracted MSA objects, i.e. (gta, viral)
"""

def split_aln(prof_aln, label1, label2):
	records1 = []
	records2 = []
	for record in prof_aln:
		if record.id == label1:
			records1.append(record)
		elif record.id == label2:
			records2.append(record)
	msa1 = msa(records1)
	msa2 = msa(records2)
	return msa1, msa2


def cmdline_parser():
        """
        creates an argparse instance
        """
        parser = argparse.ArgumentParser(description=""" """)
        parser.add_argument("-g", "--gta",
                            help="""gta sequences""",
                            dest="gta",
                            required=True)
        parser.add_argument("-v", "--viral",
                            help="""viral sequences""",
                            dest="viral",
                            required=True)
        parser.add_argument("-q", "--query",
                            help="""query sequence""",
                            dest="query",
                            required=True)
        return parser


def main():
    parser = cmdline_parser()
    args = parser.parse_args()
    buff_gta = 'temp_gta.aln'
    buff_vir = 'temp_vir.aln'

    # create alignments for reference sets
    aln_gta = muscle_aln(args.gta, buff_gta)
    aln_vir = muscle_aln(args.viral, buff_vir)

    # label records in reference MSAs as 'gta' or 'viral'
	label_seqs(aln_gta, 'gta', buff_gta)
	label_seqs(aln_vir, 'viral', buff_vir)

	# profile align reference GTA with reference Virus
    prof_aln = profile_aln(buff_gta, buff_vir, 'temp_prof.pln')

    # extract GTA and Virus reference sets out of the profile alignment
    ref_gta, ref_vir = split_aln(prof_aln, 'gta', 'viral')

    buff_query = 'temp_query.aln'
    
    # align sequences in query set
    aln_query = muscle_aln(args.query, buff_query)

    # label sequences in query MSA as 'query'
    label_seqs(aln_query, 'query', buff_query)

    # profile align query MSA with reference GTA 
    prof_aln = profile_aln(buff_gta, buff_query, 'temp_prof.pln')

    # extract query MSA that has been profile aligned with reference GTA
    ref_gta_prof_query, query_prof_ref_gta = split_aln(prof_aln, 'gta', 'query')



if __name__ == '__main__':
    main()






