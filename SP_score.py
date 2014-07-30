__author__ = "Dan Birnbaum"
__email__ = "danpbirnbaum@gmail.com"
__version__ = "0.1"


#--- standard library imports
import argparse
import sys
import os

#--- third-party imports
from Bio import AlignIO
import SubstitutionMatrix
#from Thompson_score import do_seqlogo, extract_sig_cols


def sim_score(col1, col2, subst):
    """
    computes similarity between corresponding columns of two alignments
    """
    score = 0
    for a in col1:
        for b in col2:
            score += subst[(a, b)]
    return score


def main():
    """
    The main function
    """
    parser = cmdline_parser()
    args = parser.parse_args()
    t = args.threshold

    gta_align = AlignIO.read(args.gta, "fasta")
    viral_align = AlignIO.read(args.viral, "fasta") 

    scores = sim_score(gta_align, viral_align)

    # extract the position indicies for columns with scores exceeding the user-defined threshold
    sig_cols = [i for i, j in enumerate(scores) if j > np.percentile(scores, t)]

    # parse alignments by extracting only signature columns
    
    # parsed_gta = extract_sig_cols(aln_gta, sig_cols)
    # parsed_viral = extract_sig_cols(aln_viral, sig_cols)

    #parsed_gta_fasta = AlignIO.write(parsed_gta, "test-gta.fasta", "fasta")
    #parsed_viral_fasta = AlignIO.write(parsed_viral, "test-viral.fasta", "fasta")

    #do_seqlogo("test-gta.fasta", args.output2, sig_cols)
    #do_seqlogo("test-viral.fasta", args.output1, sig_cols)



if __name__ == 'main':
    main()
