import argparse
import sys
import os
import math
from collections import defaultdict

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment as msa
from Bio.SeqRecord import SeqRecord as sr
from Bio.Seq import Seq
from position_scores import gap_list, comp_dict1
import numpy as np
from SubstitutionMatrix import SubstitutionMatrix
from sp_score import sim_score
from Zvelibil_score import stereo_score


def do_seqlogo(input_fasta, output_seqlogo, col_list):
    """
    create sequence logo of given fasta file and then label it with col_list
    """

    import weblogolib
    import corebio
    from weblogolib.colorscheme import ColorScheme
    from weblogolib.colorscheme import ColorGroup
    #assigning colors
    hydrophobicity = ColorScheme([
        ColorGroup("RKDENQ", "blue"),
        ColorGroup("SGHTAP", "green"),
        ColorGroup("YVMCLFIW",  "black"),
        ColorGroup("O",  "dark orange")
    ])

    #Here i changed this so that the column begins with 1 not 0
    col_list2 = [i+1 for i in col_list]
    fin = open(input_fasta)
    seqs = weblogolib.read_seq_data(fin)
    seqs.alphabet = corebio.seq.Alphabet('ACDEFGHIKLMNPQRSTVWYO')
    data = weblogolib.LogoData.from_seqs(seqs)
    data.composition = 'equiprobable'
    options = weblogolib.LogoOptions()
    options.yaxis_scale = 4.5
    options.yaxis_tic_interval = 2.25
    options.annotate = col_list2
    options.stack_width = 30
    options.stacks_per_line = 35
    options.scale_width = 'No'
    options.color_scheme = hydrophobicity
    format_logo = weblogolib.LogoFormat(data, options)
    fout = open(output_seqlogo, 'w')
    pdf_logo = weblogolib.pdf_formatter(data, format_logo)
    fout.write(pdf_logo)
    fout.close()


def cmdline_parser():
    """
    creates an argparse instance
    """
    parser = argparse.ArgumentParser(description=""" """)
    parser.add_argument("-v", "--viral",
                        help="""viral alignments""",
                        dest="viral",
                        required=True)
    parser.add_argument("-g", "--gta",
                        help="""GTA alignments""",
                        dest="gta",
                        required=True)
    parser.add_argument("-q", "--query",
                        help="""query sequence""",
                        dest="query",
                        required=True)
    parser.add_argument("-l", "--lcquery",
                        help='lc profile aligned with query alignment',
                        dest='lcquery', 
                        required=True)
    parser.add_argument("-o1",
                        dest="output1",
                        help="output image file for viral")
    parser.add_argument("-o2",
                        help="output image file for GTA",
                        dest="output2")
    parser.add_argument("-t",
                        dest="threshold",
                        type=float,
                        required=True,
                        help="threshold percentile value, between 0 and 100")
    return parser

def euclid(cg1, cg2):
    dis = 0.0
    dis = sum(np.power((cg1-cg2), 2))
    dis = math.sqrt(dis)
    return dis


def consensus(col, subst, alphabet):
    N = len(col)
    consensus = np.array([0]*len(alphabet)) # initialize consensus vector
    for res in col:
        # build vector representation for each amino acid in the column
        v_res = np.array([subst[(let, res)] for let in alphabet])
        consensus += v_res    
    consensus /= N
    return consensus


def thompson_score(col1, col2, subst, alphab):
    cons1 = consensus(col1, subst, alphab)
    cons2 = consensus(col2, subst, alphab)
    the_score = euclid(cons1, cons2)
    return the_score


def extract_sig_cols(aln, sig_cols):
    cols = []
    for x in sig_cols:
        cols.append(aln[:, x])

    zipped_cols = zip(*cols)
    parsed_seqs = ["".join(zipped_cols[i]) for i in range(len(zipped_cols))]

    records = [sr(Seq(s)) for s in parsed_seqs]
    wl = msa(records)   # MSA object containting only the signature columns
    return wl


def scoring_distr(aln):
    d = defaultdict(lambda: defaultdict(float))
    N = len(aln[:, 0]) # number of sequences in alignment
    for x in range(0, aln.get_alignment_length()):
        for res in aln[:, x]: 
            d[x][res] += 1
        for res in d[x]:
            d[x][res] /= N
    return d


def score(query, reference, col_map, sig_cols):
    scores = [] 
    for seq in query:
        the_score = 0
        for x in sig_cols:
            q_x = col_map[x]
            the_score += reference[x][seq[q_x]] 
        scores.append(the_score)
    scores = [x/len(sig_cols) for x in scores] 
    return scores

def main():
    """
    The main function
    """
    parser = cmdline_parser()
    args = parser.parse_args()
    t = args.threshold

    matrix =                   [[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,0],
                               [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,0],
                               [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,0],
                               [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,0],
                               [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,0],
                               [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0],
                               [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,0],
                               [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,0],
                               [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0],
                               [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,0],
                               [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,0],
                               [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0],
                               [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,0],
                               [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,0],
                               [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,0],
                               [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0],
                               [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,0],
                               [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,0],
                               [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,0],
                               [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,0],
                               [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

    # for i in range(len(matrix)):
    #     for j in range(len(matrix[i])):
    #         if i == j and not i == len(matrix):
    #             matrix[i][j] = 1

    #get substitution matrix
    blosum62 = SubstitutionMatrix("ARNDCQEGHILKMFPSTWYV-", matrix)

    alphabet = "ARNDCQEGHILKMFPSTWYV-"
    scores = []

    #read in the alignment
    aln_viral = AlignIO.read(open(args.viral, 'r'), "fasta")
    aln_gta = AlignIO.read(open(args.gta, 'r'), "fasta")
    gta_query = AlignIO.read(open(args.lcquery, 'r'), "fasta")
    query = AlignIO.read(open(args.query, 'r'), "fasta")

    for x in range(0, aln_viral.get_alignment_length()):
        gta_col = aln_gta[:, x]
        viral_col = aln_viral[:, x]
        scores.append(sim_score(gta_col, viral_col, blosum62))
    
    # find gappy columns
    gappy_cols = gap_list(args.viral, args.gta, 0.7)

    # extract the position indicies for columns with scores exceeding the user-defined threshold
    sig_cols = [i for i, j in enumerate(scores) if j > np.percentile(scores, t)]

    # remove gappy columns from signature column list
    sig_cols = [x for x in sig_cols if x not in gappy_cols]
    print len(sig_cols)

    # parse alignments by extracting only signature columns
    # pruned_gta = extract_sig_cols(aln_gta, sig_cols)
    # pruned_viral = extract_sig_cols(aln_viral, sig_cols)

    # write merged alignments only containing signature columns to fasta files
    # pruned_gta_fasta = AlignIO.write(pruned_gta, "test-gta.fasta", "fasta")
    # pruned_viral_fasta = AlignIO.write(pruned_viral, "test-viral.fasta", "fasta")

    #do_seqlogo("test-gta.fasta", args.output2, sig_cols)
    #do_seqlogo("test-viral.fasta", args.output1, sig_cols)

    # create table of frequencies for residues in reference gta alignment
    # will use these frequencies to score query sequence
    gta_reference_table = scoring_distr(aln_gta)  # frequency table based on initial gta alignment
    gta_query_reference_table = scoring_distr(gta_query)  # gappy table after aligned with query

    # need a map from signature column indices to corresponding query columns 
    reference_to_query_col_map = comp_dict1(gta_query_reference_table, gta_reference_table)

    query_scores = score(query, gta_reference_table, reference_to_query_col_map, sig_cols)

    print sum(query_scores)/len(query_scores)


if __name__ == '__main__':
    main()
