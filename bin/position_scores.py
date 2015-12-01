#!/usr/bin/env python
#File created on May 30 2014

#here the inputs are
#1. input1: profile GTA alignment
#2 input 2: profile viral alignment
#3 input 3: query profile alignment
#calculate (make sure that the columns are conserved)
#1. list of columns that are distinguishing by comparing the viral and GTA alignments
#2 calculate the relative frequency for each of the distinguising positions


#2. a profile alignment of GTA and query
#3. calculate the sequence composition of only GTA sequences in the second alignment,
#4 match with that of distinguishing columns
#5 subset columns from this alignment that are distinguishing
# based on the column that are distinguishing, calculate the score for the query
#4. extract that column number
#5. calculate the score for the query sequence


__author__ = "Migun Shakya"
__email__ = "microbeatic@gmail.com"
__version__ = "0.1"
__license__ = "The MIT License (MIT)"


#to print matplotlib plots in image
import matplotlib
matplotlib.use('Agg')
#--- standard library imports
import argparse

#--- third-party imports
import rpy2.robjects as robjects
from Bio.Align import AlignInfo
from Bio import AlignIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import brewer2mpl
import matplotlib as mpl  # for fonts

#--- project specific imports


def cmdline_parser():
    """
    creates an argparse instance
    """
    parser = argparse.ArgumentParser(description="""calculate scores of query alignment""")
    parser.add_argument('-i1', '--viral',
                        help='viral profile alignment',
                        required=True)
    parser.add_argument('-i2', '--lc',
                        help='large cluster profile GTA alignment',
                        required=True)
    parser.add_argument('-i3', '--query',
                        help='query profile aligned alignment', required=True)
    parser.add_argument('-i4', '--lc_query',
                        help='lc profile aligned with query alignment', required=True)
    parser.add_argument('-o', '--out_file',
                        help='a csv file with sequence header of query and its score', required=True)
    return parser


def main():
    """
    The main function
    """
    parser = cmdline_parser()
    args = parser.parse_args()

########################################################################################################################
#identify significant columns from viral and lc profile alignnment
########################################################################################################################
    #1 calculate columns that are gappy in both alignment
    gappy_columns = gap_list(args.viral, args.lc, 0.7)

    #2 calculate entropy for the non-gappy columns and extract that are conserved in GTA and not in viruses
    sig_col_list = sig_cols(gappy_columns, args.viral, args.lc, 2, 1.5)

    #3 Further filter the sig_col_list based on abundance of a residue
    #3A calculate dictionary of column position and abundance of each residue
    lc_dict = col_dict(args.lc)
    print "alignment of GTAs ", len(lc_dict)
    #3B subset the dictionary based on sig_col_list
    lc_dict_sig = subset_dict(lc_dict, sig_col_list)
    #3C filter the columns in GTA based on abundance (0.75)
    sig_col_list2 = sig_cols2(args.lc, sig_col_list)
    #3D subset the dictionary based on the new significant column list
    lc_dict_sig2 = subset_dict(lc_dict, sig_col_list2)
    print "significant column list %s" % sig_col_list2
####################################################################################################################
#find column position in the new profile alignment of GTA that corresponds to the original
########################################################################################################################

    #1 create dictionary for reference sequence of GTA that was profile aligned with the query
    lc2_dict = col_dict(args.lc_query)
    print "alignment of queries ", len(lc2_dict)
    #calculate columns that are similar between previous alignment and the new one
    #column_pos1 = comp_dict(lc_dict_sig2, lc2_dict, lc_dict)

    column_pos = comp_dict1(lc2_dict, lc_dict)
    column_pos1 = subset_dict(column_pos, sig_col_list2)
    print "Number of significant columns from reference data = %s" % len(sig_col_list2)
    print "Number of significant columns mapped to query alignment = %s" % len(column_pos1)

########################################################################################################################
#calculate score for query sequences and write in a file
########################################################################################################################
    gi_score_list = score_seq(args.query, sig_col_list2, column_pos1, lc_dict)

    fp = open(args.out_file, 'w')
    fp.write('\n'.join('%s, %s' % x for x in gi_score_list))
    fp.close()


def score_seq(aln_file, sig_col_list, col_map, col_abund_dict):
    """
    get the list of tuples with sequence header and score for each sequence
    """
    alignment = AlignIO.read(aln_file, 'fasta')
    aa_prop_list = []
    no_seqs = float(sum(col_abund_dict.get(sig_col_list[0]).values()))
    score_per_pos_list = []
    ncol = len(sig_col_list)
    for sequences in alignment:
        aa_prop_per_seq = []
        #working on column from each sequence
        for col1, col2 in col_map.iteritems():
            if sequences.seq[col2] in col_abund_dict.get(col1):
                #print col_abund_dict.get(col1).get(sequences.seq[col1])
                aa_abund = float(col_abund_dict.get(col1).get(sequences.seq[col2]))
                aa_prop_per_seq.append(aa_abund/no_seqs)
            else:
                aa_prop_per_seq.append(0.0)
        aa_prop_list.append((sequences.id, sum(aa_prop_per_seq)/ncol))
    return aa_prop_list


def comp_dict1(query_dict, lc_dict):
    """compare two dictionaries
    query_dict: dictionary based on the alignment of query sequence
    lc_dict: dictionary of the whole reference alignment"""

    cor_cols = []

    sorted(lc_dict, key=lc_dict.get)
    sorted(query_dict, key=query_dict.get)
    diff = len(query_dict) - len(lc_dict)
    # print diff
    # print len(lc_dict)
    k = 0
    for i in range(0, len(lc_dict)):
        if lc_dict[i] == query_dict[k]:
            cor_cols.append((i, k))
            k = k+1
        else:
            for j in range(1, diff+1):
                if lc_dict[i] == query_dict[k+j]:
                    cor_cols.append((i, k+j))
                    k = k+j
                    break
    cor_cols_dict = dict(cor_cols)
    return cor_cols_dict


def subset_dict(col_abund_dict, sig_col_list):
    """
    subsets columns abundant dictionary with only columns that are distinguishing between the two
    """
    imp_dict = dict([(i, col_abund_dict[i]) for i in sig_col_list if i in col_abund_dict])
    return imp_dict


def sig_cols2(aln_file, sig_col_list):
    """filters sig_col_list further based on the proportion
    aln_file = alignment File
    sig_col_list = list of significant columns that will be filtered further"""
    sig_col2 = []
    dict_all = col_dict(aln_file)
    #subset the dictionary with only significant positions
    sig_dict = dict([(i, dict_all[i]) for i in sig_col_list if i in dict_all])
    #further filter columns that have atleast one residue with 0.75
    for pos, prop in sig_dict.iteritems():
        residue_counts = prop.values()
        #print residue_counts
        residue_prop_list = [float(counts)/sum(residue_counts) for counts in residue_counts]
        for residue_prop in residue_prop_list:
            if residue_prop > 0.75:
                sig_col2.append(pos)
    return sig_col2


def col_dict(aln_file):
    """create dictionary of each column and number of occurence of residues"""
    dict_all = {}
    from Bio import AlignIO
    from collections import Counter

    alignment = AlignIO.read(aln_file, 'fasta')
    for col in range(0, alignment.get_alignment_length()):
        res = [i for i in alignment[:, col]]
        uni_res = Counter(res)
        uni_res_dict = dict(uni_res)
        dict_all[col] = uni_res_dict
    return dict_all


def sig_cols(gap_list, aln1, aln2, ent1_h, ent2_h):
    """
    output list of columns that conserved conserved in one compared
    to another based on the entropy for only columns that
    are not gappy in both alignments

    column score less than ent1_h and their corresponding columns with score higher than ent2_h
    in ent1_h and ent2_h respectively is reported

    columns that are conserved in ent_list1 and divergent in ent_list2 are retained
    """
    import pandas as pd

    ent_list1 = Rentropy_calc(aln1)
    ent_list2 = Rentropy_calc(aln2)

    #subset columns that do not have gap as consensus from list of entropy
    ent_nogaps1 = [(j, i) for j, i in enumerate(ent_list1) if j not in gap_list]
    ent_nogaps2 = [(j, i) for j, i in enumerate(ent_list2) if j not in gap_list]

    #convert to data frame
    ent_nogaps_df1 = pd.DataFrame.from_records(ent_nogaps1)
    ent_nogaps_df1.columns = ['position', 'score']
    ent_nogaps_df2 = pd.DataFrame.from_records(ent_nogaps2)
    ent_nogaps_df2.columns = ['position', 'score']
    #merge columns based on its position
    merged_df = pd.merge(ent_nogaps_df1, ent_nogaps_df2, on='position')
    merged_df.columns = ['position', 'ent1', 'ent2']

    con_ent1_nc_ent2 = list(merged_df[(merged_df['ent1'] > ent1_h) & (merged_df['ent2'] < ent2_h)]['position'])
    return con_ent1_nc_ent2


def Rentropy_calc(aln_file):
    """
    calculate entropy of given alignment using entropy function
    of bio3d package in R
    returns a list with entropy scores for each position
    """
    import rpy2.robjects as robjects
    r = robjects.r
    r.assign('aln_file', aln_file)
    r('require(bio3d)')
    r('read_fa = bio3d::read.fasta(aln_file, rm.dup=FALSE)')
    r('fa_entropy = entropy(read_fa)')
    fa_entropy = r('fa_entropy$H')
    Rentropy = list(fa_entropy)
    return Rentropy


def gap_list(aln_file1, aln_file2, thre):
    """outputs list of column positions that have gap as consensus (thre) in both aln_file1
    and aln_file2 alignments"""
    viral = aln_file1
    lc = aln_file2
    viral_summ = AlignInfo.SummaryInfo(viral)
    lc_summ = AlignInfo.SummaryInfo(lc)
    viral_con = viral_summ.gap_consensus(threshold=thre)
    lc_con = lc_summ.gap_consensus(threshold=thre)
    column_list = []
    for columns in range(0, len(viral_con)):
        if viral_con[columns] == lc_con[columns] == '-':
            column_list.append(columns)
    return column_list


if __name__ == '__main__':
    main()
