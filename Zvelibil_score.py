#!/usr/bin/env python
#File created on May 9 2014


__author__ = "Migun Shakya"
__email__ = "microbeatic@gmail.com"
__version__ = "0.1"
__license__ = "The MIT License (MIT)"


#--- standard library imports
import argparse
import sys
import os

#--- third-party imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
import rpy2.robjects as robjects
#--- project specific imports


def cmdline_parser():
    """
    creates an argparse instance
    """
    parser = argparse.ArgumentParser(description=""" """)
    parser.add_argument("-i1",
                        help="""viral alignments""",
                        dest="viral",
                        required=True)
    parser.add_argument("-i2",
                        help="""GTA alignments""",
                        dest="gta",
                        required=True)
    parser.add_argument("-o",
                        dest="output",
                        help="output image file")
    return parser


def main():
    """
    The main function
    """
    parser = cmdline_parser()
    args = parser.parse_args()

    #read in the alignment
    aln_viral = AlignIO.read(open(args.viral, 'r'), "fasta")
    aln_gta = AlignIO.read(open(args.gta, 'r'), "fasta")
    aln_viral_stereo = stereo_score(aln_viral)
    aln_gta_stereo = stereo_score(aln_gta)

    #convert to tuples
    viral_stereo_posi = [(i, j) for i, j in enumerate(aln_viral_stereo)]
    gta_stereo_posi = [(i, j) for i, j in enumerate(aln_viral_stereo)]

    #convert to pandas dataframe
    viral_stereo_posi_df = pd.DataFrame.from_records(viral_stereo_posi)
    gta_stereo_posi_df = pd.DataFrame.from_records(gta_stereo_posi)

    #add column headers
    viral_stereo_posi_df.columns = ['position', 'score']
    gta_stereo_posi_df.columns = ['position', 'score']

    #plot the figures
    fig = plt.figure(figsize=(30, 10))
    plt.bar(viral_stereo_posi_df['position'], viral_stereo_posi_df['score'])
    plt.bar(ta_stereo_posi_df['position'], -gta_stereo_posi_df['score'], color='r')
    plt.axis([0, 700, -1.2, 1.2])
    savefig('args.output')


def stereo_score(alignment):
    """
    function that counts the common properties between residues of each column
    """
    #dictionary with properties for each residue
    dic_prop = {'I': [1, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                'L': [1, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                'V': [1, 0, 1, 0, 0, 1, 0, 0, 0, 0],
                'C': [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                'A': [1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
                'G': [1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
                'M': [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'F': [1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                'Y': [1, 1, 0, 0, 0, 0, 1, 0, 0, 0],
                'W': [1, 1, 0, 0, 0, 0, 1, 0, 0, 0],
                'H': [1, 1, 0, 0, 0, 0, 1, 1, 0, 1],
                'K': [1, 1, 0, 0, 0, 0, 0, 1, 0, 1],
                'R': [0, 1, 0, 0, 0, 0, 0, 1, 0, 1],
                'E': [0, 1, 0, 0, 0, 0, 0, 0, 1, 1],
                'Q': [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                'D': [0, 1, 1, 0, 0, 0, 0, 0, 1, 1],
                'N': [0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
                'S': [0, 1, 1, 0, 1, 0, 0, 0, 0, 0],
                'T': [1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
                'P': [0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
                'B': [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                'Z': [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                'X': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                '-': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]}
    score_list = []
    for i in range(0, alignment.get_alignment_length()):
        #extract the unique residues in the alignment
        column = ''.join(set(alignment[:,  i]))
        stereo_list = []
        #loop through each residue
        for res in range(0, len(column)):
            #replace the residue with list of properties
            residue = column[res]
            #append the properties list to a
            stereo_prop = dic_prop.get(residue)
            stereo_list.append(stereo_prop)
        #number of common properties
        count_stereo = sum(len(set(i)) == 1 for i in zip(*stereo_list))
        #add the number of properties to a list
        score_list.append(count_stereo)
    score_list_final = [float(i*0.1) for i in score_list]
    return score_list_final

if __name__ == '__main__':
    main()