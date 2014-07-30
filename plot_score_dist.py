#!/usr/bin/env python
#File created on May 19 2014


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
import brewer2mpl
import matplotlib as mpl  # for fonts
mpl.rcParams['pdf.fonttype'] = 42 #for making fonts readable by illustrator

#--- project specific imports
import tables


def cmdline_parser():
    """
    creates an argparse instance
    """
    parser = argparse.ArgumentParser(description=""" """)
    parser.add_argument('-t1', '--viral',
                        help='alphaviral scores',
                        required=True)
    parser.add_argument('-t2', '--lc',
                        help='GTA from complete genome scores',
                        required=True)
    parser.add_argument('-t3', '--sc_comp',
                        help='small cluster from complete genome',
                        required=True)
    parser.add_argument('-t4', '--sc_draft',
                        help='small cluster alignments from draft genomes',
                        required=True)
    parser.add_argument('-t5', '--lc_draft',
                        help='large cluster alignments from draft genomes',
                        required=True)
    parser.add_argument('-t6', '--prophage',
                        help='prophage alignments from complete genomes',
                        required=True)
    parser.add_argument('-o', '--output',
                        help='output pdf files',
                        required=True)

    return parser


def main():
    """
    The main function
    """
    parser = cmdline_parser()
    args = parser.parse_args()

    viral_dist = col2list(args.viral, 1)
    lc_dist = col2list(args.lc, 1)
    prophage_dist = col2list(args.prophage, 1)
    sc_comp_dist = col2list(args.sc_comp, 1)
    sc_draft_dist = col2list(args.sc_draft, 1)
    lc_draft_dist = col2list(args.lc_draft, 1)

    fig = plt.figure()
    axis1 = fig.add_subplot(111)
    bmap = brewer2mpl.get_map('Set1', 'qualitative', 6)
    box_colors = bmap.mpl_colors

    width = 0.14
    N = 10
    ind = np.arange(N)

    hist, bins = np.histogram(viral_dist, bins=10, range=(0, 1))
    plt.bar(ind+width, hist.astype(np.float32), width, color=box_colors[0], label='alpha-proteobacteria viral refseq and marine viral metagenome')

    hist, bins = np.histogram(lc_dist, bins=10, range=(0, 1))
    bins = bins+width
    plt.bar(ind+width+width, hist.astype(np.float32), width, color=box_colors[1], label='8 other GTAs in surrounding genomic region ')

    hist, bins = np.histogram(prophage_dist, bins=10, range=(0, 1))
    bins = bins+width+width
    plt.bar(ind+width+width+width, hist.astype(np.float32), width, color=box_colors[2], label='prophages (PhiSpy) without GTA homologs')

    hist, bins = np.histogram(sc_comp_dist, bins=10, range=(0, 1))
    bins = bins+width+width+width
    plt.bar(ind+width+width+width+width, hist.astype(np.float32), width, color=box_colors[3], label='less than 8 other GTAs in surrounding genomic region')

    hist, bins = np.histogram(lc_draft_dist, bins=10, range=(0, 1))
    bins = bins+width+width+width+width
    plt.bar(ind+width+width+width+width+width, hist.astype(np.float32), width, color=box_colors[4], label='8 other GTAs in surrounding genomic region')

    hist, bins = np.histogram(sc_draft_dist, bins=10, range=(0, 1))
    bins = bins+width+width+width+width
    plt.bar(ind+width+width+width+width+width+width, hist.astype(np.float32), width, color=box_colors[5], label='less than 8 other GTAs in surrounding genomic region')


    plt.xticks(ind, np.arange(0, 1.1, 0.1))
    #plt.yscale('symlog', nonposy='clip', basey=2)
    plt.xlabel('Score distribution')
    plt.ylabel('Number of homologs')
    #for shading region
    #plt.axhspan(0, 2, facecolor='0.5', alpha=0.25)
    legend = axis1.legend(loc='upper left', shadow=False, prop={'size': 10})

    fig = axis1.get_figure()
    fig.savefig(args.output)


def frange(x, y, jump):
    """
    function to generate floats
    """
    while x < y:
        yield x
        x += jump


def col2list(filename, colnum):
    """
    this function converts a column from tab delimited file to a list

    and returns floats
    """
    data_list = []
    for line in open(filename, 'r'):
        data_list.append(line.split(',')[colnum].rstrip())
    data_float_list =[float(i) for i in data_list]
    return data_float_list


if __name__ == '__main__':
    main()
