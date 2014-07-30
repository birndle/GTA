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
from SP_score import sim_score
from SubstitutionMatrix import SubstitutionMatrix


class Thompson:

    def __init__(self, gta_aln, vir_aln, query, gta_query, thresh):

        self.viral = AlignIO.read(open(vir_aln, 'r'), "fasta")
        self.gta = AlignIO.read(open(gta_aln, 'r'), "fasta")
        self.gta_query = AlignIO.read(open(gta_query, 'r'), "fasta")
        self.query = AlignIO.read(open(query, 'r'), "fasta")
        self.t = thresh

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

        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if i == j and not i == len(matrix):
                    matrix[i][j] = 1

        #get substitution matrix
        self.subst = SubstitutionMatrix("ARNDCQEGHILKMFPSTWYV-", matrix)
        self.alphabet = "ARNDCQEGHILKMFPSTWYV-"



    def do_seqlogo(self, input_fasta, output_seqlogo, col_list):
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


    def euclid(self, cg1, cg2):
        dis = 0.0
        dis = sum(np.power((cg1-cg2), 2))
        dis = math.sqrt(dis)
        return dis


    def consensus(self, col, subst, alphabet):
        N = len(col)
        consensus = np.array([0]*len(alphabet)) # initialize consensus vector
        for res in col:
            # build vector representation for each amino acid in the column
            v_res = np.array([subst[(let, res)] for let in alphabet])
            consensus += v_res    
        consensus /= N
        return consensus


    def thompson():
        scores = []
        for x in range(0, self.viral.get_alignment_length()):
            gta_col = self.gta[:, x]
            viral_col = self.viral[:, x]
            cons1 = self.consensus(gta_col, self.subst, self.alphabet)
            cons2 = self.consensus(viral_col, self.subst, self.alphab)
            scores.append(self.euclid(cons1, cons2))
        return scores


    def extract_sig_cols(self, aln, sig_cols):
        cols = []
        for x in sig_cols:
            cols.append(aln[:, x])

        zipped_cols = zip(*cols)
        parsed_seqs = ["".join(zipped_cols[i]) for i in range(len(zipped_cols))]

        records = [sr(Seq(s)) for s in parsed_seqs]
        wl = msa(records)   # MSA object containting only the signature columns
        return wl


    def scoring_distr(self, aln):
        d = defaultdict(lambda: defaultdict(float))
        N = len(aln[:, 0]) # number of sequences in alignment
        for x in range(0, aln.get_alignment_length()):
            for res in aln[:, x]: 
                d[x][res] += 1
            for res in d[x]:
                d[x][res] /= N
        return d


    def score(self, query, reference, col_map, sig_cols):
        scores = [] 
        for seq in query:
            the_score = 0
            for x in sig_cols:
                q_x = col_map[x]
                the_score += reference[x][seq[q_x]] 
            scores.append(the_score)
        scores = [x/len(sig_cols) for x in scores] 
        return scores


    def sp(self):
        scores = []
        for x in range(0, self.viral.get_alignment_length()):
            gta_col = self.gta[:, x]
            viral_col = self.viral[:, x]
            scores.append(sim_score(gta_col, viral_col, self.subst))
        return scores

    def execute_sp(self):
        
        scores = self.sp()

        # find gappy columns
        gappy_cols = gap_list(self.viral, self.gta, 0.7)

        # extract the position indicies for columns with scores exceeding the user-defined threshold
        sig_cols = [i for i, j in enumerate(scores) if j > np.percentile(scores, self.t)]

        # remove gappy columns from signature column list
        sig_cols = [x for x in sig_cols if x not in gappy_cols]

        # create table of frequencies for residues in reference gta alignment
        # will use these frequencies to score query sequence
        gta_reference_table = self.scoring_distr(self.gta)  # frequency table based on initial gta alignment
        gta_query_reference_table = self.scoring_distr(self.gta_query)  # gappy table after aligned with query

        # need a map from signature column indices to corresponding query columns 
        reference_to_query_col_map = comp_dict1(gta_query_reference_table, gta_reference_table)

        return self.score(self.query, gta_reference_table, reference_to_query_col_map, sig_cols)        


    def execute_t(self):
        """
        The main function
        """

        scores = self.thompson()

        # find gappy columns
        gappy_cols = gap_list(self.viral, self.gta, 0.7)

        # extract the position indicies for columns with scores exceeding the user-defined threshold
        sig_cols = [i for i, j in enumerate(scores) if j > np.percentile(scores, self.t)]

        # remove gappy columns from signature column list
        sig_cols = [x for x in sig_cols if x not in gappy_cols]

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
        gta_reference_table = self.scoring_distr(self.gta)  # frequency table based on initial gta alignment
        gta_query_reference_table = self.scoring_distr(self.gta_query)  # gappy table after aligned with query

        # need a map from signature column indices to corresponding query columns 
        reference_to_query_col_map = comp_dict1(gta_query_reference_table, gta_reference_table)

        query_scores = self.score(self.query, gta_reference_table, reference_to_query_col_map, sig_cols)

        return query_scores
