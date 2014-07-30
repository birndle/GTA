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
from Zvelibil_score import stereo_score


class Scoring:

    def __init__(self, gta_aln, vir_aln, query, gta_query, thresh):

        self.viral = AlignIO.read(open(vir_aln, 'r'), "fasta")
        self.gta = AlignIO.read(open(gta_aln, 'r'), "fasta")
        self.gta_query = AlignIO.read(open(gta_query, 'r'), "fasta")
        self.query = AlignIO.read(open(query, 'r'), "fasta")
        self.t = thresh
        self.w_gta = []
        self.w_vir = []

        matrix =                   [[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,0,0],
                                   [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,0,0],
                                   [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,0,0],
                                   [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,0,0],
                                   [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,0,0],
                                   [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,0],
                                   [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,0,0],
                                   [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,0,0],
                                   [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0],
                                   [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,0,0],
                                   [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,0,0],
                                   [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,0],
                                   [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,0,0],
                                   [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,0,0],
                                   [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,0,0],
                                   [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0],
                                   [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,0,0],
                                   [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,0,0],
                                   [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,0,0],
                                   [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,0,0],
                                   [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                                   [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                matrix[i][j] = float(matrix[i][j])
                if i == j and not i == len(matrix):
                    matrix[i][j] = 1

        #get substitution matrix
        self.subst = SubstitutionMatrix("ARNDCQEGHILKMFPSTWYVX-", matrix)
        self.alphabet = "ARNDCQEGHILKMFPSTWYVX-"


    def henikoff_weighting(self, aln):
        L = aln.get_alignment_length()
        N = len(aln[:, 0])
        weights = [0]*N
        for x in range(0, L):
            col = aln[:, x]
            d = defaultdict(int)
            for res in col:
                d[res] += 1
            k = len(d)
            for i in range(len(col)):
                res = col[i]
                n = d[res]
                weights[i] += 1.0/(k*n)
        weights = [float(w)/L for w in weights]
        return weights


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


    def consensus(self, col, w):
        N = len(col)
        consensus = np.array([0]*len(self.alphabet)) # initialize consensus vector
        for i in range(len(col)):
            res = col[i]
            # build vector representation for each amino acid in the column
            v_res = np.array([self.subst[(let, res)] for let in self.alphabet])*w[i]
            consensus = consensus + v_res   
        consensus /= N
        return consensus

    def t_variance(self, col, w):
        centroid = self.consensus(col, w)
        var = 0
        N = len(col)
        for i in range(len(col)):
            res = col[i]
            weight = w[i]
            x = self.consensus([res], [1])
            var += self.euclid(x, centroid)*weight
        return var/N


    def thompson(self):
        scores = []
        for x in range(0, self.viral.get_alignment_length()):
            gta_col = self.gta[:, x]
            viral_col = self.viral[:, x]
            cons1 = self.consensus(gta_col, self.w_gta)
            cons2 = self.consensus(viral_col, self.w_vir)
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


    def scoring_distr(self, aln, weights):
        d = defaultdict(lambda: defaultdict(float))
        N = len(aln[:, 0]) # number of sequences in alignment
        for x in range(0, aln.get_alignment_length()):
            col = aln[:, x]
            for j in range(len(col)): 
                res = col[j]
                d[x][res] += weights[j]
            # for res in d[x]:
            #     d[x][res] /= N
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
            col1 = self.gta[:, x]
            col2 = self.viral[:, x]
            score = 0
            for i in range(len(col1)):
                for j in range(len(col2)):
                    a = col1[i]
                    b = col2[j]
                    w_i = self.w_gta[i]
                    w_j = self.w_vir[j]
                    score += (w_i + w_j)*self.subst[(a, b)]
            scores.append(score)
        return scores

    
    def get_sig_cols(self, scoring_fcn):
        self.w_gta = self.henikoff_weighting(self.gta)
        self.w_vir = self.henikoff_weighting(self.viral)
        scores = getattr(self, scoring_fcn)()

        # find gappy columns
        gappy_cols = gap_list(self.viral, self.gta, 0.7)

        # extract the position indicies for columns with scores exceeding the user-defined threshold
        # sig_cols = [i for i, j in enumerate(scores) if j > np.percentile(scores, self.t)]
        sig_cols = [i for i, j in enumerate(scores) if j > np.percentile(scores, self.t)]

        # remove gappy columns from signature column list
        sig_cols = [x for x in sig_cols if x not in gappy_cols]

        # get web logos
        # wl_gta = self.extract_sig_cols(self.gta, sig_cols)
        # wl_vir = self.extract_sig_cols(self.viral, sig_cols)

        # wl_gta_fasta = AlignIO.write(wl_gta, "gta%f.fasta" % (self.t), "fasta")
        # wl_vir_fasta = AlignIO.write(wl_vir, "vir%f.fasta" % (self.t), "fasta")
        # self.do_seqlogo("gta%f.fasta" % (self.t), "gta%f.pdf" % (self.t), sig_cols)
        # self.do_seqlogo("vir%f.fasta" % (self.t), "vir%f.pdf" % (self.t), sig_cols)

        return sig_cols

    def execute(self, scoring_fcn):
        
        self.w_gta = self.henikoff_weighting(self.gta)
        self.w_vir = self.henikoff_weighting(self.viral)
        scores = getattr(self, scoring_fcn)()

        # find gappy columns
        gappy_cols = gap_list(self.viral, self.gta, 0.7)

        # extract the position indicies for columns with scores exceeding the user-defined threshold
        sig_cols = [i for i, j in enumerate(scores) if j > np.percentile(scores, self.t)]
        # remove gappy columns from signature column list
        sig_cols = [x for x in sig_cols if x not in gappy_cols]

        # create table of frequencies for residues in reference gta alignment
        # will use these frequencies to score query sequence
        gta_reference_table = self.scoring_distr(self.gta, self.w_gta)  # frequency table based on initial gta alignment
        gta_query_reference_table = self.scoring_distr(self.gta_query, self.w_gta)  # gappy table after aligned with query

        # need a map from signature column indices to corresponding query columns 
        reference_to_query_col_map = comp_dict1(gta_query_reference_table, gta_reference_table)

        return self.score(self.query, gta_reference_table, reference_to_query_col_map, sig_cols)        

    def get_map(self):
        gta_reference_table = self.scoring_distr(self.gta, self.w_gta)  # frequency table based on initial gta alignment
        gta_query_reference_table = self.scoring_distr(self.gta_query, self.w_gta)  # gappy table after aligned with query

        # need a map from signature column indices to corresponding query columns 
        reference_to_query_col_map = comp_dict1(gta_query_reference_table, gta_reference_table)

        return reference_to_query_col_map 


    def zvelibil(self):
        z = stereo_score(self.gta)
        return z

