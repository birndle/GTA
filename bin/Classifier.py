import argparse
import numpy 
import Queue
import Vectorizer
import Sequence
import math
import numpy as np
import re
from collections import defaultdict

from Bio import SeqIO, AlignIO

from sklearn import cross_validation
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


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
                            required=False)
        parser.add_argument("-k",
                            help="length of k-grams to build feature vectors",
                            dest="k",
                            required=False)
        parser.add_argument("-c",
                            help="length of k-grams to build feature vectors",
                            dest="c",
                            required=False)
        parser.add_argument("-a", "--align",
                            help="alignment of GTAs in training set",
                            dest="aln",
                            required=False)
        parser.add_argument("-o", "--output",
                            help="name of output file for writing information",
                            dest="out",
                            required=False)
        parser.add_argument("-n",
                            help="number of features",
                            dest="n",
                            required=False)

        return parser

class Classifier:

    def __init__(self, gta, viral, other=[], dist_matrix=None, n=1000):
        self.gta_seqs = gta
        # print '# of GTAs in training set: %d' % len(self.gta_seqs)
        self.viral_seqs = viral
        # print '# of Viruses in training set: %d' % len(self.viral_seqs)
        self.other_seqs = other

        self.training_seqs = self.gta_seqs + self.viral_seqs + self.other_seqs

        self.X = []
        self.y = np.array([1]*len(self.gta_seqs) + [0]*len(self.viral_seqs) + [2]*len(self.other_seqs))
        self.sparse_feat_list = []
        self.model = None
        
        self.label = {}
        self.label[2] = 'other'
        self.label[1] = 'gta'
        self.label[0] = 'virus'

        self.v = None
        self.n = n
        
        self.feats = None
        self.kmers = []
        self.gta_feats = None
        self.virus_feats = None

        self.gta_kmers = []
        self.virus_kmers =[]
        self.other_kmers = []

        self.matrix = dist_matrix
        self.weights = None


    def get_weights(self):
        for sample in self.training_seqs:
            print sample.id


    def classify(self, queries):
        predictions = []
        names = []
        test_set = []
        
        for rec in queries:
            name = rec.id   # get gene identifier number
            name = name.split('|')[1]
            names.append(name)
            test_seq = str(rec.seq)  # convert query to string
            s = Sequence.Sequence(test_seq) # convert query string to Sequence object 
            s.get_kmers(4, True) # count up k-mers in query sequence
            test, blah = self.v.vectorize(s) # convert query sequence to feature vector
            for j in self.sparse_feat_list: # prune sparse features from query sequence
                del test[j]
            test = np.array(test)
            test = test[self.feats] # extract only significant features
            test_set.append(test)
        
        test_set = np.array(test_set)
        predict = self.model.predict(test_set)
        predict = [self.label[guess] for guess in predict]
        hyperp_dist = self.model.decision_function(test_set)
        hyperp_dist = [x[0] for x in hyperp_dist]
        
        return (names, predict, hyperp_dist)

    def feature_visualization(self, aln, seqs, kmers, output='', write=False):
        
        # Define mapping from column positions in raw sequence to corresponding
        # position in the multiple sequence alignment
        def get_alignment_map(alignment, raw):
            maps = []
            for i in range(len(alignment)):
                raw_seq = list(raw[i].seq)
                aligned_seq = list(alignment[i].seq)
                # print aligned_seq
                x = 0
                col_map = {}
                for y in range(len(aligned_seq)):
                    if not aligned_seq[y] == '-':
                        col_map[x] = y
                        x +=1
                maps.append(col_map)
            return maps

        aln.sort(key=lambda x: x.id)
        seqs.sort(key=lambda x:x.id)
        col_maps = get_alignment_map(self.aln, seqs)

        all_cols = [] # will hold indicies for every column in the MSA containing at least one kmer occurence
        kmer_idx = []
        for i in range(len(seqs)):
            s = seqs[i]
            a = aln[i]
            col_map = col_maps[i]

            seq = str(s.seq)
            feat_idx = []
            for kmer in kmers:
                for match in re.finditer(kmer, seq):
                    idx = range(match.start(), match.end())
                    feat_idx.extend(idx)
            feat_idx.sort()

            unq_feat_idx = list(set(feat_idx)) 
            
            # stores the index of every position in the current sequence containing kmer
            unq_feat_idx.sort()
            unq_feat_idx = [col_map[i] for i in unq_feat_idx]
            all_cols.extend(unq_feat_idx)
            kmer_idx.append(unq_feat_idx)

        all_cols = list(set(all_cols))
        
        if write:
            for j in range(len(aln)):
                a = aln[j]
                personal_kmer_idx = kmer_idx[j]

                aln_seq_list = list(a.seq)
                edited_aln_seq = [aln_seq_list[i] if i in personal_kmer_idx else aln_seq_list[i].lower() if i in all_cols else '-' if aln_seq_list[i] == '-' else 'X' for i in range(len(aln_seq_list))]
                new_seq = ''.join(edited_aln_seq)

                # WRITE NEW_SEQ TO FASTA FILE, VIEW IN JALVIEW

                f = open(output, 'a')
                f.write("> %s\n" % (a.id))
                f.write(new_seq)
                f.write("\n")
                f.close()

        return all_cols


    def get_training_set(self):
        varied_kmers = True
        hashed=False
        self.v = Vectorizer.Vectorizer(4, 2**16, varied_kmers, hashed)
        
        if hashed:
            for rec in self.training_seqs:
                seq = str(rec.seq)
                s = Sequence.Sequence(seq)
                vec = v.vectorize(s)
                self.X.append(vec)
        else:
            bag = set()
            training_seqs = []

            # Compile bag of all unique kmers occuring in the reference sequences
            for rec in self.training_seqs:
                seq = str(rec.seq)
                s = Sequence.Sequence(seq)
                training_seqs.append(s)
                sack = s.get_kmers(4, True)
                bag.update(sack.keys())

            # Give this raw bag of words to the vectorizer
            self.v.set_bag(list(bag))

            # Vectorize the samples in our training sequences.
            # And get list of kmers corresponding to frequency counts in our feature vectors
            for examp in training_seqs:
                vec, kmers = self.v.vectorize(examp)
                self.X.append(vec)

        print 'Number of features pre-pruning: %d' % len(self.X[0])

        # Prune out sparse features from vectors and from list of kmers
        self.sparse_feat_list = self.v.get_sparse_features(self.X, 3)
        for i in self.sparse_feat_list:
            del kmers[i]
            for examp in self.X:
                del examp[i]

        self.X = np.array(self.X)

        print'Number of features post-pruning: %d' % len(self.X[0])

        # self.X, self.feats, self.gta_feats, self.virus_feats = self.v.select_features(self.X, self.y, n=self.n)
        
        self.feats = self.v.select_features(self.X, self.y, n=self.n, num_classes=3)
        kmers = np.array(kmers)

        # self.gta_kmers = kmers[self.gta_feats]
        # self.virus_kmers = kmers[self.virus_feats]
        # self.other_kmers = kmers[self.other_feats]

        self.kmers = kmers[self.feats]
        new_X = []
        for vec in self.X:
            vec = list(vec[self.feats])
            new_X.append(vec)
        self.X = new_X
        
        print'Number of features post-selection: %d' % len(self.X[0])
        # self.v.select_features(self.X, self.y, method='unan')
        return


    def learn_kNN_model(self, k):
        self.model = KNeighborsClassifier(n_neighbors=k).fit(self.X, self.y)
        SV_idx = self.model.support_
        SVs = [self.training_seqs[i] for i in SV_idx]
        return SVs


    # NEED WEIGHTING
    def learn_SVM_model(self, c, multi=False):
        if multi:
            from sklearn.multiclass import OneVsRestClassifier
            model = svm.SVC(kernel='linear', C=c, class_weight='auto')
            OVA_model = OneVsRestClassifier(model).fit(self.X, self.y)
            self.model = OVA_model
        else:
            self.model = svm.SVC(kernel='linear', C=c, class_weight='auto').fit(self.X, self.y, sample_weight=self.weights)
            SV_idx = self.model.support_
            SVs = [self.training_seqs[i] for i in SV_idx]
            return SVs

    def CV(self):
        kf = cross_validation.KFold(len(self.y), n_folds=4, shuffle=True)
        scores = cross_validation.cross_val_score(self.model, self.X, self.y, cv=kf)
        print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
        return scores

# NEEDS WEIGHTING
def conservation_score(aln, feat_cols):
    from Conservation import shannon

    general = []
    selected = []

    for x in range(len(aln[0].seq)):
        column = aln[:, x]
        column_lst = list(column)
        # Ignore gappy columns
        # if column_lst.count('-')/float(len(column_lst)) < 0.3:
        if column_lst.count('-') == 0:
            print column
            ent = shannon(column, stereo=True)
            general.append(ent)
            if x in feat_cols:
                selected.append(ent)
    print len(general), len(selected)
    print "Average entropy over all columns: %f" % (sum(general)/len(general))
    print "Average entropy over selected columns: %f" % (sum(selected)/len(selected))

def main():
    parser = cmdline_parser()
    args = parser.parse_args()
    gta = list(SeqIO.parse(args.gta, "fasta"))
    viral = list(SeqIO.parse(args.viral, "fasta"))
    aligned_gta = list(SeqIO.parse(args.aln, "fasta"))
    model = Classifier(gta, viral, other=aligned_gta, n=500)
    model.get_training_set()
    model.learn_SVM_model(1.0, multi=True)

    # feat_idx = model1.feature_visualization(aligned_gta, model1.gta_seqs, model1.gta_kmers)
    # conservation_score(AlignIO.read(args.aln, "fasta"), feat_idx)
    return
    
    


if __name__ == '__main__':
    main()
