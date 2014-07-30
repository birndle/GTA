from Classifier import Classifier
import argparse
import MyDB
import parse_dists
from peewee_script import insert_new_rows

from Bio import SeqIO

def cmdline_parser():
        """
        creates an argparse instance
        """
        parser = argparse.ArgumentParser(description=""" """)

        parser.add_argument("-g", "--gta",
                            help="""gta sequences""",
                            dest="gta",
                            required=False)
        parser.add_argument("-v", "--viral",
                            help="""viral sequences""",
                            dest="viral",
                            required=False)
        parser.add_argument("-q", "--queries",
                            help="""query sequence""",
                            dest="queries",
                            required=False)
        parser.add_argument("-c", dest="C", required=False)
        return parser


model = None

if __name__ == "__main__":
    parser = cmdline_parser()
    args = parser.parse_args()
    gta = list(SeqIO.parse(args.gta, "fasta"))
    viral = list(SeqIO.parse(args.viral, "fasta"))
    model = Classifier(gta, viral)
    queries = args.queries.split(',')

    for query in queries:
        query_seqs = list(SeqIO.parse(query, "fasta"))
        gene_num = int(query[query.find('orfg')+4])
        if not model:
            # dist_matrix = parse_dists.get_dist_matrix(gene_num)
            model = Classifier(gta, viral)
            model.get_training_set()
            # model.get_weights()
            SVs = model.learn_SVM_model(float(args.C))
            model.feature_visualization()

        # results = zip(*model.classify(query_seqs))
        # insert_new_rows(results, args.homo)



