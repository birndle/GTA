from Classifier import Classifier
import argparse
from Bio import SeqIO

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
        parser.add_argument("-q", "--queries",
                            help="""query sequence""",
                            dest="queries",
                            required=False)
        parser.add_argument("-c", dest="C", required=False)
        parser.add_argument("-n",
                            help="number of features",
                            dest="n",
                            required=False)
        return parser


if __name__ == "__main__":
    parser = cmdline_parser()
    args = parser.parse_args()
    gta = list(SeqIO.parse(args.gta, "fasta"))
    viral = list(SeqIO.parse(args.viral, "fasta"))
    model = Classifier(gta, viral, n=int(args.n))
    model.get_training_set()
    SVs = model.learn_SVM_model(float(args.C))
    model.CV()

