import argparse
import sys
import Scoring
import numpy as np

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
    parser.add_argument("-q1", "--test",
                        help="""query sequence""",
                        dest="test",
                        required=True)
    parser.add_argument("-l1", "--lc_test",
                        help='lc profile aligned with query alignment',
                        dest='lc_test', 
                        required=True)
    # parser.add_argument("-q2", "--test_vir",
    #                     help="""query sequence""",
    #                     dest="test_vir",
    #                     required=True)
    # parser.add_argument("-l2", "--lc_test_vir",
    #                     help='lc profile aligned with query alignment',
    #                     dest='lc_test_vir', 
    #                     required=True)
    parser.add_argument("-o1",
                        dest="output1",
                        help="output image file for viral")
    parser.add_argument("-o2",
                        help="output image file for GTA",
                        dest="output2")
    parser.add_argument("-t",
                        dest="threshold",
                        type=float,
                        required=False,
                        help="threshold percentile value, between 0 and 100")
    parser.add_argument("-s", "--score",
                        dest="score_name",
                        required=True,
                        help="which metric would you like to use?")
    return parser

# return 1 for gta, 0 for virus
def knn(obsv, gta_clust, vir_clust, k):
    neighbors = []
    for memb in gta_clust:
        dist = euclid(obsv, memb)
        neighbors.append((dist, 1))
    for memb in vir_clust:
        dist = euclid(obsv, memb)
        neighbors.append((dist, 0))
    neighbors.sort()
    return neighbors


def classify(s, args):
    # g = Scoring.Scoring(args.gta, args.viral, args.test_gta, args.lc_test_gta, args.threshold)

    # list of arrays 
    gta_cluster = []
    vir_cluster = []
    sig_cols = s.get_sig_cols(args.score_name)
    # print sig_cols

    for col in sig_cols:
        gta_col = s.gta[:, col]
        vir_col = s.viral[:, col]
        gta_centr = s.consensus(gta_col, s.w_gta)
        vir_centr = s.consensus(vir_col, s.w_vir)
        gta_cluster.append(gta_centr)
        vir_cluster.append(vir_centr)

    col_map = s.get_map()

    classes = []
    for seq in s.query:
        votes = []
        for i in range(len(sig_cols)):
            x = sig_cols[i]
            q_x = col_map[x]
            v = s.consensus([seq[q_x]], [1])
            dist_gta = (s.euclid(v, gta_cluster[i]), 1)
            dist_vir = (s.euclid(v, vir_cluster[i]), 0)
            vote = min(dist_gta, dist_vir)[1]
            votes.append(vote)
        if sum(votes) >= len(votes)/2:
            classes.append(1)
        else:
            classes.append(0)
    # print classes
    return float(sum(classes))/float(len(classes))*100




def execute(args):
    g = Scoring.Scoring(args.gta, args.viral, args.test, args.lc_test, args.threshold)
    return classify(g, args)


def main():
    parser = cmdline_parser()
    args = parser.parse_args()
    print args.test, execute(args)

if __name__ == '__main__':
    main()









