import argparse
import sys
import Thompson
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
    parser.add_argument("-q1", "--test_gta",
                        help="""query sequence""",
                        dest="test_gta",
                        required=True)
    parser.add_argument("-l1", "--lc_test_gta",
                        help='lc profile aligned with query alignment',
                        dest='lc_test_gta', 
                        required=True)
    parser.add_argument("-q2", "--test_vir",
                        help="""query sequence""",
                        dest="test_vir",
                        required=True)
    parser.add_argument("-l2", "--lc_test_vir",
                        help='lc profile aligned with query alignment',
                        dest='lc_test_vir', 
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
                        required=False,
                        help="threshold percentile value, between 0 and 100")
    parser.add_argument("-s", "--score",
                        dest="score",
                        required=True,
                        help="which metric would you like to use?")
    return parser


def thompson(gta, viral, query, gta_query, threshold):
	t = Thompson.Thompson(gta, viral, query, gta_query, threshold)
	return t.execute_t()


def sp(gta, viral, query, gta_query, threshold):
	sp = Thompson.Thompson(gta, viral, query, gta_query, threshold)
	return sp.execute_sp()
	

def main():
    parser = cmdline_parser()
    args = parser.parse_args()
    best = 0
    optimal = 0
    for t in range(0, 100, 1):
    	score = 0
    	gta_scores = getattr(sys.modules[__name__], args.score)(args.gta, args.viral, args.test_gta, args.lc_test_gta, t)
    	vir_scores = getattr(sys.modules[__name__], args.score)(args.gta, args.viral, args.test_vir, args.lc_test_vir, t)	
    	gta = np.array(gta_scores)
    	vir = np.array(vir_scores)
    	score = np.average(gta) - np.average(vir)
    	if score > best:
    		optimal = t
    		best = score
    print optimal, best

if __name__ == '__main__':
    main()
