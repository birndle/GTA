import argparse

def cmdline_parser():
	parser = argparse.ArgumentParser(description=""" """)
	parser.add_argument("-f",                        
						help="""filename""",
                        dest="file",
                        required=True)
	return parser

def parse_txt(filename):
	f = open(filename, 'r')
	homol = ''
	data = []
	for line in f:
		if line.startswith('>'):
			words = line.split(' ')
			homol = words[0][1:-7]
		else:
			words = line.split('\t')
			id_data = words[0].split('|')
			gi = int(id_data[1])
			prediction = words[1][:-1]
			entry = (homol, gi, prediction)
			data.append(entry)
	return data


def create_delimited_data_file(out, data):
	f = open(out, 'a')
	for entry in data:
		f.write('%s\t%d\t%s\n' % entry)
	f.close()


def main():
	parser = cmdline_parser()
	args = parser.parse_args()
	data = parse_txt(args.file)
	create_delimited_data_file('data.csv', data)


if __name__ == '__main__':
    main()