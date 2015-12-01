from parse_predictions import parse_txt, cmdline_parser
import pymysql

class MyDB:

	def __init__(self, host, user, passwd, dbname):
		self.conn, self.cursor = self.login(host, user, passwd, dbname)
		self.table = ''

	def login(self, host, user, passwd, dbname):
		conn = pymysql.connect(host=host, user=user, passwd=passwd, db=dbname)
		curr = conn.cursor()
		return [conn, curr]


	""" Creates a new table in the MySQL database
	fields: a list of tuples specifyihng the columns you wish to add to the table
	i.e. (field_name, var_type [,col_length]) 

	Note: col_length should be an integer surrounded by parentheses

	e.g. ('gid', 'INT', '(15)')
	"""

	def create_table(self, table_name, fields):
		query = 'CREATE TABLE %s (' % table_name
		for field in fields:
			query += ('%s '*len(field) % field)
			query = query[:-1] + ', '
		query = query[:-2]
		query += ');'
		
		self.cursor.execute(query)


	def populate_table_from_delim_file(self, data_file):
		query = 'LOAD DATA LOCAL INFILE \'%s\' INTO TABLE %s;' % (data_file, self.table)
		# print query
		self.cursor.execute(query)

	def close_connection(self):
		self.cursor.close()
		self.conn.close()

	""" data is a list of tuples where the first value in the tuple is the data to be 
	added to the updated table and the remaining values in the tuple are qualifying statements
	that determine which row the data will be added to, i.e.:

	(0.92, 'gid', 01935969) --> UPDATE table SET col = gta WHERE gid=01935969;

	"""
	def update_col(self, col, data):
		query = 'UPDATE %s SET %s = \'%f\' WHERE %s=\'%s\';' % (self.table, col, data[0], data[1], data[2])
		self.cursor.execute(query)


	def set_table(self, table_name):
		self.table = table_name


	def add_row(self, values):
		values = tuple(values)
		query = 'INSERT INTO %s VALUES (' % self.table
		parenthetical = '\'%s\','*len(values) % values
		parenthetical = parenthetical[:-1]
		q = query + parenthetical + ');'
		self.cursor.execute(q)


	def update_row(self, val, ID):
		query = 'UPDATE %s SET 1000_best_BNS_selected_feats = \'%s\' WHERE gid = \'%d\';' % (self.table, val, ID)
		self.cursor.execute(query)


	def add_col(self, column_descip):
		query = 'ALTER TABLE %s ADD %s' % (self.table, column_descip)
		self.cursor.execute(query)

	def execute(self):
		query = 'ALTER TABLE MLPredictions CHANGE prediction all_feats_minus_sparse CHAR(10)'
		self.cursor.execute(query)

def create_delimited_data_file(out, data):
	f = open(out, 'a')
	for entry in data:
		f.write('%s\t%d\t%s\n' % entry)
	f.close()


def main():
	# data = parse_txt('svm-predictions.txt')
	# create_delimited_data_file('data.csv', data)
	db = mydb(host='localhost', user='mshakya', passwd='f001r60', dbname='mshakya')
	# db.create_table('MLPredictions', 
	# [('homolog', 'CHAR', '(7)'), ('gid', 'INT'), ('prediction', 'CHAR', '(10)')])
	db.set_table('MLPredictions')
	# dat = open('data.csv', 'r')
	# for line in dat:
	# 	vals = line.split()
	# 	db.update_row(vals[2], int(vals[1]))
	db.add_col('Confidence_1000_feats FLOAT')
	db.close_connection()
	# db.execute()






if __name__ == '__main__':
    main()


