'''Converts tables in SQLite into json files
	Attributes:
		database: path to sqlite database
		query: query execution command in string format
		
'''

import json
import os
import sqlite3
import pprint

class SQLToJSON():

    def __init__(self, database, query):
        self.database = database
        self.query = query

    def db(self):
        return sqlite3.connect(self.database)

    # returns all the table names in a sqlite database
    def table(self):
    	cursor = self.db().cursor()
    	cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    	tables = cursor.fetchall()
    	table_names = []
    	for table_name in tables:
        	table_names.append(table_name[0])
    	cursor.connection.close()
    	return table_names

    # one : return as one json file or not
    def query_table(self, table, one=True):
        cur = self.db().cursor()
        query = self.query + table
        cur.execute(query)
        r = [dict((cur.description[i][0], value)
                  for i, value in enumerate(row)) for row in cur.fetchall()]
        cur.connection.close()
        return (r if r else None) if one else r


# if __name__ == '__main__':
#     database = './cache/SabioRk.sqlite'
#     query = "select * from "
#     collection_dir = './cache/SabioRk/'
#     os.makedirs(os.path.dirname(collection_dir), exist_ok=True)

#     temp = SQLToJSON(database, query)
#     tables = temp.table()

#     for table in tables:
#         file_name = os.path.join(collection_dir + table + '.json')
#         result = SQLToJSON(database, query).query_table(table)
#         with open(file_name, "w") as f:
#                 f.write(json.dumps(result, indent=4))

