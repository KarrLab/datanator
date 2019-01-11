
import datetime
import io
import json
from datanator.core import data_source
from datanator.util import molecule_util
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm

Base = sqlalchemy.ext.declarative.declarative_base()


class ymdb(data_source.HttpDataSource):

    def load_content(self):
        self.FULL_DB_URL = 'http://www.ymdb.ca/system/downloads/current/ymdb.json.zip'
        self.directory = 'YMDB/ymdb.json'
        with open(self.directory) as json_file:
            self.entries = json.load(json_file)
        return self.entries

    def fill_sqlite(self, file):

        # def main():
        # 	entries = load_content()
        # 	for ymdb in entries:
        # 		print(ymdb['ymdb_id'])
