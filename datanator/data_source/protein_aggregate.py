import json
from datanator.util import mongo_util
from datanator.core import query_pax


class ProteinAggregate:

    def __init__(self, username=None, password=None, server=None, authSource='admin',
    	database='datanator', max_entries=float('inf'), verbose=True):

    	self.max_entries = max_entries
    	self.verbose = verbose	
    	self.mongo_manager = mongo_util.MongoUtil(MongoDB=server, username=username,
    		password=password, authSource=authSource, db=database)

    def load_content(self):
    	_, _, col_uniprot = self.mongo_manager.con_db('uniprot')
    	query = {}
    	projection = {'status': 0, '_id': 0}