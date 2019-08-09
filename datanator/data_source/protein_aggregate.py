import json
from datanator.util import mongo_util


class ProteinAggregate:

    def __init__(self, username=None, password=None, server=None, authSource='admin',
    	database='datanator'):
    	
    	self.mongo_manager = mongo_util.MongoUtil(MongoDB=server, username=username,
    		password=password, authSource=authSource, db=database)

     

