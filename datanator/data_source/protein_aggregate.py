import json
from datanator.util import mongo_util
from datanator.core import query_pax, query_kegg_orthology


class ProteinAggregate:

    def __init__(self, username=None, password=None, server=None, authSource='admin',
    	database='datanator', max_entries=float('inf'), verbose=True):

    	self.max_entries = max_entries
    	self.verbose = verbose	
    	self.mongo_manager = mongo_util.MongoUtil(MongoDB=server, username=username,
    		password=password, authSource=authSource, db=database)
    	self.pax_manager = query_pax.QueryPax(MongoDB=server, db=database,
                 collection_str='pax', verbose=verbose, max_entries=max_entries, username=username,
                 password=password, authSource=authSource)

    def load_content(self, aggregate='protein'):
    	_, _, col_uniprot = self.mongo_manager.con_db('uniprot')
    	_, _, col_protein = self.mongo_manager.con_db(aggregate)
    	query = {}
    	projection = {'status': 0, '_id': 0}
    	docs = col_uniprot.find_one(filter=query, projection=projection)
    	count = col_uniprot.count_documents({})
    	for i, doc in enumerate(docs):
    		if i == max_entries:
    			break
    		if self.verbose and i % 10 == 0:
    			print('Processing doc {} of {} ...'.format(i, count))
    		uniprot_id = doc['uniprot_id']
    		abundances = self.pax_manager.get_abundance_from_uniprot(uniprot_id)
    		doc['abundances'] = abundances
    		col_protein.update_one({'uniprot_id': doc['uniprot_id']},
    			doc,
    			upsert=True)