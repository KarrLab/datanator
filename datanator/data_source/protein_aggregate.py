import json
from datanator.util import mongo_util
from datanator.core import (query_pax, query_kegg_orthology,
                           query_taxon_tree)
import datanator.config.core
from pymongo.collation import Collation, CollationStrength
import pymongo

class ProteinAggregate:

    def __init__(self, username=None, password=None, server=None, authSource='admin',
                 src_database='datanator', max_entries=float('inf'), verbose=True,
                 collection='protein', destination_database='datanator'):
        '''
                Args:
                        src_database (:obj: `str`): name of database in which source collections reside
                        destination_database (:obj: `str`): name of database to put the aggregated collection
        '''

        self.max_entries = max_entries
        self.verbose = verbose
        self.mongo_manager = mongo_util.MongoUtil(MongoDB=server, username=username,
                                                  password=password, authSource=authSource, db=src_database)
        self.pax_manager = query_pax.QueryPax(MongoDB=server, db=src_database,
                                              collection_str='pax', verbose=verbose, max_entries=max_entries, username=username,
                                              password=password, authSource=authSource)
        self.kegg_manager = query_kegg_orthology.QueryKO(server=server, database=src_database,
                                                         verbose=verbose, max_entries=max_entries, username=username,
                                                         password=password, authSource=authSource)
        self.taxon_manager = query_taxon_tree.QueryTaxonTree(collection_str='taxon_tree', 
                verbose=verbose, max_entries=max_entries, username=username, MongoDB=server, 
                password=password, db=src_database, authSource=authSource)
        self.client, self.db, self.col = mongo_util.MongoUtil(MongoDB=server, username=username,
                                                              password=password, authSource=authSource,
                                                              db=destination_database).con_db(collection)

    def copy_uniprot(self):
        '''
            Copy relevant information from uniprot collection
        '''
        _, _, col_uniprot = self.mongo_manager.con_db('uniprot')
        query = {}
        projection = {'status': 0, '_id': 0}
        docs = col_uniprot.find(filter=query, projection=projection)
        count = col_uniprot.count_documents({})
        self.col.insert_many(docs)
        collation = Collation(locale='en', strength=CollationStrength.SECONDARY)
        self.col.create_index([("uniprot_id", pymongo.ASCENDING)], background=True, collation=collation)

    def load_abundance_from_pax(self):
        '''
            Load protein abundance data but interating from pax collection.
        '''
        _, _, col_pax = self.mongo_manager.con_db('pax')
        query = {}
        projection = {'ncbi_id': 1, 'species_name': 1,
                    'observation': 1, 'organ': 1}
        docs = col_pax.find(filter=query, projection=projection, batch_size=5)
        count = col_pax.count_documents(query)
        progress = 193
        for i, doc in enumerate(docs[progress:]):            
            species_name = doc['species_name']
            taxon_id = doc['ncbi_id']
            organ = doc['organ']
            if i == self.max_entries:
                break
            if self.verbose and i % 1 == 0:
                print('Loading abundance info {} of {} ...'.format(
                    i + progress, min(count, self.max_entries)))
            for j, obs in enumerate(doc['observation']):
                if j == self.max_entries:
                    break
                if self.verbose and j % 100 == 0 and i % 1 == 0:
                    print('  Loading observation info {} of {} ...'.format(
                        j, len(doc['observation'])))
                try:
                    uniprot_id = obs['protein_id']['uniprot_id']
                    ordered_locus_name = obs['protein_id']['string_id']
                    abundance = obs['abundance']
                    dic = {'organ': organ, 'abundance': abundance}
                    self.col.update_one({'uniprot_id': uniprot_id},
                                        {'$set': {'ncbi_taxonomy_id': taxon_id,
                                                  'species_name': species_name,
                                                  'ordered_locus_name': ordered_locus_name},
                                         '$push': {'abundances': dic} })
                except TypeError:
                    continue



    def load_ko(self):
        '''
                Load ko number for uniprot_id if such information
                exists
        '''
        query = {}
        projection = {'uniprot_id': 1, 'gene_name': 1}
        docs = self.col.find(filter=query, projection=projection, batch_size=1000)
        count = self.col.count_documents(query)
        progress = 228330
        for i, doc in enumerate(docs[progress:]):
            if i == self.max_entries + 2:  # for testing script
                break
            if self.verbose and i % 10 == 0:
                print('Loading KO info {} of {} ...'.format(
                    i + progress, min(count, self.max_entries)))

            gene_name = doc['gene_name']
            ko_number = self.kegg_manager.get_ko_by_name(gene_name)
            if ko_number != None:
                self.col.update_one({'uniprot_id': doc['uniprot_id']},
                                    {'$set': {'ko_number': ko_number}})

    def load_taxon(self):
        '''
        	Load taxon ancestor information
        '''
        query = {'ncbi_taxonomy_id': {'$exists': True}}
        projection = {'ncbi_taxonomy_id': 1, 'uniprot_id': 1}
        docs = self.col.find(filter=query, projection=projection, batch_size=1000)
        count = self.col.count_documents(query)
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            if self.verbose and i % 10 == 0:
                print('Loading taxon info {} of {} ...'.format(
                    i, min(count, self.max_entries)))
            taxon_id = doc['ncbi_taxonomy_id']
            anc_id, anc_name = self.taxon_manager.get_anc_by_id([taxon_id])
            self.col.update_one({'uniprot_id': doc['uniprot_id']},
            					{'$set': {'ancestor_name': anc_name[0],
            							'ancestor_taxon_id': anc_id[0]} })

def main():
    src_db = 'datanator'
    des_db = 'datanator'
    collection_str = 'protein'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    port = datanator.config.core.get_config(
    )['datanator']['mongodb']['port']        
    manager = ProteinAggregate(username=username, password=password, server=server, 
                               authSource='admin', src_database=src_db,
                               verbose=True, collection=collection_str, destination_database=des_db)
    # manager.copy_uniprot()
    # manager.load_abundance_from_pax()
    manager.load_ko()
    manager.load_taxon()

    # collation = Collation(locale='en', strength=CollationStrength.SECONDARY)
    # manager.collection.create_index([("uniprot_id", pymongo.ASCENDING),
    #                        ("ancestor_taxon_id", pymongo.ASCENDING)], background=True, collation=collation)
    # manager.collection.create_index([("ko_number", pymongo.ASCENDING),
    #                        ("ncbi_taxonomy_id", pymongo.ASCENDING)], background=True, collation=collation)

if __name__ == '__main__':
	main()