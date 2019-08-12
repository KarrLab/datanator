import json
from datanator.util import mongo_util
from datanator.core import (query_pax, query_kegg_orthology,
                           query_taxon_tree)


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

    def load_abundance(self):
        '''
                Load protein abundance data from paxDB
        '''
        _, _, col_uniprot = self.mongo_manager.con_db('uniprot')
        query = {}
        projection = {'status': 0, '_id': 0}
        docs = col_uniprot.find(filter=query, projection=projection)
        count = col_uniprot.count_documents({})
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            if self.verbose and i % 10 == 0:
                print('Loading abundance info {} of {} ...'.format(
                    i, min(count, self.max_entries)))
            uniprot_id = doc['uniprot_id']
            abundances = self.pax_manager.get_abundance_from_uniprot(
                uniprot_id)
            if abundances != []:
                doc['abundances'] = abundances[1:]
                doc['ncbi_taxonomy_id'] = abundances[0]['ncbi_taxonomy_id']
                doc['species_name'] = abundances[0]['species_name']
            self.col.update_one({'uniprot_id': doc['uniprot_id']},
                                {'$set': doc},
                                upsert=True)

    def load_ko(self):
        '''
                Load ko number for uniprot_id if such information
                exists
        '''
        query = {}
        projection = {'uniprot_id': 1, 'gene_name': 1}
        docs = self.col.find(filter=query, projection=projection)
        count = self.col.count_documents(query)
        for i, doc in enumerate(docs):
            if i == self.max_entries + 1:  # for testing script
                break
            if self.verbose and i % 10 == 0:
                print('Loading KO info {} of {} ...'.format(
                    i, min(count, self.max_entries)))

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
        docs = self.col.find(filter=query, projection=projection)
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
