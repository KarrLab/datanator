import json
from datanator.util import mongo_util
from datanator.core import (query_pax, query_kegg_orthology,
                           query_taxon_tree, query_protein)
import datanator.config.core
from pymongo.collation import Collation, CollationStrength
import pymongo
import os

class ProteinAggregate:

    def __init__(self, username=None, password=None, server=None, authSource='admin',
                 src_database='datanator', max_entries=float('inf'), verbose=True,
                 collection='protein', destination_database='datanator', cache_dir=None):
        '''
                Args:
                        src_database (:obj: `str`): name of database in which source collections reside
                        destination_database (:obj: `str`): name of database to put the aggregated collection
        '''

        self.max_entries = max_entries
        self.verbose = verbose
        self.cache_dir = cache_dir
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
        self.protein_manager = query_protein.QueryProtein(username=username, password=password, 
            server=server, collection_str='protein', max_entries=max_entries, database=src_database)
        self.client, self.db, self.col = mongo_util.MongoUtil(MongoDB=server, username=username,
                                                              password=password, authSource=authSource,
                                                              db=destination_database).con_db(collection)
        self.bad_kinlawid = [24416,24417,24418,24419,24420,24421,24422,24423]

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
        progress = 0
        for i, doc in enumerate(docs[progress:]):            
            species_name = doc['species_name']
            taxon_id = doc['ncbi_id']
            organ = doc['organ']
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
                                         '$push': {'abundances': dic} }, upsert=True)
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
        progress = 0
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

    def load_unreviewed_abundance(self):
        '''
            Load abundance info for proteins that are not reviewed in uniprot
        '''
        query_uniprot = {'abundances': {'$exists': True}}
        projection = {'uniprot_id': 1}
        reviewed = []
        docs = self.col.find(filter=query_uniprot, projection=projection)
        count = self.col.count_documents(query_uniprot)
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                pass
            if self.verbose and i % 1000 == 0:
                print('Getting reviewed protein with abundance {} of {} ...'.format(
                    i, min(count, self.max_entries)))
            reviewed.append(doc['uniprot_id'])

        total = self.pax_manager.collection.distinct('observation.protein_id.uniprot_id')
        unreviewed = list(set(total) - set(reviewed))

        count = len(unreviewed)
        for i, _id in enumerate(unreviewed):
            doc = {}            
            if i == self.max_entries:
                break
            if self.verbose and i % 100 == 0:
                print('Loading unreviewed protein abundance info {} of {} ...'.format(
                    i, min(count, self.max_entries)))
            abundances = self.pax_manager.get_abundance_from_uniprot(_id)
            if abundances != []:
                doc['abundances'] = abundances[1:]
                doc['ncbi_taxonomy_id'] = abundances[0]['ncbi_taxonomy_id']
                doc['species_name'] = abundances[0]['species_name']
                doc['ordered_locus_name'] = abundances[0]['ordered_locus_name']

            self.col.update_one({'uniprot_id': _id},
                                {'$set': doc}, upsert=True) 

    def load_kinlaw_from_sabio(self):
        '''
            load kinlaw_id from sabio_rk collection based on uniprot_id or
            protein name if uniprot_id information is not present
        '''
        _, _, col_sabio = self.mongo_manager.con_db('sabio_rk_new')
        projection = {'enzyme': 1, 'kinlaw_id': 1, 'taxon': 1}
        docs = col_sabio.find({}, projection=projection)
        count = col_sabio.count_documents({})
        collation = Collation('en', strength=CollationStrength.SECONDARY)
        progress = 0
        for i, doc in enumerate(docs[progress:]):
            if i == self.max_entries:
                break
            if self.verbose and i % 50 == 0:
                print('Processing Kinetics information doc {} out of {}'.format(i+progress,min(count, self.max_entries)))

            kinlaw_id = doc['kinlaw_id']
            if kinlaw_id in self.bad_kinlawid:
                continue

            enzyme = doc.get('enzyme')
            if enzyme == None or len(enzyme) > 1 :
                with open(self.cache_dir, 'a+') as f:
                    f.write('\n  There are more than 1 or no enzyme in kinetic law {}'.format(kinlaw_id))
                continue

            subunits = enzyme[0]['subunits']
            
            taxon = doc['taxon']

            if len(subunits) == 0:
                name = enzyme[0]['name']
                if name != None:
                    results = self.protein_manager.get_id_by_name(name)
                    for result in results:
                        query = {'uniprot_id': result['uniprot_id']}
                        self.col.update_one(query, {'$push': {'kinetics': {'kinlaw_id': kinlaw_id, 'ncbi_taxonomy_id': taxon} } }, 
                            collation=collation)
                else:
                    with open(self.cache_dir, 'a+') as f:
                        f.write('\n  Enzyme in kinetic law with ID {} has no name or uniprot_id'.format(kinlaw_id))                    

            else:
                proteins = []
                for subunit in subunits:
                    proteins.append(subunit['uniprot'])
                query = {'uniprot_id': {'$in': proteins}}
                projection = {'uniprot_id': 1}
                protein_docs = self.col.find(filter=query, projection=projection, collation=collation)
                for protein_doc in protein_docs:
                    self.col.update_one({'uniprot_id': protein_doc['uniprot_id']},
                                        {'$push': {'kinetics': {'kinlaw_id': kinlaw_id, 'ncbi_taxonomy_id': taxon} } }, 
                                        upsert=True, collation=collation)

    def load_bad_kinlaw(self):
        '''
            Load kinlaw IDs whose enzymes have names as "1" or "2"
        '''
        for _id in self.bad_kinlawid:
            self.col.update_many({'$and': [{'$text': {'$search': "\"lipopolysaccharide\""}}, {'ncbi_taxonomy_id': 10116}]}, 
                {'$push': {'kinetics': {'kinlaw_id': _id, 'ncbi_taxonomy_id': 10116}}})



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
    path = os.path.join(os.path.dirname(__file__), "../../logs/protein_aggregate.txt")   
    manager = ProteinAggregate(username=username, password=password, server=server, 
                               authSource='admin', src_database=src_db,
                               verbose=True, collection=collection_str, destination_database=des_db,
                               cache_dir=path)

    # manager.copy_uniprot()
    # manager.load_abundance_from_pax()
    # manager.load_ko()
    # manager.load_taxon()
    # manager.load_unreviewed_abundance()
    # manager.load_kinlaw_from_sabio()

    # collation = Collation('en', strength=CollationStrength.SECONDARY)
    # manager.collection.create_index([("ncbi_taxonomy_id", pymongo.ASCENDING),
    #                        ("ancestor_taxon_id", pymongo.ASCENDING), ("ko_number", pymongo.ASCENDING)], 
    #                        background=True)
    # manager.collection.create_index([("uniprot_id", pymongo.ASCENDING),
    #                        ("ancestor_taxon_id", pymongo.ASCENDING)], background=True, collation=collation)
    # manager.collection.create_index([("kinetics", pymongo.ASCENDING)], background=True)

if __name__ == '__main__':
	main()