from datanator_query_python.query import (query_protein, query_metabolites, 
                                         query_metabolites_meta, query_sabiork_old)
from datanator_query_python.config import config as config_mongo
from karr_lab_aws_manager.elasticsearch_kl import index_setting_file
from datanator.util import mongo_util
from karr_lab_aws_manager.elasticsearch_kl import util as es_util
import json
from pathlib import Path


class MongoToES(es_util.EsUtil):

    def __init__(self, profile_name=None, credential_path=None,
                config_path=None, elastic_path=None,
                cache_dir=None, service_name='es', index=None, max_entries=float('inf'), verbose=False):
        ''' Migrate data from mongodb to elasticsearch service on AWS

            Args:
                profile_name (:obj:`str`): AWS profile to use for authentication
                credential_path (:obj:`str`): directory for aws credentials file
                config_path (:obj:`str`): directory for aws config file
                elastic_path (:obj:`str`): directory for file containing aws elasticsearch service variables
                cache_dir (:obj:`str`): temp directory to store json for bulk upload
                service_name (:obj:`str`): aws service to be used
        '''
        super().__init__(profile_name=profile_name, credential_path=credential_path,
                config_path=config_path, elastic_path=elastic_path,
                cache_dir=cache_dir, service_name=service_name, max_entries=max_entries, verbose=verbose)
        self.index = index

    def data_from_mongo_protein(self, server, db, username, password, verbose=False,
                                readPreference='nearest', authSource='admin', projection={'_id': 0},
                                query={}):
        ''' Acquire documents from protein collection in datanator

            Args:
                server (:obj:`str`): mongodb ip address
                db (:obj:`str`): database name
                username (:obj:`str`): username for mongodb login
                password (:obj:`str`): password for mongodb login
                verbose (:obj:`bool`): display verbose messages
                readPreference (:obj:`str`): mongodb readpreference
                authSource (:obj:`str`): database login info is authenticating against
                projection (:obj:`str`): mongodb query projection
                query (:obj:`str`): mongodb query filter

            Returns:
                (:obj:`tuple`): tuple containing:

                    docs (:obj:`pymongo.Cursor`): pymongo cursor object that points to all documents in protein collection;
                    count (:obj:`int`): number of documents returned
        '''
        protein_manager = query_protein.QueryProtein(server=server, database=db,
                 verbose=verbose, username=username, authSource=authSource,
                 password=password, readPreference=readPreference)
        docs = protein_manager.collection.find(filter=query, projection=projection)
        count = protein_manager.collection.count_documents(query)
        return (count, docs)

    def data_from_mongo_metabolite(self, server, db, username, password, verbose=False,
                                readPreference='nearest', authSource='admin', projection={'_id': 0},
                                query={}):
        ''' Acquire documents from metabolite (ecmdb/ymdb) collection in datanator

            Args:
                server (:obj:`str`): mongodb ip address
                db (:obj:`str`): database name
                username (:obj:`str`): username for mongodb login
                password (:obj:`str`): password for mongodb login
                verbose (:obj:`bool`): display verbose messages
                readPreference (:obj:`str`): mongodb readpreference
                authSource (:obj:`str`): database login info is authenticating against
                projection (:obj:`str`): mongodb query projection
                query (:obj:`str`): mongodb query filter

            Returns:
                (:obj:`tuple`): tuple containing:

                    ecmdb_docs (:obj:`pymongo.Cursor`): pymongo cursor object that points to all documents in ecmdb collection;
                    ecmdb_count (:obj:`int`): number of documents returned in ecmdb;
                    ymdb_docs (:obj:`pymongo.Cursor`): pymongo cursor object that points to all documents in ymdb collection;
                    ymdb_count (:obj:`int`): number of documents returned in ymdb
        '''
        metabolite_manager = query_metabolites.QueryMetabolites( 
                 MongoDB=server, db=db,
                 verbose=verbose, username=username,
                 password=password, readPreference=readPreference,
                 authSource=authSource)
        ecmdb_docs = metabolite_manager.collection_ecmdb.find(filter=query,projection=projection)
        ecmdb_count = metabolite_manager.collection_ecmdb.count_documents(query)
        ymdb_docs = metabolite_manager.collection_ymdb.find(filter=query, projection=projection)
        ymdb_count = metabolite_manager.collection_ymdb.count_documents(query)

        return ecmdb_docs, ecmdb_count, ymdb_docs, ymdb_count

    def data_from_mongo_metabolites_meta(self, server, db, username, password, verbose=False,
                                readPreference='nearest', authSource='admin', projection={'_id': 0},
                                query={}):
        ''' Acquire documents from metabolites_meta collection in datanator

            Args:
                server (:obj:`str`): mongodb ip address
                db (:obj:`str`): database name
                username (:obj:`str`): username for mongodb login
                password (:obj:`str`): password for mongodb login
                verbose (:obj:`bool`): display verbose messages
                readPreference (:obj:`str`): mongodb readpreference
                authSource (:obj:`str`): database login info is authenticating against
                projection (:obj:`str`): mongodb query projection
                query (:obj:`str`): mongodb query filter

            Returns:
                (:obj:`tuple`): tuple containing:

                    docs (:obj:`pymongo.Cursor`): pymongo cursor object that points to all documents in protein collection;
                    count (:obj:`int`): number of documents returned
        '''
        manager = query_metabolites_meta.QueryMetabolitesMeta(MongoDB=server, db=db,
                 collection_str='metabolites_meta', verbose=verbose, username=username,
                 password=password, authSource=authSource, readPreference=readPreference)
        docs = manager.collection.find(filter=query, projection=projection)
        for doc in docs:
            if doc['InChI_Key'] is None:
                continue    
            sim_com = doc['similar_compounds']
            tmp = []
            for compound in sim_com:
                inchi_key = list(compound.keys())[0]
                score = list(compound.values())[0]
                new_dic = {'inchikey': inchi_key, 'similarity_score': score}
                tmp.append(new_dic)
            doc['similar_compounds'] = tmp
            yield doc     

    def data_from_mongo_sabiork(self, server, db, username, password, verbose=False,
                                readPreference='nearest', authSource='admin', projection={'_id': 0},
                                query={}):
        ''' Acquire documents from protein collection in datanator

            Args:
                server (:obj:`str`): mongodb ip address
                db (:obj:`str`): database name
                username (:obj:`str`): username for mongodb login
                password (:obj:`str`): password for mongodb login
                verbose (:obj:`bool`): display verbose messages
                readPreference (:obj:`str`): mongodb readpreference
                authSource (:obj:`str`): database login info is authenticating against
                projection (:obj:`str`): mongodb query projection
                query (:obj:`str`): mongodb query filter

            Returns:
                (:obj:`tuple`): tuple containing:

                    docs (:obj:`pymongo.Cursor`): pymongo cursor object that points to all documents in protein collection;
                    count (:obj:`int`): number of documents returned
        '''
        sabio_manager = query_sabiork_old.QuerySabioOld(MongoDB=server, db=db,
                 verbose=verbose, username=username, authSource=authSource,
                 password=password, readPreference=readPreference)
        docs = sabio_manager.collection.find(filter=query, projection=projection)
        count = sabio_manager.collection.count_documents(query)
        return (count, docs)

    def data_from_mongo_sabiork_rxn_entries(self, server, db, username, password, verbose=False,
                                readPreference='nearest', authSource='admin', projection={'_id': 0},
                                query={}):
        ''' Acquire documents from protein collection in datanator

            Args:
                server (:obj:`str`): mongodb ip address
                db (:obj:`str`): database name
                username (:obj:`str`): username for mongodb login
                password (:obj:`str`): password for mongodb login
                verbose (:obj:`bool`): display verbose messages
                readPreference (:obj:`str`): mongodb readpreference
                authSource (:obj:`str`): database login info is authenticating against
                projection (:obj:`str`): mongodb query projection
                query (:obj:`str`): mongodb query filter

            Returns:
                (:obj:`tuple`): tuple containing:

                    docs (:obj:`pymongo.Cursor`): pymongo cursor object that points to all documents in protein collection;
                    count (:obj:`int`): number of documents returned
        '''
        mongo_manager = mongo_util.MongoUtil(MongoDB=server, username=username,
                                            password=password, authSource=authSource, db=db)
        _, _, collection = mongo_manager.con_db('sabio_reaction_entries')
        docs = collection.find(filter=query, projection=projection)
        count = collection.count_documents(query)
        return (count, docs)

    def data_from_mongo_rna_halflife_entries(self, server, db, username, password, verbose=False,
                                readPreference='nearest', authSource='admin',
                                query={}, collection_str='rna_halflife'):
        ''' Acquire documents from protein collection in datanator

            Args:
                server (:obj:`str`): mongodb ip address
                db (:obj:`str`): database name
                username (:obj:`str`): username for mongodb login
                password (:obj:`str`): password for mongodb login
                verbose (:obj:`bool`): display verbose messages
                readPreference (:obj:`str`): mongodb readpreference
                authSource (:obj:`str`): database login info is authenticating against
                projection (:obj:`str`): mongodb query projection
                query (:obj:`str`): mongodb query filter

            Returns:
                (:obj:`tuple`): tuple containing:

                    docs (:obj:`pymongo.Cursor`): pymongo cursor object that points to all documents in protein collection;
                    count (:obj:`int`): number of documents returned
        '''
        mongo_manager = mongo_util.MongoUtil(MongoDB=server, username=username,
                                            password=password, authSource=authSource, db=db)
        _, _, collection = mongo_manager.con_db(collection_str)
        docs = collection.find(filter=query)
        count = collection.count_documents(query)
        return (count, docs)


def main():
    conf = config_mongo.Config()
    username = conf.USERNAME
    password = conf.PASSWORD
    server = conf.SERVER
    authDB = conf.AUTHDB
    db = 'datanator'
    manager = MongoToES(verbose=True, profile_name='es-poweruser', credential_path='~/.wc/third_party/aws_credentials',
                config_path='~/.wc/third_party/aws_config', elastic_path='~/.wc/third_party/elasticsearch.ini')

    filter_dir = '/root/host/karr_lab/karr_lab_aws_manager/karr_lab_aws_manager/elasticsearch_kl/filters/autocomplete_filter.json'
    analyzer_dir = '/root/host/karr_lab/karr_lab_aws_manager/karr_lab_aws_manager/elasticsearch_kl/analyzers/auto_complete.json'
    
    # old_index = 'protein'
    # new_index = 'protein_something'
    # _, _, _, = manager.migrate_index(old_index, new_index)
    # _, _, _, = manager.migrate_index(new_index, old_index)


    # # data from "protein" collection
    # index_name = 'protein'
    # _ = manager.delete_index(index_name)
    # count, docs = manager.data_from_mongo_protein(server, db, username, password, authSource=authDB)
    # mappings_dir = '/root/host/karr_lab/karr_lab_aws_manager/karr_lab_aws_manager/elasticsearch_kl/mappings/protein.json'
    # index_manager = index_setting_file.IndexUtil(filter_dir=filter_dir, analyzer_dir=analyzer_dir, mapping_properties_dir=mappings_dir)
    # setting_file = index_manager.combine_files(_filter=True, analyzer=True, mappings=True)
    # _ = manager.create_index_with_file(index_name, setting_file)
    # _ = manager.data_to_es_bulk(docs, count=count, index=index_name, _id='uniprot_id')
    
    # data from "ecmdb" and "ymdb" collection
    # ecmdb_docs, ecmdb_count, ymdb_docs, ymdb_count = manager.data_from_mongo_metabolite(server, 
    #                                                 db, username, password, authSource=authDB)
    # ecmdb = 'ecmdb'
    # _ = manager.delete_index(ecmdb)
    # ecmdb_mappings_dir = '/root/host/karr_lab/karr_lab_aws_manager/karr_lab_aws_manager/elasticsearch_kl/mappings/ecmdb.json'
    # ecmdb_index_manager = index_setting_file.IndexUtil(filter_dir=filter_dir, analyzer_dir=analyzer_dir, mapping_properties_dir=ecmdb_mappings_dir)
    # ecmdb_setting_file = ecmdb_index_manager.combine_files(_filter=True, analyzer=True, mappings=True)
    # _ = manager.create_index_with_file(ecmdb, ecmdb_setting_file)
    # status = manager.data_to_es_bulk(ecmdb_docs, index=ecmdb, count=ecmdb_count, _id='m2m_id')

    # ymdb = 'ymdb'
    # _ = manager.delete_index(ymdb)
    # ymdb_mappings_dir = '/root/host/karr_lab/karr_lab_aws_manager/karr_lab_aws_manager/elasticsearch_kl/mappings/ymdb.json'
    # ymdb_index_manager = index_setting_file.IndexUtil(filter_dir=filter_dir, analyzer_dir=analyzer_dir, mapping_properties_dir=ymdb_mappings_dir)
    # ymdb_setting_file = ymdb_index_manager.combine_files(_filter=True, analyzer=True, mappings=True)
    # _ = manager.create_index_with_file(ymdb, ymdb_setting_file)    
    # status = manager.data_to_es_bulk(ymdb_docs, index=ymdb, count=ymdb_count, _id='ymdb_id')

    # # data from "metabolites_meta" collection
    # index_name = 'metabolites_meta'
    # _ = manager.create_index(index_name)
    # docs = manager.data_from_mongo_metabolites_meta(server, db, username, password, authSource=authDB)
    # mappings_dir = '/root/host/karr_lab/karr_lab_aws_manager/karr_lab_aws_manager/elasticsearch_kl/mappings/metabolites_meta.json'
    # index_manager = index_setting_file.IndexUtil(filter_dir=filter_dir, analyzer_dir=analyzer_dir, mapping_properties_dir=mappings_dir)
    # setting_file = index_manager.combine_files(_filter=True, analyzer=True, mappings=True)
    # _ = manager.create_index_with_file(index_name, setting_file)
    # _ = manager.data_to_es_single(5225, docs, index_name, _id='InChI_Key')

    # data from "sabio_rk_old" collection
    # count, docs = manager.data_from_mongo_sabiork(server, db, username, password, authSource=authDB)
    # index_name = 'sabio_rk'
    # mappings_dir = '/root/host/karr_lab/karr_lab_aws_manager/karr_lab_aws_manager/elasticsearch_kl/mappings/sabio_rk.json'
    # index_manager = index_setting_file.IndexUtil(filter_dir=filter_dir, analyzer_dir=analyzer_dir, mapping_properties_dir=mappings_dir)
    # setting_file = index_manager.combine_files(_filter=True, analyzer=True, mappings=True)
    # _ = manager.create_index_with_file(index_name, setting_file)
    # _ = manager.data_to_es_bulk(docs, index=index_name, count=count, _id='kinlaw_id')

    # # data from "sabio_reaction_entries" collection
    # index_name = 'sabio_reaction_entries'
    # _ = manager.delete_index(index_name)
    # count, docs = manager.data_from_mongo_sabiork_rxn_entries(server, db, username, password, authSource=authDB)
    # r = manager.create_index(index_name)
    # _ = manager.data_to_es_bulk(docs, index=index_name, count=count, _id='rxn_id')

    # data from "rna_halflife" collection
    # count, docs = manager.data_from_mongo_rna_halflife_entries(server, db, username, password, authSource=authDB)
    # index_schema_path = str(Path('~/host/karr_lab/datanator/datanator/data_source/schema/rna_halflife.json').expanduser())
    # with open(index_schema_path) as json_file:
    #     index_schema = json.load(json_file)
    # _ = manager.create_index('rna_halflife', mappings=index_schema)
    # status = manager.data_to_es_bulk(docs, index='rna_halflife', count=count, _id='_id')

    # data from "taxon_tree" collection
    # index_name = 'taxon_tree'
    # count, docs = manager.data_from_mongo_rna_halflife_entries(server, db, username, password, authSource=authDB, collection_str=index_name)
    # mappings_dir = '/root/host/karr_lab/karr_lab_aws_manager/karr_lab_aws_manager/elasticsearch_kl/mappings/taxon_tree.json'
    # index_manager = index_setting_file.IndexUtil(filter_dir=filter_dir, analyzer_dir=analyzer_dir, mapping_properties_dir=mappings_dir)
    # setting_file = index_manager.combine_files(_filter=True, analyzer=True, mappings=True)
    # _ = manager.create_index_with_file(index_name, setting_file)
    # _ = manager.data_to_es_bulk(docs, index=index_name, count=count, _id='tax_id')

    # r = manager.index_health_status()
    # print(r.content.decode('utf-8'))
    # print('Operation status: {}'.format(status.content))

if __name__ == "__main__":
    main()