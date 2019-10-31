from datanator_query_python.query import (query_protein, query_metabolites, 
                                         query_metabolites_meta)
from datanator_query_python.config import config as config_mongo
from karr_lab_aws_manager.elasticsearch_kl import util as es_util


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


# def main():
#     conf = config_mongo.Config()
#     username = conf.USERNAME
#     password = conf.PASSWORD
#     server = conf.SERVER
#     authDB = conf.AUTHDB
#     db = 'datanator'
#     manager = MongoToES(verbose=True, profile_name='es-poweruser', credential_path='~/.wc/third_party/aws_credentials',
#                 config_path='~/.wc/third_party/aws_config', elastic_path='~/.wc/third_party/elasticsearch.ini')
    
#     # # data from "protein" collection
#     # count, docs = manager.data_from_mongo_protein(server, db, username, password, authSource=authDB)
#     # status_code = manager.data_to_es_bulk(count, docs, 'protein', _id='uniprot_id')
#     # manager.index_settings('protein', 0) 
    
#     # # data from "ecmdb" and "ymdb" collection
#     # ecmdb_docs, ecmdb_count, ymdb_docs, ymdb_count = manager.data_from_mongo_metabolite(server, 
#     #                                                 db, username, password, authSource=authDB)
#     # status_code_0 = manager.data_to_es_bulk(ecmdb_count, ecmdb_docs, 'ecmdb', _id='m2m_id')
#     # status_code_1 = manager.data_to_es_bulk(ymdb_count, ymdb_docs, 'ymdb', _id='ymdb_id')
#     # manager.index_settings('ecmdb', 0) 
#     # manager.index_settings('ymdb', 0) 

#     # data from "metabolites_meta" collection
#     docs = manager.data_from_mongo_metabolites_meta(server, db, username, password, authSource=authDB)
#     status_code = manager.data_to_es_single(5225, docs, 'metabolites_meta', _id='InChI_Key')
#     manager.index_settings('metabolites_meta', 0)

#     print(status_code)   

# if __name__ == "__main__":
#     main()