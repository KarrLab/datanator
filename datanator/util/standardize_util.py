'''standardize key to a uniform nomenclature
'''
from datanator.util import mongo_util
from bson.objectid import ObjectId


class StandardizeUtil(mongo_util.MongoUtil):

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf')):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        super(StandardizeUtil, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB,replicaSet=replicaSet, db=db,
                                              verbose=verbose, max_entries=max_entries)

    '''Process sabio_rk documents
    '''

    def standardize_sabio(self):
        client, _, collection = self.con_db('sabio_rk')

        # update_many has too many limitations
        i = 0
        for doc in collection.find().limit(self.max_entries):
            if self.verbose and i % 100 == 0:
                print('processing {} out of {} documents'.format(
                    i, min(self.max_entries, collection.count())))
            if 'id' in doc['resource'][0]:
                sorted_resource_list = sorted(
                    doc['resource'], key=lambda i: i['namespace'], reverse=True)
                id_list = [d['id'] for d in sorted_resource_list]
                if len(id_list) == 3:  # no kegg.reaction:
                    doc['resource'][0] = {'sabiork_reation_id': id_list[0]}
                    doc['resource'][1] = {'pubmed_id': id_list[1]}
                    doc['resource'][2] = {'ec_code_id': id_list[2]}
                else:  # with kegg.reaction
                    doc['resource'][0] = {'sabiork_reation_id': id_list[0]}
                    doc['resource'][1] = {'pubmed_id': id_list[1]}
                    doc['resource'][2] = {'kegg_id': id_list[2]}
                    doc['resource'][3] = {'ec_code_id': id_list[3]}

            if 'sequence' in doc['enzymes'][2]['subunit'][0]:
                for unit in doc['enzymes'][2]['subunit']:
                    unit['canonical_sequence'] = unit['sequence']

            collection.save(doc)
            i += 1
        client.close()

    '''Process metabolite documents
    '''

    def standardize_metabolite(self):
        client_ecmdb, _, collection_ecmdb = self.con_db('ecmdb')
        collection_ecmdb.update_many({},
                                     {'$rename': {'smiles': 'smiles_structure',
                                                  'inchi': 'inchi_structure',
                                                  'inchikey': 'inchi_key',

                                                  }
                                      }
                                     )
        client_ecmdb.close()

        client_ymdb, _, collection_ymdb = self.con_db('ymdb')
        collection_ymdb.update_many({},
                                    {'$rename': {'smiles': 'smiles_structure',
                                                 'inchi': 'inchi_structure',
                                                 'inchikey': 'inchi_key',

                                                 }
                                     }
                                    )
        client_ymdb.close()
