from datanator.core import query_nosql
import re

class MetabolitesMeta(query_nosql.DataQuery):

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf'), username = None, password = None,
                 authSource = 'admin'):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.replicaSet = replicaSet
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries

        super(MetabolitesMeta, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet,
                                              db=db, verbose=verbose, max_entries=max_entries, username = username,
                                              password = password, authSource = 'admin')
        self.frequency = 100

    def load_content(self):
        ecmdb_fields = ['m2m_id', 'inchi', 'synonyms.synonym']
        ecmdb_list = self.get_metabolite_fields(
            fields=ecmdb_fields, collection_str='ecmdb')

        ymdb_fields = ['ymdb_id', 'inchi', 'synonyms.synonym']
        ymdb_list = self.get_metabolite_fields(
            fields=ymdb_fields, collection_str='ymdb')

        collection_name = 'metabolites_meta'
        client, _, collection = self.con_db(collection_name)

        for doc in ecmdb_list:
            collection.update_one({'inchi': doc['inchi']},
                                  { '$set': doc},
                                  upsert=True)

        for doc in ymdb_list:
            collection.update_one({'inchi': doc['inchi']},
                                  { '$set': doc},
                                  upsert=True)

        # for doc in self.doc_feeder(collection_str=collection_name, query={}, projection={'inchi'}):
        #     kinlaw_id = self.find_rxn_id(inchi = doc['inchi'])
        #     rxn_participants = self.find_reaction_participants(kinlaw_id)
        #     collection.update_one({'inchi': doc['inchi']},
        #                           {'$set': {'kinlaw_id': kinlaw_id,
        #                            'reaction_participants': rxn_participants}},
        #                           upsert=False)
        

        client.close()


    def get_metabolite_fields(self, fields=None, collection_str=None):
        '''Get fields of interest from metabolite collection: ecmdb or ymdb
                Args:
                        fileds: list of fields of interest
                        collection_str: collection in which query will be done
                Return:
                        list of dictionaries
                                [{'filed1': value1, 'field2': value2, ... }, ..., {}]
        '''
        projection = {}
        for field in fields:
            projection[field] = 1
        projection['_id'] = 0
        cursor = self.doc_feeder(collection_str=collection_str, query={},
                                 projection=projection)
        i = 0
        dict_list = []
        for doc in cursor:
            if i == self.max_entries:
                break
            if i % self.frequency == 0:
                print('Getting fields of interest from {}'.format(collection_str))
            doc['inchi'] = self.parse_inchi(inchi=doc['inchi'])
            dict_list.append(doc)
            i += 1

        return dict_list

    def find_rxn_id(self, inchi=None):
        '''Find reactions' kinlaw_id in sabio_rk given inchi structures
                Return:
                        List of kinlaw_id that contain the given inchi string
        '''
        substrate = 'reaction_participant.substrate.structure.inchi'
        product = 'reaction_participant.product.structure.inchi'
        c = 'sabio_rk'
        # regular expressions are weird
        try:
            inchi = inchi.replace('(', '\(').replace(')', '\)')
        except AttributeError:
            return -1
        regex = re.compile(inchi)
        query = {'$or': [{substrate: regex},
                         {product: regex}]}
        docs = self.doc_feeder(collection_str=c, query=query)
        kinlaw_id = []
        i = 0
        for doc in docs:
            if i == self.max_entries:
                break
            if i % self.frequency == 0:
                print('Getting kinlaw_id of inchi {}'.format(inchi))
            _id = doc['kinlaw_id']
            kinlaw_id.append(_id)
            i += 1

        return kinlaw_id

    def find_metabolite_inchi(self, doc):
        '''Find inchi structure information of metabolites in ecmdb or ymdb
                return inchi information without protonation state
        '''
        try:
            return self.parse_inchi(doc['inchi'])
        except KeyError:
            return 'No key named "inchi" in given document'

    def parse_inchi(self, inchi=None):
        '''Remove molecules's protonation state
        "InChI=1S/H2O/h1H2" = > "InChI=1S/H2O"
        '''
        # if self.verbose:
        #     print('Parsing inchi by taking out protonation state')
        try:
            inchi_neutral = inchi.split('/h')[0]
            return inchi_neutral
        except AttributeError:
            return None


def main():
    MongoDB = '35.173.159.185:27017'
    db = 'datanator'
    username = None
    password = None
    manager = MetabolitesMeta(cache_dirname=None, MongoDB=MongoDB, replicaSet = None, db=db, 
                                verbose=True, max_entries=float('inf'), 
                                username = username, password = password)

    manager.load_content()


if __name__ == '__main__':
    main()
