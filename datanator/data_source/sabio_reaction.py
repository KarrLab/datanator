from datanator_query_python.util import mongo_util, file_util
from datanator_query_python.query import query_sabiork_old, query_metabolites_meta
import datanator.config.core
import hashlib
from pymongo.collation import Collation, CollationStrength
from pymongo import ASCENDING
import os
import tempfile


class RxnAggregate(mongo_util.MongoUtil):

    def __init__(self, username=None, password=None, server=None, authSource='admin',
                 src_database='datanator', max_entries=float('inf'), verbose=True,
                 collection='sabio_reaction_entries', destination_database='datanator', cache_dir=None):
        '''
                Args:
                        src_database (:obj: `str`): name of database in which source collections reside
                        destination_database (:obj: `str`): name of database to put the aggregated collection
        '''
        super().__init__(MongoDB=server, db=destination_database, username=username, password=password,
                        authSource=authSource)
        self.col = self.db_obj[collection]
        self.query_manager = query_sabiork_old.QuerySabioOld(MongoDB=server, password=password, authSource=authSource, username=username)
        self.metabolites_meta_manager = self.client[src_database]['metabolites_meta']
        self.file_manager = file_util.FileUtil()
        self.collation = Collation(locale='en', strength=CollationStrength.SECONDARY)
        self.verbose = verbose
        self.max_entries = max_entries

    def fill_collection(self):
        projection = {'_id': 0,'resource': 1, 'reaction_participant': 1,
                    'reaction_participant': 1, 'kinlaw_id': 1, 'enzymes': 1}
        _, _, collection = self.db_obj['sabio_rk_old']
        docs = collection.find({}, projection=projection)
        count = collection.count_documents({})
        start = 0
        for i, doc in enumerate(docs[start:]):
            if self.verbose and i % 100 == 0:
                print('Processing document {} out of {}'.format(i+start, count))
            if doc.get('resource') is None:
                continue
            if i == self.max_entries:
                break
            kinlaw_id = doc['kinlaw_id']
            _, have = self.query_manager.get_rxn_with_prm([kinlaw_id])
            if len(have) == 1:
                with_prm = True
            else:
                with_prm = False
            key = 'has_poi.'+str(kinlaw_id)
            rxn_id = self.get_rxn_id(doc)
            reactants = self.create_reactants(doc)
            substrate_names, product_names = self.extract_reactant_names(doc)
            enzyme_names = self.extract_enzyme_names(doc)
            ec = self.get_ec(doc)
            self.col.update_one({'rxn_id': rxn_id},
                                {'$addToSet': {'kinlaw_id': kinlaw_id},
                                '$set': {'substrates': reactants['substrate_aggregate'],
                                        'products': reactants['product_aggregate'],
                                        'substrate_names': substrate_names,
                                        'product_names': product_names,
                                        'enzyme_names': enzyme_names,
                                        'ec-code': ec,
                                        key: with_prm}}, upsert=True)
            if i == 0:
                self.col.create_index([("rxn_id", ASCENDING)], background=True)

    def get_rxn_id(self, doc):
        resource = doc['resource']
        sr = self.file_manager.search_dict_list(resource, 'namespace', 'sabiork.reaction')
        _id = sr[0]['id']
        return int(_id)

    def get_ec(self, doc):
        resource = doc['resource']
        sr = self.file_manager.search_dict_list(resource, 'namespace', 'ec-code')
        _id = sr[0]['id']
        return _id
    
    def create_reactants(self, doc):
        result = {}
        substrate_aggregate = doc['reaction_participant'][3]['substrate_aggregate']
        product_aggregate = doc['reaction_participant'][4]['product_aggregate']
        result['substrate_aggregate'] = substrate_aggregate
        result['product_aggregate'] = product_aggregate

        return result

    def extract_reactant_names(self, doc):
        """Extract compound information from doc dictionary
        
        Args:
            doc (:obj:`dict`): sabio_rk_old document

        Returns:
            (:obj:`tuple`): substrates and products names [[],[],...,[]], [[],[],...,[]] 
        """
        substrates = doc['reaction_participant'][0]['substrate']
        products = doc['reaction_participant'][1]['product']

        def extract_names(compound, side='substrate'):
            """Extract names of compound
            
            Args:
                compound (:obj:`dict`): compound information
                side (:obj:`str`, optional): substrate or product. Defaults to 'substrate'.
            
            Returns:
                (:obj:`list`): list of names for the compound
            """
            name = compound[side+'_name']
            syn = compound[side+'_synonym']
            syn.append(name)
            return syn

        def iter_compound(compounds, side='substrate'):
            """Iterate through compound list
            
            Args:
                compounds (:obj:`list`): a list of compounds
                side (:obj:`str`, optional): substrate or product. Defaults to 'substrate'.
            """
            result = []
            for compound in compounds:
                names = extract_names(compound, side=side)
                result.append(names)
            return result

        substrate_names = iter_compound(substrates, side='substrate')
        product_names = iter_compound(products, side='product')

        return substrate_names, product_names

    def extract_enzyme_names(self, doc):
        """Extract enzyme names
        
        Args:
            doc (:obj:`dict`): sabio_rk_old document

        Returns:
            (:obj:`list`): list of enzyme names
        """
        result = []
        enzymes = doc['enzymes'][0]['enzyme']
        for enzyme in enzymes:
            enzyme_name = enzyme['enzyme_name']
            syn = enzyme['enzyme_synonym']
            if syn is None:
                syn = enzyme_name
            else:
                syn.append(enzyme_name)
            result.append(syn)
        if len(enzymes) == 1:
            return syn
        else:
            return result
        
    def label_existence(self, start=0):
        """Label reactant's existence in metabolites collections.
        """
        docs = self.metabolites_meta_manager.find({})
        count = self.metabolites_meta_manager.count_documents({})
        for i, doc in enumerate(docs[start:]):
            inchi_key = doc.get('InChI_Key')
            if inchi_key is None:
                continue
            if i == self.max_entries:
                break
            if i % 50 == 0 and self.verbose:
                print("Process metabolite {} out of {} ...".format(i+start, count))
            con_0 = {'products': inchi_key}
            con_1 = {'substrates': inchi_key}
            query = {'$or': [con_0, con_1]}
            self.col.update_many(query,
                                 {'$addToSet': {'in_metabolites': inchi_key}}, upsert=False,
                                 collation=self.collation)

    def hash_null_reactants(self, start=0):
        """(https://github.com/KarrLab/datanator/issues/50)
           (https://github.com/KarrLab/datanator_rest_api/issues/116)

        Args:
            start (:obj:`int`, optional): Start of document. Defaults to 0.
        """
        con_0 = {'substrates': None}
        con_1 = {'products': None}
        con_2 = {'substrates': []}
        con_3 = {'products': []}
        query = {'$or': [con_0, con_1, con_2, con_3]}
        docs = self.col.find(query, skip=start, no_cursor_timeout=True, batch_size=10)
        count = self.col.count_documents(query)
        sabio_rk_old = self.db_obj['sabio_rk_old']
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break            
            if i % 10 == 0 and self.verbose:
                print("Process entry {} out of {} ...".format(i+start, count))
            products_hashed = doc['products']
            substrates_hashed = doc['substrates']

            if products_hashed == []:
                products_hashed = [None]
            if substrates_hashed == []:
                substrates_hashed = [None]        

            products_name = [x[0] for x in doc['product_names']]
            substrates_name = [x[0] for x in doc['substrate_names']]
            products_name_hashed = [hashlib.new('ripemd160', x.encode()).hexdigest() for x in products_name]
            substrates_name_hashed = [hashlib.new('ripemd160', x.encode()).hexdigest() for x in substrates_name]

            products_hashed = [x if x is not None else products_name_hashed[i] for i, x in enumerate(products_hashed)]
            substrates_hashed = [x if x is not None else substrates_name_hashed[i] for i, x in enumerate(substrates_hashed)]

            sabio_rk_old.update_many({"resource": {"namespace": "sabiork.reaction", "id": str(doc['rxn_id'])}},
                                     {"$set": {"reaction_participant.3.substrate_aggregate": substrates_hashed,
                                               "reaction_participant.4.product_aggregate": products_hashed}},
                                     upsert=False)

            self.col.update_one({'_id': doc['_id']},
                                {'$set': {'products': products_hashed,
                                          'substrates': substrates_hashed}},
                                upsert=False)



def main():
    cache_dirname = tempfile.mkdtemp()
    cache_dir = os.path.join(cache_dirname, 'logs.txt')
    src_db = 'datanator'
    des_db = 'datanator'
    collection_str = 'sabio_reaction_entries'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']      
    src = RxnAggregate(username=username, password=password, server=server, 
                        authSource='admin', src_database=src_db,
                        verbose=True, collection=collection_str, destination_database=des_db,
                        cache_dir=cache_dir)
    # src.fill_collection()
    # src.label_existence()

    src.hash_null_reactants()

if __name__ == '__main__':
    main()