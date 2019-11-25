from datanator.util import mongo_util, file_util, chem_util
import datanator.config.core
from pymongo.collation import Collation, CollationStrength
from pymongo import ASCENDING
import os
import tempfile


class SabioCompound:

    def __init__(self, username=None, password=None, server=None, authSource='admin',
                 src_database='datanator', dest_database=None, max_entries=float('inf'), verbose=True,
                 src_collection='sabio_compound', dest_collection=None, cache_dir=None):
        '''
                Args:
                    src_database (:obj: `str`): name of database in which source collections reside
        '''
        self.mongo_manager = mongo_util.MongoUtil(MongoDB=server, username=username,
                                                  password=password, authSource=authSource, db=src_database)
        self.file_manager = file_util.FileUtil()
        self.chem_manager = chem_util.ChemUtil()
        self.collation = Collation(locale='en', strength=CollationStrength.SECONDARY)
        self.verbose = verbose
        self.max_entries = max_entries
        self.src_collection = src_collection   
    
    def add_inchi_key(self):
        """Add inchi_key field to sabio_compound collection
        in MongoDB
        """
        query = {}
        projection = {'structures._value_inchi': 1}
        _, _, collection = self.mongo_manager.con_db(self.src_collection)
        docs = collection.find(filter=query, projection=projection)
        count = collection.count_documents(query)
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            if self.verbose and i % 100 == 0:
                print('Processing doc {} out of {}'.format(i, count))
            try:
                inchi = doc['structures'][0]['_value_inchi']
            except IndexError:
                print('Compound with id {} has no structure information'.format(doc['_id']))
            except KeyError:
                print('Compound with id {} has no structure array'.format(doc['_id']))
            inchi_key = self.chem_manager.inchi_to_inchikey(inchi)
            collection.update_one({'_id': doc['_id']},
                                 {'$set': {'inchi_key': inchi_key}})

def main():
    cache_dirname = tempfile.mkdtemp()
    cache_dir = os.path.join(cache_dirname, 'logs.txt')
    src_db = 'datanator'
    collection_str = 'sabio_compound'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']        
    src = SabioCompound(username=username, password=password, server=server, 
                        authSource='admin', src_database=src_db,
                        verbose=True, src_collection=collection_str,
                        cache_dir=cache_dir)
    src.add_inchi_key()

if __name__ == '__main__':
    main()