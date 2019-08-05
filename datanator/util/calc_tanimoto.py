import pybel
from datanator.util import mongo_util
import pymongo
import numpy as np
import multiprocessing as mp
import datanator.config.core


class CalcTanimoto(mongo_util.MongoUtil):
    '''Calculating the tanitomo similarity matrix
            given two compound collections e.g.
            ECMDB YMDB
    '''

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf'), username=None,
                 password=None, authSource='admin'):
        self.authSource = authSource
        self.username = username
        self.password = password
        self.replicaSet = replicaSet
        self.db = db
        self.MongoDB = MongoDB
        super(CalcTanimoto, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet,
                                           db=db, verbose=verbose, max_entries=max_entries, username=username,
                                           password=password, authSource=authSource)
        log_handler = pybel.ob.OBMessageHandler()
        log_handler.SetOutputLevel(0) 
        pybel.ob.obErrorLog.SetOutputLevel(0)

    def get_tanimoto(self, mol1, mol2, str_format='inchi', rounding=3):
        '''Calculates tanimoto coefficients between
                two molecules, mol1 and mol2
                Args:
                        mol1: molecule 1 in some format
                        mol2: molecule 2 in same format as molecule 1
                        str_format: format for molecular representation
                                                supported formats are provided by Pybel
                        rounding: rounding of the final results
                Return:
                        tani: rounded tanimoto coefficient
        '''
        try:           
            inchi = [mol1, mol2]
            mols = [pybel.readstring(str_format, x) for x in inchi]
            fps = [x.calcfp() for x in mols]
            return round((fps[0] | fps[1]), rounding)
        except TypeError:
            return -1

    def one_to_many(self, inchi, collection_str='metabolites_meta',
                    field='inchi', lookup='InChI_Key', num=100):
        ''' Calculate tanimoto coefficients between one
                metabolite and the rest of the 'collection_str'
                Args:
                    inchi: chosen chemical compound in InChI format
                    collection_str: collection in which comparisons are made
                    field: field that has the chemical structure
                    lookup: field that had been previous indexed
                    num: max number of compounds to be returned, sorted by tanimoto
                Returns:
                        sorted_coeff: sorted numpy array of top num tanimoto coeff
                        sorted_inchi: sorted top num inchi
        '''
        _, _, col = self.con_db(collection_str)
        coeff_np = np.empty([0])
        top_inchi = []

        np_size = 0
        projection = {field: 1, lookup: 1}
        cursor = col.find({}, projection=projection)

        count = col.count_documents({})
        total = min(count, self.max_entries)

        i = 0
        while (np_size < num):  # fill in first num tanimoto coefficients
            mol2 = cursor[i][field]
            hash2 = cursor[i][lookup]
            tanimoto = self.get_tanimoto(inchi, mol2)
            if tanimoto < 1:
                coeff_np = np.append(coeff_np, tanimoto)
                top_inchi.append(hash2)
                np_size += 1
                i += 1
            else:
                i +=1

        coeff_min = np.amin(coeff_np)
        min_index = np.argmin(coeff_np)

        i = 0
        j = 0
        for doc in cursor[num:]:  # iterate through the rest of the documents
            if i > self.max_entries:
                break
            if self.verbose and j % 200 == 0:
                print('     Calculating between given and doc {} out of {} in collection {}'.format(
                    j + num, total, collection_str))
            mol2 = doc[field]
            hash2 = doc[lookup]
            tanimoto = self.get_tanimoto(inchi, mol2)
            if tanimoto > coeff_min and tanimoto < 1:
                np.put(coeff_np, min_index, tanimoto)
                top_inchi[min_index] = hash2
                # update min coeff information
                coeff_min = np.amin(coeff_np)
                min_index = np.argmin(coeff_np)
                i += 1
                j += 1
            else:
                j += 1

        indices = np.argsort(coeff_np)
        sorted_inchi = []
        for x in (indices[::-1]):
            sorted_inchi.append(top_inchi[x])
        sorted_coeff = np.sort(coeff_np)[::-1]

        return sorted_coeff, sorted_inchi

    def many_to_many(self, collection_str1='metabolites_meta',
                     collection_str2='metabolites_meta', field1='inchi',
                     field2='inchi', lookup1='InChI_Key',
                     lookup2='InChI_Key', num=100):
        ''' Go through collection_str and assign each
                compound top 'num' amount of most similar 
                compounds
                Args:
                        collection_str1: collection in which compound is drawn
                        collection_str2: collection in which comparison is made
                        field1: field of interest in collection_str1
                        field2: filed of interest in collection_str2
                        num: number of most similar compound
                        batch_size: batch_size for each server round trip
        '''
        client = pymongo.MongoClient(
            self.MongoDB, replicaSet=self.replicaSet,
            username=self.username, password=self.password,
            authSource=self.authSource)
        db_obj = client[self.db]
        final = db_obj[collection_str1]

        projection = {'m2m_id':0,  'ymdb_id': 0, 'kinlaw_id': 0, 
                    'reaction_participants': 0, 'synonyms': 0}
        _, _, col = self.con_db(collection_str1)
        count = col.count_documents({})
        total = min(count, self.max_entries)

        ''' The rest of the code in this function is to force
            a cursor refresh every 'limit' number of documents
            because no_cursor_timeout option in pymongo's find()
            function is not working as intended
        '''
        def process_doc(doc, final, i, total = total, collection_str1 = collection_str1,
                        field1 = field1, lookup1 = lookup1, collection_str2 = collection_str2,
                        field2 = field2, lookup2 = lookup2):
            # if 'similar_compounds_corrected' in doc:
            #     if self.verbose and i % 10 ==0:
            #         print('Skipping document {} out of {} in collection {}'.format(
            #             i, total, collection_str1))
            #     return 
            if i > self.max_entries:
                return 
            if self.verbose and i % 1 == 0:
                print('Going through document {} out of {} in collection {}'.format(
                    i, total, collection_str1))
                print(doc[field1])
            compound = doc[field1]
            coeff, inchi_hashed = self.one_to_many(compound, lookup=lookup2,
                                                   collection_str=collection_str2, field=field2, num=num)
            result = []
            for a, b in zip(coeff, inchi_hashed):
                dic = {}
                dic[b] = a
                result.append(dic)

            final.update_one({lookup1: doc[lookup1]},
                             {'$set': {'similar_compounds_corrected': result}},
                             upsert=False)
 
        limit = 100    # number of documents from the cursor to be stuffed into a list
        sorted_field = lookup1 # indexed field used to sort cursor
        i = 0

        documents = list(col.find({}, projection = projection).sort(sorted_field, pymongo.ASCENDING).limit(limit))
        for doc in documents: 
            process_doc(doc, final, i)
            i += 1

        is_last_batch = False
        while not is_last_batch:
            cursor = col.find({sorted_field: {'$gt': documents[-1][sorted_field]}}, projection = projection)
            documents = list(cursor.sort(sorted_field, pymongo.ASCENDING).limit(limit))
            is_last_batch = False if len(documents) == limit else True 
            for doc in documents:
                process_doc(doc, final, i)
                i += 1

            
def main():

    db = 'datanator'
    username = datanator.config.core.get_config()['datanator']['mongodb']['user']
    password = datanator.config.core.get_config()['datanator']['mongodb']['password']
    server = datanator.config.core.get_config()['datanator']['mongodb']['server']
    port = datanator.config.core.get_config()['datanator']['mongodb']['port']
    replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
    manager = CalcTanimoto(
        MongoDB=server, replicaSet=replSet, db=db,
        verbose=True, password=password, username=username)
    manager.many_to_many(field1 = 'inchi', field2 = 'inchi')


if __name__ == '__main__':
    main()
