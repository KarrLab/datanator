import pybel
from datanator.util import mongo_util
from datanator.util import server_util
import pymongo
import numpy as np


class CalcTanimoto(mongo_util.MongoUtil):
    '''Calculating the tanitomo similarity matrix
            given two compound collections e.g.
            ECMDB YMDB
    '''

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf'), username=None, result_db = None,
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
        if mol1 == 'InChI = None' or mol2 == 'InChI = None':
        	return -1
        else:
	        inchi = [mol1, mol2]
	        mols = [pybel.readstring(str_format, x) for x in inchi]
	        fps = [x.calcfp() for x in mols]
	        return round((fps[0] | fps[1]), rounding)

    def one_to_many(self, inchi, collection_str = 'metabolites_meta', 
    				field = 'inchi_deprot', lookup = 'inchi_hashed', num = 100):
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

        while (np_size <= num):  # iterate through first 100 documents
            mol2 = cursor[np_size][field]
            hash2 = cursor[np_size][lookup]
            tanimoto = self.get_tanimoto(inchi, mol2)
            coeff_np = np.append(coeff_np, tanimoto)
            top_inchi.append(hash2)
            np_size += 1

        coeff_min = np.amin(coeff_np)
        min_index = np.argmin(coeff_np)

        i = 0
        for doc in cursor[num:]:  # iterate through the rest of the documents
            if i > self.max_entries:
                break
            mol2 = doc[field]
            hash2 = doc[lookup]
            tanimoto = self.get_tanimoto(inchi, mol2)
            if tanimoto > coeff_min:
                np.put(coeff_np, min_index, tanimoto)
                top_inchi[min_index] = hash2
                # update min coeff information
                coeff_min = np.amin(coeff_np)
                min_index = np.argmin(coeff_np)
                i += 1

        indices = np.argsort(coeff_np)
        sorted_inchi = []
        for x in (indices[::-1]):
        	sorted_inchi.append(top_inchi[x])
        sorted_coeff = np.sort(coeff_np)[::-1]

        return sorted_coeff, sorted_inchi


    def many_to_many(self, collection_str1='metabolites_meta',
                     collection_str2='metabolites_meta', field1='inchi_deprot',
                     field2='inchi_deprot', lookup1 = 'inchi_hashed', 
                     lookup2 ='inchi_hashed', batch_size = 100, num = 100,
                     no_cursor_timeout = False):
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
                    username = self.username, password = self.password,
                    authSource = self.authSource)
        db_obj = client[self.db]
        final = db_obj[collection_str1]

        projection = {field1: 1, lookup1: 1, '_id': 0}
        _, _, col = self.con_db(collection_str1)
        cursor = col.find({}, projection=projection, batch_size = batch_size, 
                        no_cursor_timeout = no_cursor_timeout)
        count = col.count_documents({})
        total = min(count, self.max_entries)
        
        i = 0
        for doc in cursor:
            if 'similar_compounds' in doc:
                continue
            if i > self.max_entries:
                break
            if self.verbose and i % 100 == 0:
                print('Going through document {} out of {} in collection {}'.format(i, total, collection_str1))
            compound = doc[field1]
            coeff, inchi_hashed = self.one_to_many(compound, lookup = lookup2,
                                             collection_str=collection_str2, field=field2, num=num)
            dic = {}
            for a, b in zip(coeff, inchi_hashed):
                dic[b] = a

            final.update_one( {lookup1: doc[lookup1]},
                            {'$set': {'similar_compounds': dic} },
                            upsert = False)
            i += 1
        cursor.close()

def main():

    db = 'datanator'
    config_file = '/root/host/karr_lab/datanator/.config/config.ini'
    username, password, server, port = server_util.ServerUtil(
        config_file=config_file).get_user_config()
    manager = CalcTanimoto(
        MongoDB=server, replicaSet=None, db=db,
        verbose=True, password=password, username=username,
        result_db = 'datanator')
    manager.many_to_many(no_cursor_timeout = True)


if __name__ == '__main__':
    main()
