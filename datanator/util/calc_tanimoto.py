import pybel
from datanator.util import mongo_util
from datanator.util import server_util


class CalcTanimoto(mongo_util.MongoUtil):
    '''Calculating the tanitomo similarity matrix 
            given two compound collections e.g.
            ECMDB YMDB
    '''

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf'), username=None,
                 password=None, authSource='admin'):
        super(CalcTanimoto, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet,
                                           db=db, verbose=verbose, max_entries=max_entries, username=username,
                                           password=password, authSource=authSource)

    def get_tanimoto(self, mol1, mol2, str_format='inchi'):
        '''Calculates tanimoto coefficients between
                two molecules, mol1 and mol2
                Args:
                        mol1: molecule 1 in some format
                        mol2: molecule 2 in same format as molecule 1
                Return:
                        tani: Tanimoto coefficient
        '''
        inchi = [mol1, mol2]
        mols = [pybel.readstring(str_format, x) for x in inchi]
        fps = [x.calcfp() for x in mols]
        return round((fps[0] | fps[1]), 3)
