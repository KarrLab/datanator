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

    def get_tanimoto(self, mol1, mol2):
        _, _, col_ecmdb = self.con_db('ecmdb')
        _, _, col_ymdb = self.con_db('ymdb')
