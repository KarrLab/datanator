import tempfile
'''
Cron jobs for updating mongoDB datanator database
'''
import pymongo
from datanator.util import (mongo_util, index_collection)
from datanator.data_source import (corum_nosql, intact_nosql, metabolite_nosql,
                                   pax_nosql, sabio_rk_nosql, sabio_rk, uniprot_nosql,
                                   sqlite_to_json)
import os


if __name__ == '__main__':
    cache_dirname = tempfile.mkdtemp()
    MongoDB = 'mongodb://mongo:27017'
    db = 'aggregate_test'
    replicaSet = 'rs0'
    verbose = True
    index_manager = index_collection.IndexCollection(
        cache_dirname, MongoDB, replicaSet, db, verbose)

    '''Fill corum collection and index
    '''
    fill_corum = corum_nosql.CorumNoSQL(
        cache_dirname, MongoDB, db, replicaSet='rs0', verbose=verbose)
    fill_corum.load_content()
    index_manager.index_corum('corum')

    '''Fill ECMDB collection
    '''
    fill_ecmdb = metabolite_nosql.MetaboliteNoSQL(
        cache_dirname, 'ecmdb', MongoDB, db, verbose=verbose)
    fill_ecmdb.write_to_json()
    index_manager.index_strdb('ecmdb')

    '''Fill YMDB collection
    '''
    fill_ymdb = metabolite_nosql.MetaboliteNoSQL(
        cache_dirname, 'ymdb', MongoDB, db, verbose=verbose)
    fill_ymdb.write_to_json()
    index_manager.index_strdb('ymdb')

    '''Fill INTACT complex and interaction collections
    '''
    fill_intact = intact_nosql.IntActNoSQL(
        cache_dirname, MongoDB, db, replicaSet='rs0', verbose=verbose)
    fill_intact.load_content()
    index_manager.index_strdb('intact_interaction')
    index_manager.index_intact_complex('intact_complex')

    '''Fill PAX collection
    '''
    fill_pax = pax_nosql.PaxNoSQL(cache_dirname, MongoDB, db, verbose=verbose)
    fill_pax.load_content()
    index_manager.index_pax('pax')

    '''Fill Uniprot collection
    '''
    fill_uniprot = uniprot_nosql.UniprotNoSQL(MongoDB, db)
    fill_uniprot.load_uniprot()
    index_manager.index_uniprot('uniprot')

    # '''Fill sabio_rk collection
    # 	Two Steps:
    # 			Build SQL database
    # 			Parse local .sqlite files into corresponding json files
    # 			Parse SQL into MongoDB
    # '''
    # fill_sabio_sql = sabio_rk.SabioRk(cache_dirname=cache_dirname).load_content()
    # print('passed sabio sql')
    # write_sabio_json(cache_dirname)
    # fill_sabio_nosql = sabio_rk_nosql.SabioRkNoSQL(
    #     db, MongoDB, cache_dirname, None, verbose=verbose)
    # file_names, file_dict = fill_sabio_nosql.load_json()
    # fill_sabio_nosql.make_doc(file_names, file_dict)
    # index_manager.index_sabio('sabio_rk')


    shutil.rmtree(cls.cache_dirname)


def write_sabio_json(self, cache_dirname):
    database = os.path.join(cache_dirname, 'SabioRk.sqlite')
    query = "select * from "
    collection_dir = cache_dirname + '/SabioRk/'
    os.makedirs(os.path.dirname(cache_dirname), exist_ok=True)

    temp = sqlite_to_json.SQLToJSON(database, query)
    tables = temp.table()

    i = 0
    for table in tables:
        file_name = os.path.join(cache_dirname + table + '.json')
        result = sqlite_to_json.SQLToJSON(
            database, query).query_table(table, False)
        with open(file_name, "w") as f:
            f.write(json.dumps(result, indent=4))
        print('Written {} json'.format(i))
        i += 1
