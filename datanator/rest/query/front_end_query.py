import datanator.config.core
from datanator.core import query_nosql
from datanator.util import chem_util
from datanator.util import file_util
from datanator.util import mongo_util
import sys
# print("hello")

# sys.path.insert(0, '~/karr_lab/datanator_frontend2/datanator')

class QueryFrontEnd:
    def __init__(self):
        db = 'datanator'
        username = datanator.config.core.get_config()[
            'datanator']['mongodb']['user']
        password = datanator.config.core.get_config(
        )['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config(
        )['datanator']['mongodb']['server']
        port = datanator.config.core.get_config(
        )['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config(
        )['datanator']['mongodb']['replSet']
        self.dataquery_manager = query_nosql.DataQuery(MongoDB=MongoDB, replicaSet=replSet, db=db,
                                                       username=username, password=password)
        self.metabolitesmeta_manager = query_nosql.QueryMetabolitesMeta(MongoDB=MongoDB, replicaSet=replSet, db=db,
                                                                        username=username, password=password)
        self.sabio_manager = query_nosql.QuerySabio(MongoDB=MongoDB, replicaSet=replSet, db=db,
                                                    username=username, password=password)
        self.taxontree_manager = query_nosql.QueryTaxonTree(MongoDB=MongoDB, replicaSet=replSet, db=db,
                                                            username=username, password=password)
        self.mongoutil_manager = mongo_util.MongoUtil(MongoDB=MongoDB, replicaSet=replSet, db=db,
                                                      username=username, password=password)
        self.chem_manager = chem_util.ChemUtil()
        self.file_manager = file_util.FileUtil()

    def string_query(self, string):
        query = string
        results = self.dataquery_manager.find_text(query, collection='ecmdb')
        return(results)

    def inchi_query_metabolite(self, string):
        # query for inchi string
        inchi_deprot = self.chem_manager.simplify_inchi(inchi=string)
        inchi_hashed = self.chem_manager.hash_inchi(inchi=inchi_deprot)
        ids = self.metabolitesmeta_manager.get_ids_from_hash(inchi_hashed)

        list_jsons = []
        if ids['m2m_id'] != None:
            query = {'m2m_id': ids['m2m_id']}
            projection = {'_id': 0}
            _, _, col = self.mongoutil_manager.con_db('ecmdb')
            doc = col.find_one(filter=query, projection=projection)
            list_jsons.append(doc)

        if ids['ymdb_id'] != None:
            query = {'ymdb_id': ids['ymdb_id']}
            projection = {'_id': 0}
            _, _, col = self.mongoutil_manager.con_db('ymdb')
            doc = col.find_one(filter=query, projection=projection)
            list_jsons.append(doc)
        return(list_jsons)

    def inchi_query_organism(self, string, organism):
        ''' Find metabolite (defined by string) concentration
            in e.coli and yeast, return documents
            Args:
                string: inchi value
                organism: concentration value missing in this organism
        '''

        tree = self.taxontree_manager.get_anc_id_by_name([organism])
        for entry in tree:
            print(entry)
        anc, dist = self.taxontree_manager.get_common_ancestor(
            organism, "Escherichia coli")
        distance_e = dist[0]
        anc, dist = self.taxontree_manager.get_common_ancestor(
            organism, "Saccharomyces cerevisiae")
        distance_y = dist[0]

        inchi_deprot = self.chem_manager.simplify_inchi(inchi=string)
        inchi_hashed = self.chem_manager.hash_inchi(inchi=inchi_deprot)
        ids = self.metabolitesmeta_manager.get_ids_from_hash(inchi_hashed)
        list_jsons = []

        if ids['m2m_id'] != None:
            query = {'m2m_id': ids['m2m_id']}
            projection = {'_id': 0}
            _, _, col = self.mongoutil_manager.con_db('ecmdb')
            doc = col.find_one(filter=query, projection=projection)
            doc["taxon_distance"] = distance_e
            list_jsons.append(doc)

        if ids['ymdb_id'] != None:
            query = {'ymdb_id': ids['ymdb_id']}
            projection = {'_id': 0}
            _, _, col = self.mongoutil_manager.con_db('ymdb')
            doc = col.find_one(filter=query, projection=projection)
            doc["taxon_distance"] = distance_y
            list_jsons.append(doc)
        return(list_jsons)

    def molecule_name_query(self, string, organism):

        response = None
        # try:
        inchi = self.metabolitesmeta_manager.get_metabolite_inchi([string])[
            0]
        response = self.inchi_query_organism(inchi, organism, string)
        # except Exception as e:
        # print(e)

        return(response)
