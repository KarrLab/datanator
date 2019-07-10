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
        self.db = query_nosql.DataQuery(MongoDB=MongoDB, replicaSet=replSet, db=db,
                                                       username=username, password=password)
        self.metab_db = query_nosql.QueryMetabolitesMeta(MongoDB=MongoDB, replicaSet=replSet, db=db,
                                                                        username=username, password=password)
        self.sabio_manager = query_nosql.QuerySabio(MongoDB=MongoDB, replicaSet=replSet, db=db,
                                                    username=username, password=password)
        self.tax_db = query_nosql.QueryTaxonTree(MongoDB=MongoDB, replicaSet=replSet, db=db,
                                                            username=username, password=password)
        self.mongoutil_manager = mongo_util.MongoUtil(MongoDB=MongoDB, replicaSet=replSet, db=db,
                                                      username=username, password=password)
        self.chem_manager = chem_util.ChemUtil()
        self.file_manager = file_util.FileUtil()

    def get_ecmdb_entries(self, m2m_ids, organism):
        query = { 'm2m_id': {"$in": m2m_ids} } 
        projection = {'_id': 0}
        cursor = self.db.doc_feeder(collection_str='ecmdb', query=query ,projection=projection)
        list_e_coli = []
        for doc in cursor:
            if doc['concentrations']:
                if len(doc['concentrations']['concentration']) > 0:
                    list_e_coli.append(doc)
        return(list_e_coli)

    def get_ymdb_entries(self, ymdb_ids, organism):
        query = { 'ymdb_id': {"$in": ymdb_ids} } 
        projection = {'_id': 0}
        cursor2 = self.db.doc_feeder(collection_str='ymdb', query=query ,projection=projection)
        list_yeast = []
        for doc in cursor2:
            if doc['concentrations']:
                if len(doc['concentrations']['concentration']) > 0:
                    list_yeast.append(doc)
        return(list_yeast)

    def get_conc_ids(self, list_names):
        m2m_ids = []
        ymdb_ids = []
        for entry in self.metab_db.get_metabolite_inchi(list_names):
            if entry['m2m_id']:
                m2m_ids.append(entry['m2m_id'])
            if entry['ymdb_id']:
                ymdb_ids.append(entry['ymdb_id'])
        return(m2m_ids, ymdb_ids)
    """
    def get_generic_names(self, molecule_name):
        raw, result = self.metab_db.get_metabolite_similar_compounds([molecule_name], num = 5, threshold = 0.6)
        list_names = []
        list_scores = []
        for key in result[0]:
            list_names.append(key)
            list_scores.append(result[0][key])
        return(list_names, list_scores)
    """

    def get_generic_concs(self, molecule_name, organism):
        raw, result = self.metab_db.get_metabolite_similar_compounds([molecule_name], num = 30, threshold = 0.6)
        print(raw)
        print(result)
        
        list_names = []
        list_scores = []
        for key in result[0]:
            if key != 'None':
                list_names.append(key)
                list_scores.append(result[0][key])
        print(list_names)

        m2m_ids = []
        ymdb_ids = []
        id_to_score = {}
        list_db_id = self.metab_db.get_metabolite_inchi(list_names)
        for i in range(len(list_db_id)):
            entry = list_db_id[i]
            if entry['m2m_id']:
                m2m_ids.append(entry['m2m_id'])
                id_to_score[entry['m2m_id']] = list_scores[i]
            if entry['ymdb_id']:
                ymdb_ids.append(entry['ymdb_id'])
                id_to_score[entry['ymdb_id']] = list_scores[i]

        #m2m_ids, ymdb_ids = self.get_conc_ids(list_names)
        print(m2m_ids)
        print(ymdb_ids)
        ecmdb_data = self.get_ecmdb_entries(m2m_ids, organism)
        ymdb_data = self.get_ymdb_entries(ymdb_ids, organism)
        
        
        for entry in ecmdb_data:
            entry["tanitomo_similarity"] = id_to_score[entry["m2m_id"]]
        for entry in ymdb_data:
            entry["tanitomo_similarity"] = id_to_score[entry["ymdb_id"]]

        return(ecmdb_data, ymdb_data)

    def molecule_name_query(self, molecule_name, organism, abstract_default=False):

        #list_metabolites = ["ATP", "GTP"]
        list_metabolites = [molecule_name]
        print(self.metab_db.get_metabolite_inchi(list_metabolites))
        m2m_ids, ymdb_ids = self.get_conc_ids(list_metabolites)
        ecmdb_data = self.get_ecmdb_entries(m2m_ids, organism)
        ymdb_data = self.get_ymdb_entries(ymdb_ids, organism)


        for entry in ecmdb_data:
            entry["tanitomo_similarity"] = 1
        for entry in ymdb_data:
            entry["tanitomo_similarity"] = 1


        response = []
        if (not (ecmdb_data or ymdb_data)) or abstract_default:
            new_ecmdb_data, new_ymdb_data  = self.get_generic_concs(molecule_name, organism)
            ecmdb_data = ecmdb_data + new_ecmdb_data
            ymdb_data = ymdb_data + new_ymdb_data
        else:
            pass
        response.append(ecmdb_data)
        response.append(ymdb_data)


        anc, dist = self.tax_db.get_common_ancestor(organism, "Escherichia coli")
        for doc in ecmdb_data:
            doc["taxon_distance"] = dist[0]
        anc, dist = self.tax_db.get_common_ancestor(organism, "Saccharomyces cerevisiae")
        for doc in ymdb_data:
            doc["taxon_distance"] = dist[0]



        result_id, result_name = self.tax_db.get_anc_by_name([organism])
        result_name[0].append(organism)
        response.append(result_name)
        #response = self.inchi_query(inchi, organism, string)
        return(response)

