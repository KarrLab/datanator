import json  
import pandas as pd 
from datanator_query_python.config import config
from datanator_query_python.util import mongo_util


class QuerySubcellularProteinLocalization(mongo_util.MongoUtil):

    def __init__(self, dataset, directory, max_entries, 
                 query_collection_uniprot, query_collection_taxon_tree,
                 MongoDB, db, username, password):
        super().__init__(MongoDB=MongoDB,
                         db=db,
                         username=username,
                         password=password)
        

        # local dataset to parse 
        self.dataset = dataset
        self.directory = directory
        self.max_entries = max_entries

        # query for UniProt and canon ancestors
        self.query_collection_uniprot = self.db_obj[query_collection_uniprot]
        self.query_collection_taxon_tree = self.db_obj[query_collection_taxon_tree]


    def create_genotype_obj(self):
        """
        queries the taxonomy information for homo sapien for genotype object in json file
        creates the object for genotype
        """
        query = {"tax_id": 9606}
        
        result = self.query_collection_taxon_tree.find_one(query) 
        array_canon_ancestors = []

        # creates an object for each canon ancestor
        for i in range(len(result["canon_anc_ids"])):
            obj_canon_ancestor = {}
            obj_canon_ancestor["ncbi_taxonomy_id"] = result["canon_anc_ids"][i]
            obj_canon_ancestor["name"] = result["canon_anc_names"][i]
            array_canon_ancestors.append(obj_canon_ancestor)

        # initialize a dictionary for genotype 
        obj_genotype = {}

        # taxon object in genotype object
        obj_taxon_in_genotype = {}
        obj_taxon_in_genotype["ncbi_taxonomy_id"] = 9606
        obj_taxon_in_genotype["name"] = "Homo sapiens"
        obj_taxon_in_genotype["canon_ancestors"] = array_canon_ancestors
        obj_genotype["taxon"] = obj_taxon_in_genotype
        
        return obj_genotype # returns one singular object with taxon information

    
    def create_entity_obj(self, gene_name):
        """
        queries gene name to return entity object for each gene
        
        Args:
            gene_name(:obj:'str'): second column of subcellular data, name of gene
        """
        query = {"gene_name": gene_name, "ncbi_taxonomy_id": 9606}
        result = self.query_collection_uniprot.find_one(query)

        # initialize an entity object
        obj_entity = {}
        array_identifier = []
        obj_entity["type"] = "protein"
        obj_entity["name"] = result["protein_name"]
        # initialize identifier object in entity
        obj_in_identifier = {}
        obj_in_identifier["namespace"] = "uniprot_id"
        obj_in_identifier["value"] = result["uniprot_id"]
        array_identifier.append(obj_in_identifier)

        obj_entity["identifiers"] = array_identifier

        return obj_entity

    
    def createJSON(self):
        data = pd.read_csv(self.dataset, delimiter='\t')
        data = data.where(pd.notnull(data), None)


        for i in range(self.max_entries):
            d = {}
            
            genotype = self.create_genotype_obj()
            entity = self.create_entity_obj(gene_name=data.iloc[i, 1])

            d["entity"] = entity
            d["genotype"] = genotype

            d["values"] = []
            for column_name in data.columns[1:]:
                if column_name != "SeqID" and column_name != "PSortVersion":
                    values_dict = {}
                    values_dict["type"] = column_name
                    values_dict["value"] = data[column_name].iloc[i]
                    d["values"].append(values_dict)

            d["source"] = [
                {
                    "namespace": "Protein Atlas",
                    "value": "subcellular_location"
                }
            ]
            d["schema_version"] = "2.0"

            with open(self.directory+"/{}.json".format(data.iloc[i,1]), "w+") as f:
                json.dump(d, f)




def main():
    conf = config.Justin()
    username = conf.USERNAME 
    password = conf.PASSWORD 
    MongoDB = conf.SERVER 

    # location of dataset and where files will be created
    dataset = "datanator/docs/protein_localization/subcellular_localization/subcellular_location.tsv"
    directory = "./datanator/docs/protein_localization/subcellular_localization/subcellularJSONfiles"
    max_entries = 100

    # database and collection being queried
    db = "datanator"
    query_collection_taxon_tree = "taxon_tree"
    query_collection_uniprot = "uniprot"


    src = QuerySubcellularProteinLocalization(dataset, directory, max_entries, 
                                              query_collection_uniprot, query_collection_taxon_tree,
                                              MongoDB, db, username, password)
    
    src.createJSON()

if __name__ == "__main__":
    main()