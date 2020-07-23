from datanator_query_python.util import mongo_util
from datanator_query_python.config import config


class TransformMetabolitesMeta(mongo_util.MongoUtil):
    def __init__(self, MongoDB=None,
                       db=None,
                       username=None,
                       password=None,
                       max_entries=float('inf')):
        super().__init__(MongoDB=MongoDB,
                         db=db,
                         username=username,
                         password=password)
        self.max_entries = max_entries
        self.entity_col = self.db_obj["entity"]
        self.identifier_col = self.db_obj["identifier"]
        self.obs_col = self.db_obj["observation"]

    def process_docs(self, skip=0):
        query = {}
        projection = {"_id": 0}
        docs = self.client["datanator-test"]["metabolites_meta"].find(filter=query, projection=projection)
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            query = {"identifiers": {"$elemMatch": {'namespace': 'inchikey', 'value': doc["InChI_Key"]}}}
            entity = self.build_entity(doc)
            ob_y, ob_e = self.build_obs(doc)
            self.entity_col.update_one(query,
                                       {"$set": {"type": entity["type"],
                                                 "name": entity["name"],
                                                 "description": entity["description"],
                                                 "schema_version": "2.0"},
                                        "$addToSet": {"identifiers": {"$each": entity["identifiers"]},
                                                      "synonyms": {"$each": entity["synonyms"]},
                                                      "related": {"$each": entity["related"]},
                                                      "similarity": {"$each": entity["similarity"]},
                                                      "structures": {"$each": entity["structures"]}}},
                                        upsert=True)

            if ob_y != {}:
                query = {"$and": [{"identifier": {'namespace': 'inchikey', 'value': doc["InChI_Key"]}},
                                {"source": {"$elemMatch": ob_y["source"][0]}}]}
                self.obs_col.update_one(query,
                                        {"$set": {"genotype": ob_y["genotype"],
                                                  "entity": ob_y["entity"],
                                                  "schema_version": "2.0",
                                                  "identifier": ob_y["entity"]["identifiers"][0]},
                                         "$addToSet": {"values": {"$each": ob_y["values"]},
                                                       "source": {"$each": ob_y["source"]}}},
                                         upsert=True)
            if ob_e != {}:
                query = {"$and": [{"identifier": {'namespace': 'inchikey', 'value': doc["InChI_Key"]}},
                                {"source": {"$elemMatch": ob_e["source"][0]}}]}
                self.obs_col.update_one(query,
                                        {"$set": {"genotype": ob_e["genotype"],
                                                  "entity": ob_e["entity"],
                                                  "schema_version": "2.0",
                                                  "identifier": ob_e["entity"]["identifiers"][0]},
                                         "$addToSet": {"values": {"$each": ob_e["values"]},
                                                       "source": {"$each": ob_e["source"]}}},
                                         upsert=True)

    def build_entity(self, obj):
        """Build entity object from obj.
        Go into entity collection

        Args:
            obj (:obj:`Obj`): source object.

        Return:
            (:obj:`Obj`), e.g.
            {
                "entity": {
                    "type": "metabolite",
                    "name": "2-Ketobutyric acid",
                    "synonyms": [],
                    "identifiers": [{}... {}],
                    "related": [{"namespace": "", "value": ""}]
                }
            }        
        """
        entity = {}
        entity["structures"] = []
        entity["type"] = "metabolite"
        entity["name"] = obj.get("name")
        entity["synonyms"] = obj.get("synonyms", [])
        entity["identifiers"] = []
        entity["identifiers"].append({"namespace": "inchikey",
                                    "value": obj.get("InChI_Key")})
        entity["identifiers"].append({"namespace": "biocyc_id",
                                    "value": obj.get("biocyc_id")})
        entity["identifiers"].append({"namespace": "cas_registry_number",
                                    "value": obj.get("cas_registry_number")})
        entity["identifiers"].append({"namespace": "chemical_formula",
                                    "value": obj.get("chemical_formula")})
        entity["identifiers"].append({"namespace": "chebi_id",
                                    "value": obj.get("chebi_id")})
        entity["identifiers"].append({"namespace": "chemspider_id",
                                    "value": obj.get("chemspider_id")})
        entity["identifiers"].append({"namespace": "hmdb_id",
                                    "value": obj.get("hmdb_id")})
        entity["structures"].append({"format": "inchi",
                                    "value": obj.get("inchi")})
        entity["identifiers"].append({"namespace": "kegg_id",
                                    "value": obj.get("kegg_id")})
        entity["identifiers"].append({"namespace": "m2m_id",
                                    "value": obj.get("m2m_id")})
        entity["identifiers"].append({"namespace": "pubchem_compound_id",
                                    "value": obj.get("pubchem_compound_id")})
        entity["structures"].append({"format": "smiles",
                                    "value": obj.get("smiles")})
        entity["identifiers"].append({"namespace": "ymdb_id",
                                    "value": obj.get("ymdb_id")})
        entity["identifiers"].append({"namespace": "average_molecular_weight",
                                    "value": obj.get("average_molecular_weight")})
        entity["similarity"] = []
        for sim in obj['similar_compounds']:
            entity["similarity"].append({"namespace": "inchikey",
                                        "value": sim["inchikey"],
                                        "similarity_score": sim["similarity_score"]})
        entity["description"] = obj.get("description")
        entity["related"] = []
        if obj["pathways"] is None:
            return entity
        pathways = obj["pathways"]["pathway"]
        if isinstance(pathways, list):
            for path in pathways:
                related = {"description": path.get("name"),
                            "namespace": "kegg_map_id",
                            "value": path.get("kegg_map_id")}
                entity["related"].append(related)
                self.identifier_col.update_one({"$and": [{"namespace": "kegg_map_id"},
                                                        {"value": path.get("kegg_map_id")}]},
                                                {"$set": related},
                                                upsert=True)
        else:
            path = pathways
            related = {"description": path.get("name"),
                        "namespace": "kegg_map_id",
                        "value": path.get("kegg_map_id")}
            entity["related"].append(related)
            self.identifier_col.update_one({"$and": [{"namespace": "kegg_map_id"},
                                                    {"value": path.get("kegg_map_id")}]},
                                            {"$set": related},
                                            upsert=True)

        return entity

    def build_obs(self, obj):
        """Build observation objects from obj.
        Go into observations collection.

        Args:
            obj (:obj:`Obj`): source object.

        Return:
            obj(:obj:`Obj`)
            {
                "entity": "weight",
                "value": [],
                "source": {}, ...
            }
        """
        cell_locs = obj.get("cellular_locations", [])
        entity = {"type": "metabolite",
                  "name": obj.get("name"),
                  "identifiers": [{"namespace": "inchikey",
                                   "value": obj.get("InChI_Key")}],
                  "schema_version": "2.0"}
        values_e = []
        values_y = []
        if cell_locs != [] and cell_locs is not None: 
            for loc in cell_locs:
                if isinstance(loc["cellular_location"], str):
                    if "YMDB" in loc["reference"]:
                        values_y.append({"type": "cellular_location",
                                        "value": loc["cellular_location"]})
                    else:
                        values_e.append({"type": "cellular_location",
                                        "value": loc["cellular_location"]})
                else:
                    for l in loc["cellular_location"]:
                        if "YMDB" in loc["reference"]:
                            values_y.append({"type": "cellular_location",
                                            "value": l})
                        else:
                            values_e.append({"type": "cellular_location",
                                            "value": l})
        for prop in obj["property"]:
            if obj.get("ymdb_id") is not None:
                values_y.append({"type": prop["kind"],
                                "value": prop["value"]})
            else:
                values_e.append({"type": prop["kind"],
                                "value": prop["value"]})    
        ob_y = {}
        ob_e = {}
        if values_y != []:
            genotype = {"taxon": {"ncbi_taxonomy_id": 4932,
                                "name": "Saccharomyces cerevisiae",
                                "canon_ancestors": [{"ncbi_taxonomy_id": 131567,
                                                    "name": "cellular organisms"},
                                                    {"ncbi_taxonomy_id": 2759,
                                                    "name": "Eukaryota"},
                                                    {"ncbi_taxonomy_id": 4751,
                                                    "name": "Fungi"},
                                                    {"ncbi_taxonomy_id": 4890,
                                                    "name": "Ascomycota"},
                                                    {"ncbi_taxonomy_id": 4891,
                                                    "name": "Saccharomycetes"},
                                                    {"ncbi_taxonomy_id": 4892,
                                                    "name": "Saccharomycetales"},
                                                    {"ncbi_taxonomy_id": 4893,
                                                    "name": "Saccharomycetaceae"},
                                                    {"ncbi_taxonomy_id": 4930,
                                                    "name": "Saccharomyces"}]}}
            ob_y = {"entity": entity,
                    "genotype": genotype,
                    "values": values_y,
                    "source": [{"namespace": "YMDB", "value": obj.get("ymdb_id"), "level": "secondary"}],
                    "schema_version": "2.0"}
        if values_e != []:
            genotype = {"taxon": {"ncbi_taxonomy_id": 562,
                                "name": "Escherichia coli",
                                "canon_ancestors": [{"ncbi_taxonomy_id": 131567,
                                                    "name": "cellular organisms"},
                                                    {"ncbi_taxonomy_id": 2,
                                                    "name": "Bacteria"},
                                                    {"ncbi_taxonomy_id": 1224,
                                                    "name": "Proteobacteria"},
                                                    {"ncbi_taxonomy_id": 1236,
                                                    "name": "Gammaproteobacteria"},
                                                    {"ncbi_taxonomy_id": 91347,
                                                    "name": "Enterobacterales"},
                                                    {"ncbi_taxonomy_id": 543,
                                                    "name": "Enterobacteriaceae"},
                                                    {"ncbi_taxonomy_id": 561,
                                                    "name": "Escherichia"}]}}
            ob_e = {"entity": entity,
                    "genotype": genotype,
                    "values": values_e,
                    "source": [{"namespace": "ECMDB", "value": obj.get("m2m_id"), "level": "secondary"}],
                    "schema_version": "2.0"}
        return ob_y, ob_e
        