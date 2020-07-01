from datanator_query_python.util import mongo_util
from datanator_query_python.config import config
import copy


class TMC(mongo_util.MongoUtil):
    def __init__(self,
                 MongoDB=None,
                 db=None,
                 des_col=None,
                 username=None,
                 password=None,
                 max_entries=float('inf')):
        super().__init__(MongoDB=MongoDB,
                         db=db,
                         username=username,
                         password=password)
        self.max_entries = max_entries
        self.des = self.db_obj[des_col]

    def process_docs(self, skip=0):
        query = {}
        projection = {"_id": 0}
        docs = self.client["datanator-test"]["metabolite_concentrations"].find(filter=query, projection=projection)
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            sets = self.build_conc_observation(doc)
            for _set in sets:
                if _set[0] != {}:
                    print(_set[0])
                    query = {"$and": [{"identifier": {'namespace': 'inchikey', 'value': doc["inchikey"]}},
                                      {"source": {"$elemMatch": _set[0]["source"][0]}}]}
                    self.des.update_one(query,
                                        {"$set": {"genotype": _set[0]["genotype"],
                                                  "entity": _set[0]["entity"],
                                                  "schema_version": "2.0",
                                                  "identifier": _set[0]["entity"]["identifiers"][0],
                                                  "environment": _set[0]["environment"]},
                                        "$addToSet": {"values": {"$each": _set[0]["values"]},
                                                      "source": {"$each": _set[0]["source"]}}},
                                        upsert=True)
                if _set[1] != {}:
                    query = {"$and": [{"identifier": {'namespace': 'inchikey', 'value': doc["inchikey"]}},
                                      {"source": {"$elemMatch": _set[1]["source"][0]}}]}
                    self.des.update_one(query,
                                        {"$set": {"genotype": _set[1]["genotype"],
                                                  "entity": _set[1]["entity"],
                                                  "schema_version": "2.0",
                                                  "identifier": _set[1]["entity"]["identifiers"][0]},
                                        "$addToSet": {"values": {"$each": _set[1]["values"]},
                                                      "source": {"$each": _set[1]["source"]}}},
                                        upsert=True)

    def build_conc_observation(self, obj):
        """Build concentration observation object(s) from documents in
        metabolite_concentrations collection.

        Args:
            obj (:obj:`Obj`): Document.

        Return:
            (:obj:`Iter`)
        """
        identifier = {"namespace": "inchikey",
                      "value": obj["inchikey"]}
        entity = {"name": obj["metabolite"],
                  "type": "metabolite",
                  "identifiers": [identifier]}
        schema_version = "2.0"        
        for conc in obj["concentrations"]:
            genotype = {}
            environment = {}
            ob_val = {}
            values = []
            aff_values = []
            if conc.get("affinities") is not None:                
                for aff in conc.get("affinities"):
                    entity["structure"] = {"format": "EC",
                                           "value": aff.get("ec_number")}
                    aff_values.append({"type": "k_i",
                                        "substrate": entity,
                                        "value": aff.get("k_i"),
                                        "units": aff.get("unit")})
                    aff_values.append({"type": "k_m",
                                        "substrate": entity,
                                        "value": aff.get("k_m"),
                                        "units": aff.get("unit")})
                error = 0
            else: # has env information
                genotype["growthPhase"] = conc.get("growth_status")
                environment["media"] = conc.get("growth_media")
                environment["growthSystem"] = conc.get("growth_system")
                error = conc.get("error")

            unit = conc.get('unit')
            value = conc.get("concentration")
            
            #unit conversion
            if unit == "M":
                unit = unit
                value = float(value)
                error = float(error)                
            else:
                unit = "M"
                value = float(value) / 1000
                error = float(error) / 1000
            ob_val["units"] = unit
            ob_val["value"] = value
            ob_val["uncertainty"] = error
            ob_val["type"] = "concentration"
            values.append(ob_val)
            canon_ancestors = []
            for _id, name in zip(conc.get("canon_anc_ids"), conc.get("canon_anc_names")):
                canon_ancestors.append({"ncbi_taxonomy_id": _id,
                                        "name": name})
            genotype["taxon"] = {"ncbi_taxonomy_id": conc.get("ncbi_taxonomy_id"),
                                "name": conc.get("species_name"),
                                "canon_ancestors": canon_ancestors}
            source = conc.get("reference")
            if source.get("text"):
                source["text"] = source["text"].replace("'", "”")
            source["value"] = source["id"]           
            source.pop("id")
            source = [source]
            if aff_values != []: 
                aff_obs = {"entity": entity,
                           "values": aff_values,
                           "genotype": genotype,
                           "source": source,
                           "schema_version": schema_version,
                           "identifier": identifier}
            else:
                aff_obs = {}
            new_ent = copy.deepcopy(entity)
            new_ent.pop("structure", None)            
            conc_obs = {"entity": new_ent,
                        "values": values,
                        "genotype": genotype,
                        "environment": environment,
                        "source": source,
                        "schema_version": schema_version,
                        "identifier": identifier}            

            yield (conc_obs, aff_obs)