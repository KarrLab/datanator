from datanator_query_python.util import mongo_util
from datanator_query_python.config import config
import copy
import numpy as np
from numpy.lib.arraysetops import isin


class Transform(mongo_util.MongoUtil):
    def __init__(self,
                 MongoDB=None,
                 db=None,
                 des_col=None,
                 username=None,
                 password=None,
                 max_entries=float('inf'),
                 verbose=True):
        super().__init__(MongoDB=MongoDB,
                         db=db,
                         username=username,
                         password=password)
        self.max_entries = max_entries
        self.col = des_col
        self.db = db
        self.verbose = verbose

    def process_docs(self,
                     col, 
                     skip=0):
        """Processing documents and transform.

        Args:
            col(:obj:`str`): Name of the source collection.
        """
        query = {}
        projection = {"_id": 0}
        docs = self.client["datanator-test"][col].find(filter=query, projection=projection)
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            if i % 20 == 0 and self.verbose:
                print("Processing doc {}".format(i))
            if col == "uniprot":
                entity = self.build_uniprot_entity(doc)
                obs = self.build_uniprot_observation(doc)
                self.update_entity(entity,
                                   entity["identifiers"][0],
                                   db=self.db)
                if obs != {}:
                    self.update_observation(obs,
                                            obs["source"][0])


    def build_uniprot_entity(self, obj):
        """Build entity from uniprot collection.

        Args:
            (:obj:`Obj`): object from which entity object will be built.

        Return:
            (:obj:`Obj`): entity object.
        """
        _type = "protein"
        identifiers = []
        related = []
        genotype = {}
        name = obj.get("protein_name")
        identifiers.append({"namespace": "uniprot_id",
                            "value": obj.get("uniprot_id")})
        for o in obj.get("add_id"):
            identifiers.append({"namespace": o.get("name_space"),
                                "value": o.get("value")})
        canon_ancestors = []
        for _id, name in zip(obj.get("canon_anc_ids"), obj.get("canon_anc_names")):
            canon_ancestors.append({"ncbi_taxonomy_id": _id,
                                    "name": name})
        genotype["taxon"] = {"ncbi_taxonomy_id": obj.get("ncbi_taxonomy_id"),
                            "obj": obj.get("species_name"),
                            "canon_ancestors": canon_ancestors}
        structures = [{"format": "canonical_sequence",
                       "value": obj.get("canonical_sequence")}]
        mod = obj.get("modifications")            
        if mod is not None:
            if isinstance(mod, dict):
                if mod.get("concrete") is True and mod.get("monomeric_form_issues") == np.NaN and mod.get("pro_issues") == np.NaN:
                    identifiers.append({"namespace": "pro_id",
                                        "value": mod.get("pro_id")})
                    structures.append({"format": "processed_sequence_iubmb",
                                       "value": mod.get("processed_sequence_iubmb"),
                                       "molecular_weight": mod.get("processed_molecular_weight"),
                                       "charge": mod.get("processed_charge"),
                                       "formula": mod.get("processed_formula"),
                                       "source": [{"namespace": "pro_id",
                                                   "value": mod.get("pro_id"),
                                                   "level": "secondary"},
                                                  {"namespace": "doi",
                                                   "value": mod.get("reference")["doi"],
                                                   "level": "primary"}]})
                    structures.append({"format": "modified_sequence_abbreviated_bpforms",
                                       "value": mod.get("modified_sequence_abbreviated_bpforms"),
                                       "molecular_weight": mod.get("modified_molecular_weight"),
                                       "charge": mod.get("modified_charge"),
                                       "formula": mod.get("modified_formula"),
                                       "modification": {
                                           "description": mod.get("modifications"),
                                           "formula": mod.get("modifications_formula"),
                                           "weight": mod.get("modifications_molecular_weight"),
                                           "charge": mod.get("modifications_charge")
                                        },
                                       "source": [{"namespace": "pro_id",
                                                   "value": mod.get("pro_id"),
                                                   "level": "secondary"},
                                                  {"namespace": "doi",
                                                   "value": mod.get("reference")["doi"],
                                                   "level": "primary"}]})
                    structures.append({"format": "modified_sequence_bpforms",
                                       "value": mod.get("modified_sequence_bpforms")})
            elif isinstance(mod, list):
                for o in mod:
                    if o.get("concrete") and np.isnan(o.get("monomeric_form_issues")) and np.isnan(o.get("pro_issues")):
                        identifiers.append({"namespace": "pro_id",
                                            "value": o.get("pro_id")})
                        structures.append({"format": "processed_sequence_iubmb",
                                        "value": o.get("processed_sequence_iubmb"),
                                        "molecular_weight": o.get("processed_molecular_weight"),
                                        "charge": o.get("processed_charge"),
                                        "formula": o.get("processed_formula"),
                                        "source": [{"namespace": "pro_id",
                                                    "value": o.get("pro_id"),
                                                    "level": "secondary"},
                                                    {"namespace": "doi",
                                                    "value": o.get("reference")["doi"],
                                                    "level": "primary"}]})
                        structures.append({"format": "modified_sequence_abbreviated_bpforms",
                                        "value": o.get("modified_sequence_abbreviated_bpforms"),
                                        "molecular_weight": o.get("modified_molecular_weight"),
                                        "charge": o.get("modified_charge"),
                                        "formula": o.get("modified_formula"),
                                        "modification": {
                                            "description": o.get("modifications"),
                                            "formula": o.get("modifications_formula"),
                                            "weight": o.get("modifications_molecular_weight"),
                                            "charge": o.get("modifications_charge")
                                            },
                                        "source": [{"namespace": "pro_id",
                                                    "value": o.get("pro_id"),
                                                    "level": "secondary"},
                                                    {"namespace": "doi",
                                                    "value": o.get("reference")["doi"],
                                                    "level": "primary"}]})
                        structures.append({"format": "modified_sequence_bpforms",
                                        "value": o.get("modified_sequence_bpforms")})
        related.append({"namespace": "ec",
                        "value": obj.get("ec_number")})
        identifiers.append({"namespace": "entrez_id",
                            "value": obj.get("entrez_id")})
        identifiers.append({"namespace": "entry_name",
                            "value": obj.get("entry_name")})
        related.append({"namespace": "gene_name",
                        "value": obj.get("gene_name")})
        for n in obj.get("ko_name"):
            related.append({"namespace": "ko_name",
                            "value": n})  
        related.append({"namespace": "ko_number",
                        "value": obj.get("ko_number")})
        return {"type": _type,
                "name": name,
                "synonyms": [],
                "identifiers": identifiers,
                "related": related,
                "genotype": genotype,
                "structures": structures,
                "schema_version": "2.0"}

    def build_uniprot_observation(self, obj):
        """Build observation from uniprot collection.

        Args:
            (:obj:`Obj`): object from which observation object will be built.

        Return:
            (:obj:`Obj`): observation object.
        """
        abundances = obj.get("abundances", [])
        schema_version = "2.0"
        result = {}
        if len(abundances) == 0:
            return result
        else:
            values = []
            for a in abundances:
                values.append({"type": "protein_abundance",
                               "value": a.get("abundance"),
                               "units": "ppm",
                               "organ": a.get("organ")})
            entity = {"type": "protein",
                      "name": obj.get("protein_name"),
                      "identifiers": [{"namespace": "uniprot_id",
                                      "value": obj.get("uniprot_id")}],
                      "schema_version": schema_version}
            return {"entity": entity,
                    "values": values,
                    "source": [{"namespace": "paxdb",
                                "value": obj.get("uniprot_id")}],
                    "identifier": {"namespace": "uniprot_id",
                                   "value": obj.get("uniprot_id")},
                    "schema_version": schema_version}


def main():
    pass

if __name__ == "__main__":
    main()