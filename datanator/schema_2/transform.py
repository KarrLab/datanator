from datanator_query_python.util import mongo_util
from datanator_query_python.config import config
import copy
import numpy as np
import re


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
                     db="datanator-test", 
                     skip=0):
        """Processing documents and transform.

        Args:
            col(:obj:`str`): Name of the source collection.
            db(:obj:`Obj`): Name of database.
        """
        query = {}
        projection = {"_id": 0}
        docs = self.client[db][col].find(filter=query, projection=projection)
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
            elif col == "rna_halflife_new":
                obs = self.build_rna_observation(doc)
                for ob in obs:
                    _filter = {"$and": [{"identifier": ob["identifier"]},
                                        {"source": {"$elemMatch": ob["source"][0]}},
                                        {"environment": ob["environment"]}]}
                    self.update_observation(ob,
                                            ob["source"][0],
                                            db=self.db,
                                            query=_filter)
            elif col == "rna_modification":
                if doc.get("amino_acid") is None:
                    continue
                entity = self.build_rna_modification_entity(doc)
                obs = self.build_rna_modification_observation(doc)
                self.update_entity(entity,
                                   entity["identifiers"][0])
                for ob in obs:
                    query = {"$and": [{"identifier": ob["identifier"]},
                                       {"genotype.taxon.ncbi_taxonomy_id": ob["genotype"]["taxon"]["ncbi_taxonomy_id"]}]}
                    self.update_observation(ob,
                                            ob["source"][0],
                                            query=query)  

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
        mod = self.client["datanator-test"]["protein_modifications"].find({"uniprot_id": obj.get("uniprot_id")})            
        if mod is not None:
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

    def build_rna_observation(self, obj):
        """Build RNA observation object from rna_haflife_new collection.

        Args:
            obj(:obj:`Obj`): object to be transformed.

        Return:
            (:obj:`list` of :obj:`Obj`)
        """
        result = []
        pattern = r"^gm[\d]+"
        schema_version = "2.0"
        for life in obj.get("halflives"):
            entity = {"schema_version": schema_version}
            value = {"type": "rna_halflife"}
            environment = {}
            genotype = {}
            source = []            
            related = []
            related.append({"namespace": "uniprot_id",
                            "value": obj.get("uniprot_id")})
            related.append({"namespace": "kegg_orthology_id",
                            "value": obj.get("ko_number")})
            entity["type"] = "RNA"
            entity["name"] = obj.get("protein_names")[0]
            entity["identifiers"] = [{"namespace": "uniprot_id",
                                      "value": obj.get("uniprot_id")}]
            values = []
            if life.get("transcript_size") is not None:
                entity["structures"] = []
            if life.get("values") is None:
                for key, val in life.items():
                    if key == "unit":
                        value["units"] = val
                    elif key == "halflife":
                        value["value"] = float(val)
                    elif key == "ncbi_taxonomy_id": # protein entity includes taxon info
                        continue
                    elif key == "species_name":
                        continue
                    elif key == "reference":
                        for ref in val:
                            source.append({"namespace": list(ref.keys())[0],
                                        "value": list(ref.values())[0]})
                    elif key == "gene_name":
                        related.append({"namespace": key,
                                        "value": val})
                    elif key == "gene_symbol":
                        related.append({"namespace": key,
                                        "value": val})
                    elif key == "systematic_name": 
                        entity["identifiers"].append({"namespace": key,
                                                    "value": val})
                    elif key == "accession_id": 
                        if isinstance(val, str):
                            entity["identifiers"].append({"namespace": key,
                                                        "value": val})
                        elif isinstance(val, list):
                            for _id in val:
                                entity["identifiers"].append({"namespace": key,
                                                            "value": _id})                            
                    elif key == "variation_coefficient":
                        value["uncertainty"] = float(val)
                    elif key == "growth_medium":
                        environment["media"] = val
                    elif key == "ordered_locus_name":
                        entity["identifiers"].append({"namespace": key,
                                                    "value": val})
                    elif key == "doubling_time":
                        environment[key] = val
                    elif key == "r_squared":
                        value["uncertainty"] = float(val)
                    elif key == "standard_error":
                        value["uncertainty"] = float(val)
                    elif key == "transcript_size":
                        entity["structures"].append({"format": key,
                                                    "value": str(val)})
                    elif key == "cds_size":
                        entity["structures"].append({"format": key,
                                                    "value": str(val)})
                    elif key == "intron_size":
                        entity["structures"].append({"format": key,
                                                    "value": str(val)})
                    elif key == "genomic_size":
                        entity["structures"].append({"format": key,
                                                    "value": str(val)})
                    elif key == "intron_count":
                        entity["structures"].append({"format": key,
                                                    "value": str(val)})
                    elif key == "std":
                        value["uncertainty"] = float(val)
                    elif key == "ar_cog":
                        entity["identifiers"].append({"namespace": key,
                                                    "value": val})
                    elif key == "cog":
                        entity["identifiers"].append({"namespace": key,
                                                    "value": val})
                    elif key == "quantification_method":
                        environment["condition"] = "Quantified via {}.".format(val)
                values.append(value)
                result.append({"entity": entity,
                            "genotype": genotype,
                            "values": values,
                            "environment": environment,
                            "source": source,
                            "identifier": {"namespace": "uniprot_id",
                                            "value": obj.get("uniprot_id")},
                            "schema_version": schema_version,
                            "related": related})
            else:  #gmxxxxx
                val = life.get("values")  
                entity["identifiers"].append({"namespace": "accession_id",
                                              "value": life.get("accession_id")[:-1]})
                related.append({"namespace": "gene_symbol",
                                "value": life.get("gene_symbol")})
                source = life.get("reference")                                          
                for o in val:
                    value = {}
                    environment = {}
                    genotype = {}
                    for k, v in o.items():
                        if re.search(pattern, k):
                            value["value"] = v * 3600.0
                            genotype["cellLine"] = k
                            value["units"] = "s"
                        elif k == "note":
                            environment["condition"] = v
                        elif k == "biological_replicates":
                            environment["replicate"] = v       
                    result.append({"entity": entity,
                                "genotype": genotype,
                                "values": [value],
                                "environment": environment,
                                "source": source,
                                "identifier": {"namespace": "uniprot_id",
                                                "value": obj.get("uniprot_id")},
                                "schema_version": schema_version,
                                "related": related}) 
        return result                   

    def build_rna_modification_entity(self, obj):
        """build entity objects from rna_modification collection

        Args:
            obj (:obj:`Obj`): Object from which entity will be built.

        Return:
            (:obj:`Obj`)
        """
        entity = {"type": "RNA", "schema_version": "2.0",
                  "name": obj["definition"], "related": [],
                  "identifiers": []}
        entity["related"].append({"namespace": "amino_acid",
                                  "value": obj["amino_acid"]})
        entity["related"].append({"namespace": "aa_code",
                                  "value": obj["aa_code"]})                    
        entity["related"].append({"namespace": "aa_name",
                                  "value": obj["aa_name"]})
        entity["related"].append({"namespace": "kegg_orthology_id",
                                  "value": obj["kegg_orthology_id"]})                                  
        entity["related"].append({"namespace": "kegg_orthology_name",
                                  "value": obj["kegg_gene_name"]})
        entity["related"].append({"namespace": "kegg_pathway_id",
                                  "value": obj["kegg_pathway_id"]})
        entity["related"].append({"namespace": "kegg_pathway_name",
                                  "value": obj["kegg_pathway_name"]})
        entity["identifiers"].append({"namespace": "amino_acid",
                                      "value": obj["amino_acid"]})
        return entity

    def build_rna_modification_observation(self, obj):
        """Build observation objects from rna_modification collection.

        Args:
            obj (:obj:`Obj`): object from which observation objects will be built.

        Return:
            (:obj:`list` of :obj:`Obj`)
        """
        result = []        
        entity = self.build_rna_modification_entity(obj)
        for mod in obj["modifications"]:
            genotype = {}
            values = []
            entity["stuctures"] = []
            identifier = {"namespace": "sequence_bpforms",
                          "value": mod["sequence_bpforms"]}
            if mod["ncbi_taxonomy_id"] is None:
                genotype["taxon"] = self.build_taxon_object(mod["organism"], _format="tax_name")
            else:    
                genotype["taxon"] = self.build_taxon_object(mod["ncbi_taxonomy_id"])
            genotype["organ"] = mod["organellum"]
            entity["stuctures"].append({"format": "sequence_modomics",
                                        "value": mod["sequence_modomics"]})
            entity["stuctures"].append({"format": "sequence_bpforms",
                                        "value": mod["sequence_bpforms"]})    
            entity["stuctures"].append({"format": "sequence_iupac",
                                        "value": mod["sequence_iupac"]})
            values.append({"type": "length",
                           "value": mod["length"]}) 
            values.append({"type": "anticodon",
                           "value": mod["anticodon"]})  
            values.append({"type": "number_of_modifications",
                           "value": mod["number_of_modifications"]})
            values.append({"type": "number_of_modified_a",
                           "value": mod["number_of_modified_a"]})
            values.append({"type": "number_of_modified_u",
                           "value": mod["number_of_modified_u"]})
            values.append({"type": "number_of_modified_c",
                           "value": mod["number_of_modified_c"]})
            values.append({"type": "number_of_modified_g",
                           "value": mod["number_of_modified_g"]})
            values.append({"type": "formula",
                           "value": mod["formula"]})
            values.append({"type": "molecular_weight",
                           "value": mod["molecular_weight"]})
            values.append({"type": "charge",
                           "value": mod["charge"]})
            values.append({"type": "canonical_formula",
                           "value": mod["canonical_formula"]})
            values.append({"type": "canonical_molecular_weight",
                           "value": mod["canonical_molecular_weight"]})
            values.append({"type": "canonical_charge",
                           "value": mod["canonical_charge"]})
            values.append({"type": "bpforms_errors",
                           "value": mod["bpforms_errors"]})
            source = [{"namespace": "doi",
                       "value": mod["reference"]["doi"]}]
            result.append({"entity": entity,
                           "values": values,
                           "genotype": genotype,
                           "source": source,
                           "identifier": identifier,
                           "schema_version": "2.0"})
        return result

    def build_sabio_entity(self, obj):
        """Build entity objects from sabio_rk_old collection.

        Args:
            obj (:obj:`Obj`): object from which entity objects will be built.

        Return:
            (:obj:`Obj`)
        """
        entity = {"type": "reaction"}
        structures = [{"format": "reaction_equation", "enzyme": {"subunit": []}}]        
        identifiers = []
        related = []
        source = []
        related.append({"namespace": "kinlaw_id",
                        "value": obj.get("kinlaw_id")})
        resources = obj.get("resource", [{}])
        for r in resources:
            for key, val in r.items():
                if key == "pubmed":
                    source.append({"namespace": key,
                                   "value": val})
                elif key == "ec-code":
                    related.append({"namespace": "ec_code",
                                    "value": val})
                elif key == "sabiork.reaction":
                    identifiers.append({"namespace": "sabio_rxn_id",
                                        "value": val})
                else:
                    with open("edgecase.txt", "a+") as f:
                        f.write("Kinetic law {} has more info in resources".format(obj["kinlaw_id"]))
        enzymes_list = obj["enzymes"]
        for e in enzymes_list:
            if e.get("subunit") is not None:
                for i, su in enumerate(e.get("subunit")):
                    structures[0]["enzyme"]["subunit"].append({"namespace": "uniprot_id", 
                                                               "value": su.get("uniprot_id"),
                                                               "coefficient": su.get("subunit_coefficient")})
            elif e.get("enzyme") is not None:
                for i, en in enumerate(e.get("enzyme")):
                    structures[0]["enzyme"]["name"] = e.get("enzyme")[i]["enzyme_name"]
                    structures[0]["enzyme"]["id"] = e.get("enzyme")[i]["enzyme_id"]
            elif e.get("compartment") is not None:
                for i, en in enumerate(e.get("compartment")):
                    structures[0]["enzyme"]["compartment"] = e.get("compartment")[i]["compartment_name"]
        parameters = obj["parameter"]
        for p in parameters:
            pass

def main():
    pass

if __name__ == "__main__":
    main()