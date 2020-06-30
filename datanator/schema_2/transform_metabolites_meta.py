def build_entity(obj):
    """Build entity object from obj.

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
    entity["identifiers"].append({"namespace": "inchi",
                                  "value": obj.get("inchi")})
    entity["identifiers"].append({"namespace": "kegg_id",
                                  "value": obj.get("kegg_id")})
    entity["identifiers"].append({"namespace": "m2m_id",
                                  "value": obj.get("m2m_id")})
    entity["identifiers"].append({"namespace": "pubchem_compound_id",
                                  "value": obj.get("pubchem_compound_id")})
    entity["identifiers"].append({"namespace": "smiles",
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

def build_observation(obj):
    """Build observation objects from obj.

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
    obs = []