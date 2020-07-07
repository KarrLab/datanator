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
                 max_entries=float('inf'),
                 verbose=True):
        super().__init__(MongoDB=MongoDB,
                         db=db,
                         username=username,
                         password=password)
        self.max_entries = max_entries
        self.des = self.db_obj[des_col]
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
        docs = self.client["datanator-test"]["metabolite_concentrations"].find(filter=query, projection=projection)
        pass

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
        identifiers.append({"namespace": "uniprot",
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
        structures = {"format": "canonical_sequence_iubmb",
                      "value": obj.get("canonical_sequence")}
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