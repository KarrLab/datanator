from bioservices import *
from datanator_query_python.util import mongo_util
import requests
import simplejson


class XRef(mongo_util.MongoUtil):
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
        self.kegg = KEGG()
        self.uniprot = UniProt()
        self.ortho = self.client["datanator"]["orthodb"]
        self.uniprot_col = self.client["datantor-test"]["uniprot"]

    def get_kegg_rxn(self, _id):
        """Use bioservice to request kegg reaction information.

        Args:
            _id(:obj:`str`): Kegg reaction id.

        Return:
            (:obj:`tuple` of :obj:`list`): substrates and products in inchikey lists
        """
        agg = self.kegg.get(_id)
        obj = self.kegg.parse(agg)
        eqn = obj["EQUATION"]
        defi = obj["DEFINITION"]
        substrates = []
        products = []
        for i, (c, d) in enumerate(zip(eqn.split(" <=> "), defi.split(" <=> "))):
            if i == 0:  # substrates string
                # kegg compound id to ChEBI id
                tmp = [KEGG().parse(KEGG().get(x))["DBLINKS"]["ChEBI"] for x in c.split(" + ")]
                # chebi id to inchikey
                for x in tmp:
                    try:
                        substrates.append(ChEBI().getCompleteEntity("CHEBI:"+x[0:5]).inchiKey)
                    except AttributeError:
                        substrates.append(ChEBI().getCompleteEntity("CHEBI:"+x[0:5]).smiles)
                substrates_name = [x for x in d.split(" + ")]
            elif i == 1: # products string
                tmp = [KEGG().parse(KEGG().get(x))["DBLINKS"]["ChEBI"] for x in c.split(" + ")]
                # chebi id to inchikey
                for x in tmp:
                    try:
                        products.append(ChEBI().getCompleteEntity("CHEBI:"+x[0:5]).inchiKey)
                    except AttributeError:
                        products.append(ChEBI().getCompleteEntity("CHEBI:"+x[0:5]).smiles)
                products_name = [x for x in d.split(" + ")]
        return substrates, substrates_name, products, products_name

    def uniprot_id_to_orthodb(self, _id, cache={}):
        """Convert uniprot id to orthodb group.

        Args:
            (:obj:`str`): Uniprot ID.
            (:obj:`Obj`, optional): existing x-ref records.

        Return:
            (:obj:`Obj`): {"orthodb_id": ... "orthodb_name": ...}
        """
        if cache.get(_id) is not None:
            return cache.get(_id), cache
        else:
            try:
                u = self.uniprot.search("id:{}".format(_id),
                                        columns="database(OrthoDB)")
            except TypeError:
                return {"orthodb_id": None,
                        "orthodb_name": None}, cache                
            _list = u.split()
            if len(_list) == 2:
                return {"orthodb_id": None,
                        "orthodb_name": None}, cache
            orthodb = _list[2][:-1]
            try:
                name = self.ortho.find_one(filter={"orthodb_id": orthodb})["orthodb_name"]
            except TypeError:
                tmp = requests.get("https://dev.orthodb.org/group?id={}".format(orthodb)).json()["data"]
                if isinstance(tmp, list):
                    return {"orthodb_id": None,
                            "orthodb_name": None}, cache
                else:
                    name = tmp["name"]
            obj = {"orthodb_id": orthodb,
                   "orthodb_name": name}
            cache[_id] = obj
            return obj, cache
        
    def uniprot_to_interpro(self, _id):
        """Convert uniprot ID to interpro id.

        Args:
            (:obj:`str`): Uniprot ID.
            (:obj:`Obj`, optional): existing x-ref records.

        Return:
            (:obj:`list`): List of InterPro IDs.
        """
        u = self.uniprot.search("id:{}".format(_id),
                                columns="database(InterPro)")
        tmp = u.strip().split("\n")
        if len(tmp) != 1:
            return tmp[1].split(";")[:-1]
        else:
            return []

    def name_level(self,
                   name,
                   level=2759,
                   limit=1):
        """Given protein name and taxon level, return orthodb group ID.

        Args:
            name (:obj:`str`): Name of the protein.
            level (:obj:`int`, optional): Taxon level. Defaults to 2759.
            limit (:obj:`int`, optional): Number of data to return. Defaults to 1.

        Return:
            (:obj:`str`): Orthodb Group ID.
        """
        tmp = requests.get("https://dev.orthodb.org/search?query={}&limit={}&level={}".format(name, limit, level)).json()["data"]
        if len(tmp) == 0:
            return ""
        else:
            return tmp[0]