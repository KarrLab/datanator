from bioservices import *
import requests


class XRef:
    def __init__(self):
        self.kegg = KEGG()
        self.uniprot = UniProt()

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
            u = self.uniprot.search("id:{}".format(_id),
                                    columns="database(OrthoDB),database(KO)")
            orthodb = u.split()[4][:-1]
            r = requests.get("https://dev.orthodb.org/group?id={}".format(orthodb))
            name = r.json()["data"]["name"]
            obj = {"orthodb_id": orthodb,
                   "orthodb_name": name}
            cache[_id] = obj
            return obj, cache