from pathlib import Path
from datanator_query_python.util import mongo_util
import pandas as pd
from bioservices import *


class InsRxn(mongo_util.MongoUtil):
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
        self.kegg = KEGG()

    def parse_docs(self,
                   _dir="./docs/"):
        """Parse database with xls files from
        http://www.cecafdb.org/, stored in ./docs/

        Args:
            _dir(:obj:`str`): Directory in which xls files are stored.
        """
        results = []
        paths = Path(_dir).glob('**/*.xls')
        for i, path in enumerate(paths):
            if i == self.max_entries:
                break
            df = pd.read_excel(str(path),
                              sheet_name="sheet one",
                              skiprows=[0,1,2,3,4,7,15,16,
                                        17,18,19,20,21],
                              usecols="A:F")
            for i, row in df.iterrows():
                print(row.to_json())

    def parse_eqn(self,
                  eqn):
        """Parse rxn equations into arrays of names of 
        substrates and products.

        Args:
            eqn(:obj:`str`): Rxn equation. e.g. Succinate-mit + CoA <==> Succinyl-CoA

        Return
            (:obj:`tuple` of :obj:`list`)
        """
        pass

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
        substrates = []
        products = []
        for i, c in enumerate(eqn.split(" <=> ")):
            if i == 0:
                # kegg compound id to ChEBI id
                tmp = [KEGG().parse(KEGG().get(x))["DBLINKS"]["ChEBI"] for x in c.split(" + ")]
                # chebi id to inchikey
                substrates = [ChEBI().getCompleteEntity("CHEBI:"+x[0:5]).inchiKey for x in tmp]
            elif i == 1:
                tmp = [KEGG().parse(KEGG().get(x))["DBLINKS"]["ChEBI"] for x in c.split(" + ")]
                # chebi id to inchikey
                products = [ChEBI().getCompleteEntity("CHEBI:"+x[0:5]).inchiKey for x in tmp]
        return substrates, products