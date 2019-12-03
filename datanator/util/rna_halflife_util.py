from datanator_query_python.query import query_uniprot
from datanator.data_source import uniprot_nosql
from datanator_query_python.util import mongo_util
import math
import pandas as pd


class RnaHLUtil(mongo_util.MongoUtil):

    def __init__(self, server=None, username=None, password=None, src_db=None,
                des_db=None, protein_col=None, rna_col=None, authDB='admin', readPreference=None,
                max_entries=float('inf'), verbose=False):
        super().__init__(self, MongoDB=server, db=des_db, verbose=verbose, max_entries=max_entries,
        username=username, password=password, authSource=authDB, readPreference=readPreference)
        self.max_entries = max_entries
        self.uniprot_query_manager = query_uniprot.QueryUniprot(username=username, password=password,
                                                                server=server, authSource=authDB,
                                                                database=src_db, collection_str=protein_col)
        self.uniprot_collection_manager = uniprot_nosql.UniprotNoSQL(MongoDB=server, db=des_db, verbose=True,
        username=username, password=password, authSource=authDB, collection_str=protein_col)

    def fill_uniprot_by_oln(self, oln):
        """Fill uniprot collection using ordered locus name
        
        Args:
            oln (:obj:`str`): Ordered locus name
        """
        gene_name, protein_name = self.uniprot_query_manager.get_gene_protein_name_by_oln(oln)
        if gene_name is None and protein_name is None: # no such entry in uniprot collection
            self.uniprot_collection_manager.load_uniprot(query=True, msg=oln)
        else:
            return

    def uniprot_names(self, results, count):
        """Extract protein_name and gene_name from returned
        tuple of uniprot query function
        
        Args:
            results (:obj:`Iter`): pymongo cursor object.
            count (:obj:`int`): Number of documents found.

        Return:
            (:obj:`tuple` of :obj:`str`): gene_name and protein_name
        """
        if count == 0:
            return '', ''
        else:
            for result in results:
                gene_name = result['gene_name']
                protein_name = result['protein_name']
                return gene_name, protein_name

    def make_df(self, url, sheet_name, header=None, names=None, usecols=None,
                skiprows=None):
        """Read online excel file as dataframe

        Args:
            url (:obj:`str`): excel file url
            sheet_name (:obj:`str`): name of sheet in xlsx
            header (:obj:`int`): Row (0-indexed) to use for the column labels of the parsed DataFrame.
            names (:obj:`list`): list of column names to use
            usecols (:obj:`int` or :obj:`list` or :obj:`str`): Return a subset of the columns.

        Returns:
            (:obj:`pandas.DataFrame`): xlsx transformed to pandas.DataFrame
        """
        if self.max_entries == float('inf'):
            nrows = None
        else:
            nrows = self.max_entries
        data = pd.read_excel(url, sheet_name=sheet_name, nrows=nrows, header=header,
                            names=names, usecols=usecols, skiprows=skiprows)
        return data
    