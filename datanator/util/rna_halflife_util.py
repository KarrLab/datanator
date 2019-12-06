from datanator_query_python.query import query_uniprot
from datanator.data_source import uniprot_nosql
from datanator_query_python.util import mongo_util
import math
import pandas as pd


class RnaHLUtil(mongo_util.MongoUtil):

    def __init__(self, server=None, username=None, password=None, src_db=None,
                des_db=None, protein_col=None, rna_col=None, authDB='admin', readPreference=None,
                max_entries=float('inf'), verbose=False):
        super().__init__(MongoDB=server, db=des_db, verbose=verbose, max_entries=max_entries,
        username=username, password=password, authSource=authDB, readPreference=readPreference)
        _, _, self.rna_hl_collection = self.con_db(rna_col)
        self.max_entries = max_entries
        self.verbose = verbose
        self.uniprot_query_manager = query_uniprot.QueryUniprot(username=username, password=password,
                                                                server=server, authSource=authDB,
                                                                database=src_db, collection_str=protein_col)
        self.uniprot_collection_manager = uniprot_nosql.UniprotNoSQL(MongoDB=server, db=des_db, verbose=True,
        username=username, password=password, authSource=authDB, collection_str=protein_col)

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

    def make_df(self, url, sheet_name, header=None, names=None, usecols=None,
                skiprows=None, nrows=None, na_values=None):
        """Read online excel file as dataframe

        Args:
            url (:obj:`str`): excel file url
            sheet_name (:obj:`str`): name of sheet in xlsx
            header (:obj:`int`): Row (0-indexed) to use for the column labels of the parsed DataFrame.
            names (:obj:`list`): list of column names to use
            usecols (:obj:`int` or :obj:`list` or :obj:`str`): Return a subset of the columns.
            nrows (:obj:`int`): number of rows to parse. Defaults to None.

        Returns:
            (:obj:`pandas.DataFrame`): xlsx transformed to pandas.DataFrame
        """
        data = pd.read_excel(url, sheet_name=sheet_name, nrows=nrows, header=header,
                            names=names, usecols=usecols, skiprows=skiprows, na_values=na_values)
        return data

    def fill_uniprot_with_df(self, df, ordered_locus_name):
        """Fill uniprot colleciton with ordered_locus_name
        from excel sheet
        
        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be inserted into uniprot collection.
            Assuming df conforms to the schemas required by load_uniprot function in uniprot.py
            ordered_locus_name (:obj:`str`): name of column that stores ordered locus name information.
        """
        # df = self.make_df(url, sheet_name, usecols=usercols, skiprows=skiprows,
        # names=names)
        row_count = len(df.index)
        for index, row in df.iterrows():
            if index == self.max_entries:
                break
            if index % 10 == 0 and self.verbose:
                print("Inserting locus {} out of {} into uniprot collection.".format(index, row_count))
            oln = row[ordered_locus_name]
            self.fill_uniprot_by_oln(oln)    