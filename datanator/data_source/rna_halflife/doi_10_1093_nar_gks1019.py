import pandas as pd
from datanator.util import rna_halflife_util
import datetime
import datanator.config.core
import datetime
from pymongo.collation import Collation, CollationStrength


class Halflife(rna_halflife_util.RnaHLUtil):

    def __init__(self, cache_dir=None, server=None, src_db=None, protein_col=None,
                authDB=None, readPreference=None, username=None, password=None,
                verbose=None, max_entries=None, des_db=None, rna_col=None):
        """Init
        
        Args:
            cache_dir (:obj:`str`, optional): Cache directory for logs. Defaults to None.
            server (:obj:`str`, optional): MongoDB server address. Defaults to None.
            db (:obj:`str`, optional): Database where initial uniprot collection resides. Defaults to None.
            collection_str (:obj:`str`, optional): name of collection. Defaults to None.
            authDB (:obj:`str`, optional): MongoDB authentication database. Defaults to None.
            readPreference (:obj:`str`, optional): MongoDB read preference. Defaults to None.
            username (:obj:`str`, optional): MongoDB username. Defaults to None.
            password (:obj:`str`, optional): MongoDB password. Defaults to None.
            verbose (:obj:`bool`, optional): Wheter to display verbose messages. Defaults to None.
            max_entries (:obj:`int`, optional): Number of records to be processed. Defaults to None.
            uniprot_col_db (:obj:`int`, optional): Database to which new uniprot records will be inserted. Defaults to None.
        """
        super().__init__(server=server, username=username, password=password, src_db=src_db,
        des_db=des_db, protein_col=protein_col, rna_col=rna_col, authDB=authDB, readPreference=readPreference,
        max_entries=max_entries, verbose=verbose)
        self.collation = Collation('en', strength=CollationStrength.SECONDARY)
        self.max_entries = max_entries
        self.verbose = verbose

    def fill_uniprot(self, url, sheet_name, usercols='B:D', skiprows=[0,1,2],
                    insertion=True):
        """Fill uniprot colleciton with ordered_locus_name
        from excel sheet
        
        Args:
            url (:obj:`str`): URL for Excel sheet.
            sheet_name (:obj:`str`): sheet name within Excel.
            usecols (:obj:`int` or :obj:`list` or :obj:`str`): Return a subset of the columns.
            skiprows (:obj:`list`): rows to skip (0-indexed)
            insertion (:obj:`bool`): whether to insert new records to uniprot collection.

        Return:
            (:obj:`pandas.DataFrame`): Dataframe
        """
        df = self.make_df(url, sheet_name, usecols=usercols, skiprows=skiprows,
        names=['ordered_locus_name', 'half_life', 'r_squared'])
        row_count = len(df.index)
        if insertion:
            for index, row in df.iterrows():
                if index == self.max_entries:
                    break
                if index % 10 == 0 and self.verbose:
                    print("Inserting locus {} out of {} into uniprot collection.".format(index, row_count))
                oln = row['ordered_locus_name']
                self.fill_uniprot_by_oln(oln)
        return df

    def fill_rna_halflife(self, df, species):
        """load data into rna_halflife collection
        
        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be loaded into the database
            species (:obj:`list`): species name and ncbi_id
        """
        row_count = len(df.index)
        for i, row in df.iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            halflives = {}
            oln = row['ordered_locus_name']
            halflives['halflife'] = row['half_life'] * 60
            halflives['r_squared'] = row['r_squared']
            halflives['unit'] = 's'
            halflives['reference'] = [{'doi': '10.1093/nar/gks1019', 'pubmed_id': '23125364'}]
            halflives['growth_medium'] = 'Middlebrook 7H9 with the ADC supplement (Difco) and 0.05% Tween80, at 37 degree celcius.'
            halflives['ordered_locus_name'] = oln
            halflives['species'] = species[0]
            halflives['ncbi_taxonomy_id'] = species[1]
            gene_name, protein_name = self.uniprot_query_manager.get_gene_protein_name_by_oln(oln)
            if gene_name is not None: # record exists in uniprot collection with gene_name
                self.rna_hl_collection.update_one({'gene_name': gene_name},
                                                  {'$set': {'modified': datetime.datetime.utcnow()},
                                                   '$addToSet': {'halflives': halflives,
                                                                'protein_synonyms': protein_name}}, 
                                                   collation=self.collation, upsert=True)
            elif (gene_name is None and protein_name is not None and 
                  protein_name != 'Uncharacterized protein'): # record exists in uniprot collection with non-filler protein_name
                self.rna_hl_collection.update_one({'protein_name': protein_name},
                                                  {'$set': {'modified': datetime.datetime.utcnow(),
                                                            'gene_name': gene_name},
                                                   '$addToSet': {'halflives': halflives,
                                                                'protein_synonyms': protein_name}}, 
                                                   collation=self.collation, upsert=True)
            else:
                query = {'halflives.ordered_locus_name': oln}
                doc = self.rna_hl_collection.find_one(filter=query, collation=self.collation)
                if doc is not None:
                    self.rna_hl_collection.update_one({'halflives.ordered_locus_name': oln},
                                                    {'$set': {'modified': datetime.datetime.utcnow(),
                                                              'gene_name': gene_name},
                                                    '$addToSet': {'halflives': halflives,
                                                                  'protein_synonyms': protein_name}}, 
                                                    collation=self.collation, upsert=True)
                else:
                    doc = {'halflives': [halflives], 'modified': datetime.datetime.utcnow(),
                            'gene_name': gene_name, 'protein_name': protein_name}
                    self.rna_hl_collection.insert_one(doc)


def main():
    src_db = 'datanator'
    des_db = 'datanator'
    rna_col = 'rna_halflife'
    protein_col = 'uniprot'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    src = Halflife(server=server, src_db=src_db,
        protein_col=protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True, max_entries=float('inf'),
        des_db=des_db, rna_col=rna_col)
    url = 'https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/41/1/10.1093/nar/gks1019/2/gks1019-nar-00676-a-2012-File003.xlsx?Expires=1578425844&Signature=ZRFUxLdn4-vaBt5gQci~0o56KqyR9nJj9i32ig5X6YcfqiJeV3obEq8leHGdDxx6w~KABgewiQ66HTB7gmuG~2GL-YgxPKYSjt17WrYMkc-0ibw6TMlTvWZZfvw-lPe~wvpmVfNEXnTbP7jHyNLu9jeJ6yhoXvgIyQtzA5PbEI1fyXEgeZzOKMltmITqL3g3APsPsagCTC66rwrBT23Aghh6D314uilT2DZHCc68MH2nyV~qAhFqIQiOj-7VTEKqkDPvPYvuE2KNKXdvW23gk100YV~58ozbt8ijRz5Gr5gPtE~f1Ab5l260EIbWHJNabMRleInJQqUIDPFN4C38PQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'
    # df = src.fill_uniprot(url, 'Supplementary Table 1', insertion=False)
    # src.fill_rna_halflife(df, ['Mycobacterium tuberculosis H37Rv', 83332])
    df = src.fill_uniprot(url, 'Supplementary Table 2', skiprows=list(range(0,6)))
    src.fill_rna_halflife(df, ['Mycolicibacterium smegmatis MC2 155', 246196])

if __name__ == '__main__':
    main()