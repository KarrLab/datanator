from datanator.util import rna_halflife_util
from datanator_query_python.query import query_uniprot_org
import datetime
import datanator.config.core
from pymongo.collation import Collation, CollationStrength
import tempfile
import shutil


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
        max_entries=max_entries, verbose=verbose, cache_dir=cache_dir)
        self.collation = Collation('en', strength=CollationStrength.SECONDARY)
        self.max_entries = max_entries
        self.verbose = verbose

    def fill_rna_halflife(self, df, start=2):
        """Fill rna_halflife collection with information parsed
        from the publication
        
        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be loaded.
            start (:obj:`int`, optional): Starting row in df. Defaults to 0.
        """
        row_count = len(df.index)
        for i, row in df.iloc[start:].iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            probeset_id = row.get('probeset_id')
            gene_symbol = row.get('gene_symbol')
            accession_id = row.get('accession_id')
            halflife_0 = {'gm07029': row.get('gm07029_a1'), 'biological_replicates': 'a1', 'note': 'independent cell cultures for the same cell line', 'unit': 'hr'}
            halflife_1 = {'gm07029': row.get('gm07029_a2'), 'biological_replicates': 'a2', 'note': 'independent cell cultures for the same cell line', 'unit': 'hr'}
            halflife_2 = {'gm07029': row.get('gm07029_a3'), 'biological_replicates': 'a3', 'note': 'independent cell cultures for the same cell line', 'unit': 'hr'}
            halflife_3 = {'gm10835': row.get('gm10835_a1'), 'biological_replicates': 'a1', 'note': 'independent cell cultures for the same cell line', 'unit': 'hr'}
            halflife_4 = {'gm10835': row.get('gm10835_a2'), 'biological_replicates': 'a2', 'note': 'independent cell cultures for the same cell line', 'unit': 'hr'}
            halflife_5 = {'gm10835': row.get('gm10835_a3'), 'biological_replicates': 'a3', 'note': 'independent cell cultures for the same cell line', 'unit': 'hr'}
            halflife_6 = {'gm12813': row.get('gm12813_a1'), 'biological_replicates': 'a1', 'note': 'independent cell cultures for the same cell line', 'unit': 'hr'}
            halflife_7 = {'gm12813': row.get('gm12813_a1_dup'), 'technical_replicates': 'a1', 'note': 'separate RNA aliquots from the same cell culture', 'unit': 'hr'}
            halflife_8 = {'gm12813': row.get('gm12813_a2'), 'biological_replicates': 'a2', 'note': 'independent cell cultures for the same cell line', 'unit': 'hr'}
            halflife_9 = {'gm12813': row.get('gm12813_a2_dup'), 'technical_replicates': 'a2', 'note': 'separate RNA aliquots from the same cell culture', 'unit': 'hr'}
            # halflife_10 = {'gm12813': row.get('gm12813_a3'), 'biological_replicates': 'a3', 'note': 'independent cell cultures for the same cell line', 'unit': 'hr'}
            halflife_11 = {'gm12813': row.get('gm12813_a3_dup'), 'technical_replicates': 'a3', 'note': 'separate RNA aliquots from the same cell culture', 'unit': 'hr'}
            halflife_12 = {'gm07019': row.get('gm07019'), 'unit': 'hr'}
            halflife_13 = {'gm12812': row.get('gm12812'), 'unit': 'hr'}
            halflife_14 = {'gm12814': row.get('gm12814'), 'unit': 'hr'}
            halflife_15 = {'gm12815': row.get('gm12815'), 'unit': 'hr'}
            anova_p_value_3 = row.get('anova_3')
            fdr_3 = row.get('fdr_3')
            anova_p_value_7 = row.get('anova_7')
            fdr_7 = row.get('fdr_7')
            uniprot_manager = query_uniprot_org.QueryUniprotOrg(accession_id + ' homo sapiens')
            uniprot_id = uniprot_manager.get_uniprot_id()
            ko_number = uniprot_manager.get_kegg_ortholog()
            protein_names = uniprot_manager.get_protein_name()
            species = 'Homo sapiens'
            ncbi_taxonomy_id = 9606
            reference = [{'doi': '10.1038/srep01318'}]

            values = [halflife_0, halflife_1, halflife_2, halflife_3,
                        halflife_4, halflife_5, halflife_6, halflife_7,
                        halflife_8, halflife_9, halflife_11,
                        halflife_12, halflife_13, halflife_14, halflife_15]
            obj = {'accession_id': accession_id,
                    'probeset_id': probeset_id,
                    'values': values,
                    'anova_3': anova_p_value_3,
                    'anova_7': anova_p_value_7,
                    'false_discovery_rate_3': fdr_3,
                    'false_discovery_rate_7': fdr_7,
                    'species': species,
                    'ncbi_taxonomy_id': ncbi_taxonomy_id,
                    'reference': reference,
                    'gene_symbol': gene_symbol}

            if uniprot_id is not None:
                self.rna_hl_collection.update_one({'uniprot_id': uniprot_id},
                                                  {'$addToSet': {'halflives': obj},
                                                   '$set': {'ko_number': ko_number,
                                                            'protein_names': protein_names}},
                                                  collation=self.collation, upsert=True)
            else:
                self.rna_hl_collection.update_one({'identifier': accession_id},
                                                  {'$addToSet': {'halflives': obj}},
                                                  collation=self.collation, upsert=True)                


def main():
    src_db = 'datanator'
    des_db = 'datanator'
    rna_col = 'rna_halflife'
    protein_col = 'uniprot'
    cache_dir = tempfile.mkdtemp()
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    src = Halflife(server=server, src_db=src_db,
        protein_col=protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True, max_entries=float('inf'),
        des_db=des_db, rna_col=rna_col, cache_dir=cache_dir)
    url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fsrep01318/MediaObjects/41598_2013_BFsrep01318_MOESM2_ESM.xls'
    names = ['probeset_id', 'gene_symbol', 'accession_id', 'gm07029_a1', 'gm07029_a2', 'gm07029_a3',
            'gm10835_a1', 'gm10835_a2', 'gm10835_a3', 'gm12813_a1', 'gm12813_a1_dup', 'gm12813_a2', 'gm12813_a2_dup',
            'gm12813_a3_dup', 'gm07019', 'gm12812', 'gm12814', 'gm12815', 'irrelavent', 'anova_3',
            'fdr_3', 'anova_7', 'fdr_7']
    df_s1 = src.make_df(url, 'S3', names=names, usecols='A:W', file_type='xls')
    src.fill_rna_halflife(df_s1)
    shutil.rmtree(cache_dir)

if __name__ == '__main__':
    main()