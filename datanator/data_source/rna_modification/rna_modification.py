import pandas as pd
import numpy as np
import datetime
from datanator_query_python.util import mongo_util
import datanator.config.core
from pymongo.collation import Collation, CollationStrength


class RNAMod(mongo_util.MongoUtil):

    def __init__(self, MongoDB=None, db=None, collection_str=None, username=None,
                 password=None, authSource='admin', readPreference='nearest',
                 verbose=True, max_entries=float('inf')):
        """Init        
        """
        super().__init__(MongoDB=MongoDB, db=db, username=username, password=password,
                         authSource=authSource, readPreference=readPreference)
        self.collation = Collation('en', strength=CollationStrength.SECONDARY)
        self.collection = self.db_obj[collection_str]
        self.taxon_collection = self.client['datanator']['taxon_tree']
        self.max_entries = max_entries
        self.verbose = verbose

    def fill_trna_primary(self, file_location, sheet_name='tRNA', start_row=[0],
                        use_columns='A:H', column_names=['Amino_acid', 'aa_Code', 'aa_Name',
                        'kegg_orthology_id', 'kegg_gene_Name', 'Definition', 'kegg_Pathway_id',
                        'kegg_Pathway_name']):
        """
        """
        print('here')
        df = pd.read_excel(file_location, sheet_name=sheet_name,
                         header=0, usecols=use_columns,
                         skiprows=start_row) 
        df.columns = [x.lower() for x in column_names]
        row_count = len(df.index)
        result = []
        for i, row in df.iterrows():
            if i == self.max_entries:
                break
            if i % 50 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            obj = row.to_dict()
            obj.pop('pathways')
            result.append(obj)
        self.collection.insert_many(result)      
    
    def fill_trna_collection(self, file_location, start_row=0,
                            column_names=[],
                            reference={'doi': '10.1093/nar/gkx1030'},
                            query='amino_acid'):
        """
        Fill collection collection_str.

        Args:
            sheet_name(:obj:`str`, optional): Name of sheet in excel.
            start_row (:obj:`int`, optional): Read from csv row. Defaults to 0.
            use_columns(:obj:`str`): Indicates comma separated list of Excel column letters and column ranges (e.g. “A:E” or “A,C,E:F”). Ranges are inclusive of both sides.
            column_names(:obj:`list` of :obj:`str`): Names of columns used.
            reference(:obj:`Obj`): reference information.
            unit(:obj:`str`): Unit used for Km.
        """
        df = pd.read_csv(file_location, skiprows=start_row, header=0,
                         engine='c', sep='\t').replace({np.nan:None})
        df.columns = [x.lower() for x in column_names]
        row_count = len(df.index)
        for i, row in df.iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            row['reference'] = reference
            row['last_modified'] = datetime.datetime.utcnow()
            try:
                row['ncbi_taxonomy_id'] = self.taxon_collection.find_one({'tax_name': row['organism']},
                                                                        projection={'tax_id': 1}, collation=self.collation).get('tax_id')
            except AttributeError:
                row['ncbi_taxonomy_id'] = None
            obj = row.to_dict()
            obj.pop('amino_acid', None)
            obj.pop('kegg_orthology_id')
            self.collection.update_one({query: row[query]},
                                        {'$addToSet': {'modifications': obj}},
                                        collation=self.collation, upsert=True)

    def fill_dummy_orthodb(self):
        docs = self.collection.find({})
        for i, doc in enumerate(docs):
            bak = doc.get("kegg_orthology_name")
            orthodb_id = doc.get("kegg_gene_name", bak)
            orthodb_name = orthodb_id
            self.collection.update_one({"_id": doc["_id"]},
                                        {"$set": {"orthodb_id": orthodb_id,
                                                  "orthodb_name": orthodb_name}})
        print("Done.")


import datanator.config.core
from pathlib import Path

def main():
    db = 'datanator'
    collection_str = 'rna_modification'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    manager = RNAMod(MongoDB=MongoDB, db=db, collection_str=collection_str,
                     username=username, password=password)

    # file_location = str(Path('~/karr_lab/datanator/docs/rna_modification/rna-ortholog-groups.xlsx').expanduser())
    # manager.fill_trna_primary(file_location)
    # column_names = ['kegg_orthology_id', 'kegg_orthology_name', 'definition',
    #                 'Mamallian_subunit', 'pathways', 'kingdom']
    # manager.fill_trna_primary(file_location, sheet_name='rRNA',
    #                           column_names=column_names, use_columns='A:F')

    # column_names = ['amino_acid', 'Anticodon', 'Organism', 'Organellum', 
    #                 'kegg_orthology_id', 'Sequence_MODOMICS', 'Sequence_BpForms',
    #                 'Sequence_IUPAC', 'Length',	'Number_of_modifications', 'Number_of_modified_A',
    #                 'Number_of_modified_C',	'Number_of_modified_G',	'Number_of_modified_U',	
    #                 'Formula', 'Molecular_weight', 'Charge', 'Canonical_formula',	
    #                 'Canonical_molecular_weight', 'Canonical_charge', 'Extra_formula',
    #                 'Extra_molecular_weight', 'Extra_charge', 'BpForms_errors']
    # file_location = Path('~/karr_lab/datanator/docs/rna_modification/modomics.trna.tsv').expanduser()
    # manager.fill_trna_collection(file_location, column_names=column_names,
    #                              start_row=None)

    # column_names = ['gen_bank', 'Organism', 'Organellum', 
    #                 'kegg_orthology_name', 'subunit', 'kegg_orthology_id', 
    #                 'sequence_modomics', 'Sequence_BpForms',
    #                 'Sequence_IUPAC', 'Length',	'Number_of_modifications', 'Number_of_modified_A',
    #                 'Number_of_modified_C',	'Number_of_modified_G',	'Number_of_modified_U',	
    #                 'Formula', 'Molecular_weight', 'Charge', 'Canonical_formula',	
    #                 'Canonical_molecular_weight', 'Canonical_charge', 'Extra_formula',
    #                 'Extra_molecular_weight', 'Extra_charge', 'BpForms_errors']
    # file_location = Path('~/karr_lab/datanator/docs/rna_modification/modomics.rrna.tsv').expanduser()
    # manager.fill_trna_collection(file_location, column_names=column_names,
    #                              start_row=None, query='kegg_orthology_name')

    manager.fill_dummy_orthodb()

if __name__ == '__main__':
    main()