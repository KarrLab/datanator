from datanator_query_python.query import query_uniprot_org
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config
from pymongo.collation import Collation, CollationStrength


class Reorg:
    """Reorganize docs into categories uniprot_id.
    
    Args:
        mongo_util (): [description]
    """

    def __init__(self, cache_dirname=None, MongoDB=None, src_db='datanator',
                 verbose=False, max_entries=float('inf'), username=None, 
                 password = None, authSource='admin', readPreference='nearest',
                 des_collection='rna_halflife_new', src_collection='rna_halflife',
                 des_db='test'):
        """Init.
        
        Args:
            cache_dirname ([type], optional): [description]. Defaults to None.
            MongoDB ([type], optional): [description]. Defaults to None.
            replicaSet ([type], optional): [description]. Defaults to None.
            db (str, optional): [description]. Defaults to 'test'.
            verbose (bool, optional): [description]. Defaults to False.
            max_entries ([type], optional): [description]. Defaults to float('inf').
            username ([type], optional): [description]. Defaults to None.
            password ([type], optional): [description]. Defaults to None.
            authSource (str, optional): [description]. Defaults to 'admin'.
            readPreference (str, optional): [description]. Defaults to 'nearest'.
        """
        self.max_entries = max_entries
        self.verbose = verbose
        self.src_client, self.src_db, self.src_collection = mongo_util.MongoUtil(cache_dirname=cache_dirname, MongoDB=MongoDB, db=src_db,
                                                                                verbose=verbose, max_entries=max_entries, username=username, 
                                                                                password=password, authSource=authSource, readPreference=readPreference).con_db(collection_str=src_collection)
        self.des_client, self.des_db, self.des_collection = mongo_util.MongoUtil(cache_dirname=cache_dirname, MongoDB=MongoDB, db=des_db,
                                                                                verbose=verbose, max_entries=max_entries, username=username, 
                                                                                password=password, authSource=authSource, readPreference=readPreference).con_db(collection_str=des_collection)
        self.collation = Collation('en', strength=CollationStrength.SECONDARY)

    def helper(self, doi, start=0):
        """helper function for each publication
        
        Args:
            doi (:obj:`str`): DOI of publication.
        """
        query = {'halflives.reference.doi': doi}
        # project = {
        #     '$project': {
        #         'halflives': {
        #             '$filter': {
        #                 'input': "$halflives",
        #                 'as': "ref",
        #                 'cond': {'$eq': ['$$ref.reference', {'doi': doi}]}
        #             }
        #         }
        #     }
        # }
        # pipeline = [{'$match': query}, project]
        docs = self.src_collection.find(filter=query, collation=self.collation)
        count = self.src_collection.count_documents(query)
        return docs, count

    def fill_helper(self, doi, field_name, start=0, species=None):
        """Method to fill new collection across different dois.
        
        Args:
            doi (:obj:`str`): DOI of publications.
            field_name (:obj:`str`): Name of the field that indicates the mRNA identifier.
            start (:obj:`int`, optional): Starting document position. Defaults to 0.
            species (:obj:`str`, optional): NCBI Taxonomy name of the organism
        """
        docs, count = self.helper(doi, start=start)
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            if self.verbose and i % 50 == 0:
                print('Processing doc {} out of {} ...'.format(i, count-start))
            for subdoc in doc['halflives']:
                reference = subdoc.get('reference')[0]['doi']
                if reference != doi:
                    continue
                if doi == '10.1371/journal.pone.0059059':
                    systematic_name = doc['gene_name']
                else:
                    systematic_name = subdoc.get(field_name)
                    if isinstance(systematic_name, list):
                        systematic_name = systematic_name[0]
                if species is None:
                    species = subdoc.get('species')

                uniprot_org_manager = query_uniprot_org.QueryUniprotOrg(systematic_name+' '+species)
                uniprot_id = uniprot_org_manager.get_uniprot_id()
                ko = uniprot_org_manager.get_kegg_ortholog()
                protein_names = uniprot_org_manager.get_protein_name()
                if uniprot_id is not None:
                    self.des_collection.update_one({'uniprot_id': uniprot_id},
                                                    {'$addToSet': {'halflives': subdoc},
                                                     '$set': {'protein_names': protein_names,
                                                              'ko_number': ko}}, upsert=True, collation=self.collation)                
                else:
                    uniprot_org_manager = query_uniprot_org.QueryUniprotOrg(systematic_name)
                    uniprot_id = uniprot_org_manager.get_uniprot_id()
                    ko = uniprot_org_manager.get_kegg_ortholog()
                    protein_names = uniprot_org_manager.get_protein_name()
                    if uniprot_id is not None:
                        self.des_collection.update_one({'uniprot_id': uniprot_id},
                                                        {'$addToSet': {'halflives': subdoc},
                                                        '$set': {'protein_names': protein_names,
                                                                 'ko_number': ko}}, upsert=True, collation=self.collation)
                    else:
                        self.des_collection.update_one({'identifier': systematic_name},
                                                        {'$addToSet': {'halflives': subdoc},
                                                        '$set': {'protein_names': [protein_names],
                                                                 'ko_number': None}}, upsert=True, collation=self.collation)
                        print(systematic_name) 


    def fill_cell(self, start=0):
        """Processing 10.1016/j.cell.2013.12.026.
        
        Args:
            start (:obj:`int`, optional): Starting document position. Defaults to 0.
        """
        doi = '10.1016/j.cell.2013.12.026'        
        self.fill_helper(doi, 'systematic_name', start=start)

    def fill_mbc(self, start=0):
        """Processing 10.1091/mbc.e11-01-0028
        
        Args:
            start (:obj:`int`, optional): Starting document position. Defaults to 0.
        """
        doi = '10.1091/mbc.e11-01-0028'        
        self.fill_helper(doi, 'systematic_name', start=start, species='Saccharomyces cerevisiae')

    def fill_nar_gks(self, start=0):
        """Processing 10.1093/nar/gks1019
        
        Args:
            start (:obj:`int`, optional): Starting document position. Defaults to 0.
        """
        doi = '10.1093/nar/gks1019'        
        self.fill_helper(doi, 'ordered_locus_name', start=start, species='Mycolicibacterium smegmatis')

    def fill_nar_gkt(self, start=0):
        """Processing 10.1093/nar/gkt1150
        
        Args:
            start (:obj:`int`, optional): Starting document position. Defaults to 0.
        """
        doi = '10.1093/nar/gkt1150'        
        self.fill_helper(doi, 'ordered_locus_name', start=start, species='Escherichia coli strain K-12')

    def fill_gr_131(self, start=0):
        """Processing 10.1101/gr.131037.111
        
        Args:
            start (:obj:`int`, optional): Starting document position. Defaults to 0.
        """
        doi = '10.1101/gr.131037.111'        
        self.fill_helper(doi, 'accession_id', start=start, species='Mus musculus')

    def fill_gb_2012(self, start=0):
        """Processing 10.1186/gb-2012-13-4-r30
        
        Args:
            start (:obj:`int`, optional): Starting document position. Defaults to 0.
        """
        doi = '10.1186/gb-2012-13-4-r30'        
        self.fill_helper(doi, 'ordered_locus_name', start=start)

    def fill_s12864(self, start=0):
        """Processing 10.1186/s12864-016-3219-8
        
        Args:
            start (:obj:`int`, optional): Starting document position. Defaults to 0.
        """
        doi = '10.1186/s12864-016-3219-8'        
        self.fill_helper(doi, 'ordered_locus_name', start=start, species='Methanosarcina acetivorans')

    def fill_journal_pone(self, start=0):
        """Processing 10.1371/journal.pone.0059059
        
        Args:
            start (:obj:`int`, optional): Starting document position. Defaults to 0.
        """
        doi = '10.1371/journal.pone.0059059'        
        self.fill_helper(doi, 'ordered_locus_name', start=start, species='Lactococcus lactis subsp. lactis (strain IL1403)')

from multiprocessing import Process
import datanator.config.core

def joint_operation(src):
    src.fill_cell()
    src.fill_mbc()
    src.fill_nar_gks()
    src.fill_nar_gkt()
    src.fill_gr_131()
    src.fill_gb_2012()
    src.fill_s12864()
    src.fill_journal_pone()

def main():
    des_db = 'datanator'
    src_db = 'datanator'
    src_collection = 'rna_halflife'
    des_collection = 'rna_halflife_new'
    username = datanator.config.core.get_config()['datanator']['mongodb']['user']
    password = datanator.config.core.get_config()['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']    
    src = Reorg(MongoDB=MongoDB, src_db=src_db,
                verbose=False, username=username, 
                password=password, authSource='admin', readPreference='nearest',
                des_collection=des_collection, src_collection=src_collection,
                des_db=des_db)
    p = Process(target=joint_operation(src))
    p.start()
    p.join()

if __name__ == '__main__':
    main()