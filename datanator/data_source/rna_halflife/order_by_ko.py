from datanator_query_python.query import query_uniprot_org
from datanator_query_python.util import mongo_util
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
        docs = self.src_collection.find(filter=query, skip=start)
        count = self.src_collection.count_documents(query)
        return docs, count

    def fill_helper(self, doi, field_name, start=0):
        """Method to fill new collection across different dois.
        
        Args:
            doi (:obj:`str`): DOI of publications.
            field_name (:obj:`str`): Name of the field that indicates the mRNA identifier.
            start (:obj:`int`, optional): Starting document position. Defaults to 0.
        """
        docs, count = self.helper(doi, start=start)
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            if self.verbose and i % 50 == 0:
                print('Processing doc {} out of {} ...'.format(i, count-start))
            for subdoc in doc:
                systematic_name = subdoc.get(field_name)
                species = subdoc.get('species')
                if systematic_name is not None:
                    uniprot_org_manager = query_uniprot_org.QueryUniprotOrg(systematic_name+' '+species)
                    uniprot_id = uniprot_org_manager.get_uniprot_id()
                    if uniprot_id is not None:
                        self.des_collection.update_one({'uniprot_id': uniprot_id},
                                                       {'$addToSet': {'halflives': subdoc}}, upsert=True, collation=self.collation)
                    else:
                        continue
                else:
                    continue

    def fill_cell(self, start=0):
        """Processing 10.1016/j.cell.2013.12.026.
        
        Args:
            start (:obj:`int`, optional): Starting document position. Defaults to 0.
        """
        doi = '10.1016/j.cell.2013.12.026'        
        self.fill_helper(doi, 'systematic_name', start=start)