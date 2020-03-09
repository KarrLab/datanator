from datanator_query_python.query import query_rna_halflife, query_uniprot
from datanator_query_python.util import mongo_util


class FillGeneName(mongo_util.MongoUtil):

    def __init__(self, server=None, db='datanator', collection_str='rna_halflife', username=None, 
                password=None, authSource='admin', readPreference='nearest', verbose=False, max_entries=float('inf')):
        super().__init__(MongoDB=server, db=db, verbose=verbose, max_entries=max_entries,
                        username=username, password=password, authSource=authSource, readPreference=readPreference)
        self.client, self.db, self.collection = self.con_db(collection_str)
        self.rna_query = query_rna_halflife.QueryRNA(server=server, username=username, password=password, verbose=verbose,
                                                db=db, collection_str=collection_str, authDB=authSource, readPreference=readPreference)
        self.uniprot_query = query_uniprot.QueryUniprot(username=username, password=password, server=server, authSource=authSource,
                                                        database=db, collection_str=collection_str, readPreference=readPreference)

    def fill_with_oln(self):
        """Fill gene_name with 'ordered_locus_name' field.
        """
        con_0 = {'gene_name': None}
        con_1 = {'halflives.ordered_locus_name': {'$exists': True}}
        pass