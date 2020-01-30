from datanator_query_python.query import query_kegg_organism_code
from datanator_query_python.util import mongo_util
import requests
from bs4 import BeautifulSoup


class KeggGeneOrtholog(mongo_util.MongoUtil):

    def __init__(self, server, src_db='datanator', des_db='datanator', collection_str='uniprot',
                username=None, password=None, readPreference='nearest', authSource='admin', verbose=True,
                max_entries=float('inf')):
        super().__init__(MongoDB=server, db=des_db, verbose=verbose, max_entries=max_entries,
                        username=username, password=password, authSource=authSource,
                        readPreference=readPreference)
        self.des_client, self.des_db, self.des_collection = self.con_db(collection_str)
        self.koc_manager = query_kegg_organism_code.QueryKOC(username=username, password=password,
        server=server, authSource=authSource, collection_str='kegg_organism_code', readPreference=readPreference)
        self.endpoint = 'https://www.kegg.jp/ssdb-bin/ssdb_best?org_gene='

    def get_html(self, query):
        """Get HTML file based on org:gene_code string,
        e.g. aly:ARALYDRAFT_486312.
        
        Args:
            query (:obj:`str`): org:gene_code string.
        """
        r = requests.get(self.endpoint + query)
        soup = BeautifulSoup(r.content, 'html.parser')
        return soup
    
    def parse_html(self, soup, top_hits=10):
        """Parse out the top_hits number of gene_orthologs
        from HTML e.g. (https://www.kegg.jp/ssdb-bin/ssdb_best?org_gene=aly:ARALYDRAFT_486312).
        
        Args:
            soup (:obj:`BeautifulSoup`): BeautifulSoup object
            top_hits (:obj:`int`, optional): Number of hits needed. Defaults to 10.
        """
        objs = soup.find_all(attrs={"title": "species"}, limit=top_hits)
        for obj in objs:
            yield obj.get('value')