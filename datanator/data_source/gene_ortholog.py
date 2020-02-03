from datanator_query_python.query import query_kegg_organism_code, query_uniprot, query_kegg_orthology, query_protein
from datanator_query_python.util import mongo_util
from datanator.data_source import uniprot_nosql
import datanator.config.core
import json
import requests
from bs4 import BeautifulSoup
import re


class KeggGeneOrtholog(mongo_util.MongoUtil):

    def __init__(self, server, src_db='datanator', des_db='datanator', collection_str='uniprot',
                username=None, password=None, readPreference='nearest', authSource='admin', verbose=True,
                max_entries=float('inf')):
        super().__init__(MongoDB=server, db=des_db, verbose=verbose, max_entries=max_entries,
                        username=username, password=password, authSource=authSource,
                        readPreference=readPreference)
        self.collection_str = collection_str
        self.max_entries = max_entries
        self.verbose = verbose
        self.des_client, self.des_db, self.des_collection = self.con_db(collection_str)
        self.koc_manager = query_kegg_organism_code.QueryKOC(username=username, password=password,
        server=server, authSource=authSource, collection_str='kegg_organism_code', readPreference=readPreference,
        database=src_db)
        self.uniprot_manager = query_uniprot.QueryUniprot(username=username, password=password, server=server,
        authSource=authSource, database=src_db, collection_str='uniprot', readPreference=readPreference)
        self.kegg_manager = query_kegg_orthology.QueryKO(username=username, password=password, server=server,
        authSource=authSource, database=src_db, max_entries=max_entries, verbose=verbose, readPreference=readPreference)
        self.protein_manager = query_protein.QueryProtein(username=username, password=password, server=server,
        authSource=authSource, database=src_db, max_entries=max_entries, verbose=verbose, readPreference=readPreference)
        self.uniprot_nosql_manager = uniprot_nosql.UniprotNoSQL(MongoDB=server, db=des_db, max_entries=max_entries,
        verbose=verbose, username=username, password=password, authSource=authSource)
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
    
    def parse_html(self, soup):
        """Parse out gene_orthologs from HTML 
        (https://www.kegg.jp/ssdb-bin/ssdb_best?org_gene=aly:ARALYDRAFT_486312).
        
        Args:
            soup (:obj:`BeautifulSoup`): BeautifulSoup object
        """
        org_gene_objs = soup.find_all(attrs={"type": "checkbox"})
        value_objs = soup.find_all(string=re.compile('<->'))
        for org_gene, value_str in zip(org_gene_objs, value_objs):
            values_list = ' '.join(value_str.split()).split()
            bits = int(values_list[-4])
            identity = float(values_list[-3])
            overlap = int(values_list[-2])            
            margin = values_list[-5]
            if margin[0] == '(':
                margin = int(margin[1:-1])
                length = int(values_list[-7])
                sw_score = int(values_list[-6])
            elif margin[0] == '-':
                margin = None
                length = int(values_list[-8])
                sw_score = int(values_list[-7])
            else:
                margin = int(margin[0:-1])
                length = int(values_list[-8])
                sw_score = int(values_list[-7])                                
            org_gene_str = org_gene.get('value')
            docs, count = self.uniprot_manager.get_id_by_org_gene(org_gene_str)
            ncbi_id = self.koc_manager.get_ncbi_by_org_code(org_gene_str.split(':')[0])
            uniprot_id = []
            if count != 0:
                for doc in docs:
                    uniprot_id.append(doc['uniprot_id'])
            else:
                proteins = self.parse_gene_info(org_gene_str.split(':')[1])
                if isinstance(proteins, str):
                    self.uniprot_nosql_manager.load_uniprot(query=True, msg=proteins.split('.')[0], species=[ncbi_id])
                    doc = self.uniprot_manager.get_info_by_entrez_id(org_gene_str.split(':')[1])
                    if doc is not None:
                        uniprot_id.append(doc)                    
                elif proteins != []: 
                    for protein in proteins:
                        self.uniprot_nosql_manager.load_uniprot(query=True, msg=protein.split('.')[0], species=[ncbi_id])
                        doc = self.uniprot_manager.get_info_by_entrez_id(org_gene_str.split(':')[1])
                        if doc is not None:
                            uniprot_id.append(doc)
            yield {'org_gene': org_gene_str,
                   'uniprot_id': uniprot_id,
                   'length': length,
                   'sw_score': sw_score,
                   'margin': margin,
                   'bits': bits,
                   'identity': identity,
                   'overlap': overlap,
                   'reference': self.endpoint + org_gene_str}

    def uniprot_to_org_gene(self, uniprot_id):
        """Given uniprot_id, convert to kegg org_gene format.
        
        Args:
            uniprot_id (:obj:`str`): Uniprot ID.

        Return:
            (:obj:`str`): Kegg org_gene format.
        """
        protein_doc = self.protein_manager.get_meta_by_id([uniprot_id])[0] #id from db so won't have missing entries
        ko = protein_doc.get('ko_number')
        if ko is None:
            return 'No KO info found.'
        ncbi_id = protein_doc['ncbi_taxonomy_id']
        org_code = self.koc_manager.get_org_code_by_ncbi(ncbi_id)
        gene_id = self.kegg_manager.get_gene_ortholog_by_id_org(ko, org_code)
        return '{}:{}'.format(org_code, gene_id)

    def parse_gene_info(self, gene):
        """Use mygene.info to get protein information
        given a string of gene code.
        
        Args:
            gene (:obj:`str`): Gene information.

        Return:
            (:obj:`list` of :obj:`str`): List of protein IDs. 
        """
        endpoint = 'https://mygene.info/v3/gene/' + gene
        r = requests.get(endpoint)
        r.raise_for_status
        info = json.loads(r.content)
        accession = info.get('accession')
        if accession is None:
            return []
        else:
            return accession.get('protein', [])

    def load_data(self, skip=0, top_hits=10):
        """Loading data.
        
        Args:
            skip (:obj:`int`, optional): Beginning of the documents. Defaults to 0.
            top_hits (:obj:`int`, optional): Number of top hits to iterate through. Defaults to 10.
        """
        con_0 = {'entrez_id': {'$ne': None}}
        con_1 = {'ko_number': {'$ne': "nan"}}
        query = {'$and': [con_0, con_1]}
        projection = {'_id': 0}
        docs = self.uniprot_manager.collection.find(filter=query, projection=projection, 
                                                    collation=self.uniprot_manager.collation, skip=skip,
                                                    no_cursor_timeout=True)
        count = self.uniprot_manager.collection.count_documents(query, collation=self.uniprot_manager.collation)
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            if i % 50 == 0 and self.verbose:
                print('Processing document {} out of {} ...'.format(i+skip, count))
            uniprot_id = doc['uniprot_id']
            org_gene_code = self.uniprot_to_org_gene(uniprot_id)
            print('    Processing protein {} whose org_gene code is: {}'.format(uniprot_id, org_gene_code))
            soup = self.get_html(org_gene_code)
            results = self.parse_html(soup)
            tmp = []
            for j, result in enumerate(results):
                if j == top_hits:
                    break
                else:
                    tmp.append(result)
            self.uniprot_manager.collection.update_one({'uniprot_id': uniprot_id},
                                                       {'$set': {'orthologs': tmp}}, upsert=False, 
                                                        collation=self.uniprot_manager.collation)
        docs.close()


def main():
    des_db = 'datanator'
    collection_str = 'uniprot'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    manager = KeggGeneOrtholog(server, collection_str=collection_str, username=username,
    password=password, des_db=des_db)
    manager.load_data(skip=158673)

if __name__ == '__main__':
    main()