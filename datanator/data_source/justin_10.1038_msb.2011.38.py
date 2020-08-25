import zipfile
import pandas as pd 
from datanator_query_python.util import mongo_util
from datanator_query_python.config import config as q_conf
import os


class proteinHalfLives(mongo_util.MongoUtil):
    def __init__(self, MongoDB, db, collection, username, password,
                 path):
        super().__init__(MongoDB=MongoDB,
                         db=db,
                         username=username,
                         password=password)
        self.collection = self.db_obj[collection]
        self.MongoDB = MongoDB
        self.username = username
        self.password = password
        self.path = path


    def unzip_file(self):
        with zipfile.ZipFile(self.path) as zfile:
            zfile.extractall(self.path.replace('.zip', ''))

    
    def parse_protein(self):
        """parse data for protein half life and protein abundances
        """
        filepath = self.path.replace('.zip', '')
        data = pd.read_excel(filepath+"/Table_S4.xls")
        #print(data.columns) 

        for i in range(len(data)):
            d = {} 
            query = Query(MongoDB=self.MongoDB, db='datanator-test', entity_query_collection='uniprot', taxon_query_collection='taxon_tree', 
                          username=self.username, password=self.password)
            d['entity'] = {'type': 'protein', 
                           'name': query.query_entity(data.iloc[i,0][:6]),
                           'identifiers': [{'namespace': 'uniprot_id', 
                                            'value': data.iloc[i,0][:6]},
                                            {'namespace': 'entry_name',
                                            'value': data.iloc[i,0][7:]}]}
            d['identifier'] = {'namespace': 'uniprot_id', 'value': data.iloc[i,0][:6]}
            d['values'] = [{'type': 'half life', 'value': data['Half_Life/days'][i]*86400, 'units': 'seconds'}, 
                           {'type': 'protein abundance', 'value': data['protein abundance'][i], 'units': 'copies/cell'}]
            d['genotype'] = query.query_genotype('Mycoplasma pneumoniae')
            d['source'] = [{'namespace': 'doi', 'value': '10.1038/msb.2011.38'}]
            d['schema_version'] = '2.0'
            #print(d)
            self.collection.update_one({'type': 'protein', 
                                        'name': query.query_entity(data.iloc[i,0][:6]),
                                        'identifiers': [{'namespace': 'uniprot_id', 
                                                         'value': data.iloc[i,0][:6]},
                                                        {'namespace': 'entry_name',
                                                         'value': data.iloc[i,0][7:]}]},
                                        {'$set': d},
                                        upsert=True) 
            print('row {} has been added'.format(str(i)))


    def parse_rna(self):
        ''' parse data of RNA abundance 
        '''
        filepath = self.path.replace('.zip', '')
        data = pd.read_excel(filepath+'/Table_S4.xls')
        rows_no_gene_name = [] # list with rows which are missing a gene name in database

        for i in range(len(data)): # ROW 3 DOES NOT HAVE A GENE NAME
            d = {} 
            query = Query(MongoDB=self.MongoDB, db='datanator-test', entity_query_collection='uniprot', taxon_query_collection='taxon_tree', 
                          username=self.username, password=self.password)
            d['entity'] = {'type': 'RNA', 
                           'name': query.query_gene_name(data.iloc[i,0][:6]),
                           'identifiers': [query.query_gene_identifier(data.iloc[i,0][:6]),
                                           {'namespace': 'uniprot_id', 
                                           'value': data.iloc[i,0][:6]},
                                           {'namespace': 'entry_name',
                                           'value': data.iloc[i,0][7:]}]}
            d['identifier'] = query.query_gene_identifier(data.iloc[i,0][:6])
            d['values'] = [{'type': 'RNA abundance', 'value': data['mRNA abundance'][i], 'units': 'copies/cell'}]
            d['genotype'] = query.query_genotype('Mycoplasma pneumoniae')
            d['source'] = [{'namespace': 'doi', 'value': '10.1038/msb.2011.38'}]
            d['schema_version'] = '2.0'
            #print(d)
            if d['entity']['name'] == None: 
                rows_no_gene_name.append(i)
            else:
                self.collection.update_one({'type': 'RNA', 
                                            'name': query.query_entity(data.iloc[i,0][:6]),
                                            'identifiers': [{'namespace': 'uniprot_id', 
                                                             'value': data.iloc[i,0][:6]},
                                                            {'namespace': 'entry_name',
                                                             'value': data.iloc[i,0][7:]}]},
                                            {'$set': d},
                                            upsert=True)

            print('row {} has been added'.format(str(i)))
        
        print('Rows missing gene names: ' + str(rows_no_gene_name))
        
        
        
class Query(mongo_util.MongoUtil):
    def __init__(self, MongoDB, db, entity_query_collection, taxon_query_collection, username, password):
        super().__init__(MongoDB=MongoDB, db=db, username=username, password=password)
        self.entity_collection = self.db_obj[entity_query_collection]
        self.taxon_collection = self.db_obj[taxon_query_collection]


    def query_entity(self, uniprot_id):
        """takes the uniprot id and returns sring with protein name

        Args: 
            uniprot_id {:obj:`str`}: uniprotID of the wanted protein name

        Returns:
            protein_name (:obj:`str`}: string of the protein name
        """
        query = {'uniprot_id': uniprot_id}
        projection = {'protein_name': 1, '_id': 0}

        result = self.entity_collection.find_one(query, projection)
        return str(result['protein_name'])


    def query_genotype(self, tax_name):
        '''takes name of species and returns genotype object

        Args:
            tax_name {:obj:'str'}: name of species

        Returns:
            genotype_obj {:obj:obj}: genotype object for document
        '''
        query = {'tax_name': tax_name}
        projection = {'_id': 0, 'tax_id': 1, 'canon_anc_ids': 1, 'canon_anc_names': 1}
        result = self.taxon_collection.find_one(query, projection)
        genotype = {'taxon': {'ncbi_taxonomy_id': result['tax_id'], 'name': 'Mycoplasma pneumoniae', 'canon_ancestors': []}}

        for i in range(len(result['canon_anc_ids'])):
            genotype['taxon']['canon_ancestors'].append({'ncbi_taxonomy_id': result['canon_anc_ids'][i], 'name': result['canon_anc_names'][i]})

        return genotype

    
    def query_gene_identifier(self, uniprot_id):
        """ takes uniprot ID and returns the identifier object of the protein
        
        Args:
            uniprot_id {:obj:'str'}: uniprot ID of protein
        
        Returns:
            identifier {:obj:'obj'}: object for identifier
        """
        query = {'uniprot_id': uniprot_id}
        projection = {'_id': 0, 'ko_number': 1, 'gene_name': 1}
        result = self.entity_collection.find_one(query, projection)
        return {'namespace': 'ko_number', 'value': result['ko_number']}

    
    def query_gene_name(self, uniprot_id):
        """ takes uniprot ID and returns the gene name of the protein 
        Args:
            uniprot_id (:obj:'str'): uniprot ID of protein
        
        Return:
            gene_name (:obj:'str'): gene name of the protein
        """
        query = {'uniprot_id': uniprot_id}
        projection = {'_id': 0, 'gene_name': 1}
        result = self.entity_collection.find_one(query, projection)
        return result['gene_name']
        



def main():
    conf = q_conf.Justin()
    username = conf.USERNAME
    password = conf.PASSWORD
    MongoDB = conf.SERVER
    db = 'datanator-demo'
    collection = 'observation'
    filepath = 'datanator/docs/msb201138-sup-0003.zip' 
    src = proteinHalfLives(MongoDB=MongoDB, db=db, collection=collection, username=username, password=password, path=filepath)
    #src.unzip_file()
    #src.parse_protein()
    src.parse_rna()



if __name__ == "__main__":
    main()