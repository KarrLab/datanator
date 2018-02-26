from kinetic_datanator.core import data_source
import sqlalchemy.ext.declarative
import pandas as pd
from sqlalchemy import Column, Integer, String, Float
# from six.moves.urllib.request import urlretrieve
# import gzip
# from six import BytesIO


Base = sqlalchemy.ext.declarative.declarative_base()
test = True

class UniprotData(Base):
    """ Represents protein interactions in from the IntAct Database

    Attributes:
        Index (:obj:`int`): Index of the DB
        interactor_a (:obj:`str`): represents participant A
        interactor_b (:obj:`str`): represents participant B
        publications (:obj:`str`): resource
        interaction (:obj:`str`): interaction ID
        feature_a (:obj:`str`): binding site of participant A
        feature_b (:obj:`str`): binding site of participant B
        stoich_a (:obj:`str`): stoichiometry of participant A
        stoich_b (:obj:`str`): stoichiometry of participant B

    """

    __tablename__ = 'uniprot'
    index = Column(Integer, primary_key = True)
    uniprot_id = Column(String(255), unique = True)
    entry_name = Column(String(255))
    gene_name = Column(String(255))
    protein_name = Column(String(255))
    canonical_sequence = Column(String(255))
    length = Column(Integer)
    mass = Column(String(255))
    ec_number = Column(String(255))
    entrez_id = Column(Integer)
    status  = Column(String(255))



class Uniprot(data_source.HttpDataSource):
    """

    """
    base_model = Base
    ENDPOINT_DOMAINS = {'uniprot' : 'http://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes#',
                        'uniprot-test': 'http://www.uniprot.org/uniprot/?sort=&desc=&compress=no&query=&fil=organism:%22Rattus%20norvegicus%20(Rat)%20[10116]%22%20AND%20reviewed:yes&limit=10&force=no&preview=true&format=tab&columns=id,entry%20name,genes(PREFERRED),protein%20names,sequence,length,mass,ec,database(GeneID),reviewed'}


    def load_content(self):

        if test:
            self.write_data_to_txt()
            pand = pd.read_csv('kinetic_datanator/data_source/cache/uniprot-test.txt', delimiter = '\t')
        else:
            pand = pd.read_csv('kinetic_datanator/data_source/cache/uniprot-all.txt', delimiter = '\t')

        new_columns = ['uniprot_id', 'entry_name', 'gene_name', 'protein_name', 'canonical_sequence', 'length', 'mass',
            'ec_number', 'entrez_id', 'status']
        pand.columns = new_columns

        pand['entrez_id'] = pand['entrez_id'].str.replace(';', '')

        if not self.max_entries == float('inf'):
            pand = pand[0:self.max_entries]

        pand.to_sql(name = 'uniprot', con=self.engine, if_exists = 'append', chunksize = 1000)
        self.session.commit()


    def write_data_to_txt(self):

        database_url = self.ENDPOINT_DOMAINS['uniprot-test']
        req = self.requests_session
        response = req.get(database_url)
        response.raise_for_status()
        content = str(response.content.decode('utf-8')).split('\n')
        with open(self.cache_dirname+'/uniprot-test.txt', 'w') as text_file:
            for line in content:
                text_file.write(line+'\n')
