from kinetic_datanator.core import data_source
import sqlalchemy.ext.declarative
import pandas as pd
from sqlalchemy import Column, Integer, String, Float
from functools import lru_cache
from six.moves.urllib.request import urlretrieve
import gzip
from six import BytesIO


Base = sqlalchemy.ext.declarative.declarative_base()

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
    ENDPOINT_DOMAINS = {'uniprot' : 'http://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes#'}


    def load_content(self):

        #TODO: Figure out way to get the textfile from uniprot website

        pand = pd.read_csv('kinetic_datanator/data_source/cache/uniprot-all.txt', delimiter = '\t')

        new_columns = ['uniprot_id', 'entry_name', 'gene_name', 'protein_name', 'canonical_sequence', 'length', 'mass',
            'ec_number', 'entrez_id', 'status']
        pand.columns = new_columns

        pand['entrez_id'] = pand['entrez_id'].str.replace(';', '')

        if not self.max_entries == float('inf'):
            pand = pand[0:self.max_entries]

        pand.to_sql(name = 'uniprot', con=self.engine, if_exists = 'append', chunksize = 1000)
        self.session.commit()
