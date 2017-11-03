from bioservices import UniProt
from kinetic_datanator.core import data_source
import sqlalchemy.ext.declarative
import pandas as pd
from sqlalchemy import Column, Integer, String, Float
from functools import lru_cache

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
    name = Column(String(255))
    gene_name = Column(String(255))
    canonical_sequence = Column(String(255))
    length = Column(Integer)
    mass = Column(String(255))



class Uniprot(data_source.HttpDataSource):
    """

    """
    base_model = Base


    def load_content(self):

        u = UniProt(verbose = False)

        dt = u.get_df(self.uniprot_list, nChunk = 200)

        columns = ['Entry', 'Entry name', 'Gene names  (primary )', 'Sequence', 'Length', 'Mass' ]

        pand = dt.loc[:, columns]
        new_columns = ['uniprot_id', 'name', 'gene_name', 'canonical_sequence', 'length', 'mass']
        pand.columns = new_columns

        pand.to_sql(name = 'uniprot', con=self.engine, if_exists = 'append', chunksize = 1000)
        self.session.commit()
