from sqlalchemy import create_engine, ForeignKey, exists, Column, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref, sessionmaker
from kinetic_datanator.core import data_source


Base = declarative_base()



class GO(data_source.HttpDataSource):
    """ A local sqlite copy of the GO Ontology

    """

    base_model = Base
    ENDPOINT_DOMAINS = {
        'GO':
    }

    def load_content(self):
        pass
