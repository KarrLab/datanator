""" CHEBI Ontology Module
:Author:  Saahith Pochiraju  <saahith116@gmail.com>
:Date: 2018-04-09
:Copyright: 2018, Karr Lab
:License: MIT
"""

from sqlalchemy import create_engine, ForeignKey, exists, Column, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref, sessionmaker
from kinetic_datanator.core import data_source
import networkx
import obonet
import os


Base = declarative_base()

#TODO: Figure out if this is a support package for the website or will be integrated in to a CS tables
#TODO: Figure out whether pronto is better
#TODO: Include wc_util caching for caching graph

class Chebi(data_source.HttpDataSource):
    """ A local copy of the Chebi Ontology

    """
    base_model = Base
    ENDPOINT_DOMAINS = {'chebi_obo': "ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo"}

    def load_content(self):

        #NOTE: Estimated time for loading the graph is 150 seconds
        if os.path.isfile(self.cache_dirname+'/chebi.pkl'):
            self.graph = self.get_graph()
        else:
            self.graph = self.generate_graph()
            self.cache_graph(self.graph)

    def get_name(self, graph, id):
        return graph.nodes(data=True)[id]['name']

    def generate_graph(self):
        return obonet.read_obo(self.ENDPOINT_DOMAINS['chebi_obo'])

    def cache_graph(self, graph):
        networkx.write_gpickle(graph, self.cache_dirname+'/chebi.pkl')

    def get_graph(self):
        return networkx.read_gpickle(self.cache_dirname+'/chebi.pkl')
