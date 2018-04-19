from sqlalchemy import create_engine, ForeignKey, exists, Column, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref, sessionmaker
from kinetic_datanator.core import data_source
import os
import networkx
import obonet

Base = declarative_base()



class GO(data_source.HttpDataSource):
    """ A local sqlite copy of the GO Ontology

    """

    base_model = Base
    ENDPOINT_DOMAINS = {'GO_obo': 'http://purl.obolibrary.org/obo/go.obo'}


    def load_content(self):

        #NOTE: Estimated time for loading the graph is 150 seconds
        if os.path.isfile(self.cache_dirname+'/go.pkl'):
            self.graph = self.get_graph()
        else:
            self.graph = self.generate_graph()
            self.cache_graph(self.graph)


    def get_name(self, graph, id):
        return graph.nodes(data=True)[id]['name']

    def generate_graph(self):
        return obonet.read_obo(self.ENDPOINT_DOMAINS['GO_obo'])

    def cache_graph(self, graph):
        networkx.write_gpickle(graph, self.cache_dirname+'/go.pkl')

    def get_graph(self):
        return networkx.read_gpickle(self.cache_dirname+'/go.pkl')
