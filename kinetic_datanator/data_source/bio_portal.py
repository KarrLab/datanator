"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-23
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_source
import kinetic_datanator.config.paths
import os
import pronto
import requests
import wc_utils.config.core

config_manager = wc_utils.config.core.ConfigManager(kinetic_datanator.config.paths.core)


class BioPortal(data_source.DataSource):
    """ Loads ontologies from BioPortal

    Attributes:
        cache_dir (:obj:`str`): directory to store local copies of ontologies
        BIOPORTAL_DOWNLOAD_ENDPOINT (:obj:`str`): URL pattern to download ontologies
        DEFAULT_CACHE_DIR (:obj:`str`): default directory to store local copies of ontologies
    """

    BIOPORTAL_DOWNLOAD_ENDPOINT = 'http://data.bioontology.org/ontologies/{}/download'
    DEFAULT_CACHE_DIR = os.path.join(os.path.dirname(__file__), 'cache')

    def __init__(self, cache_dir=DEFAULT_CACHE_DIR):
        """
        Args:
            cache_dir (:obj:`str`, optional): directory to store local copies of ontologies
        """
        if not os.path.isdir(cache_dir):
            os.makedirs(cache_dir)
        self.cache_dir = cache_dir

    def load_ontology(self, id):
        """ Load ontology and download the ontology from BioPortal if neccessary

        Args:
            id (:obj:`str`): identifier of the ontology in BioPortal

        Returns:
            :obj:`pronto.Ontology`: ontology
        """
        filename = self.get_ontology_filename(id)
        if not os.path.isfile(filename):
            self.download_ontology(id)
        ontology = pronto.Ontology(filename)
        return ontology

    def get_ontology_filename(self, id):
        """ Get the local filename to store a copy of an ontology

        Args:
            id (:obj:`str`): identifier of the ontology in BioPortal
        """
        return os.path.join(self.cache_dir, id + '.obo')

    def download_ontology(self, id):
        """ Download an ontology from BioPortal

        Args:
            id (:obj:`str`): identifier of the ontology in BioPortal
        """
        response = requests.get(self.BIOPORTAL_DOWNLOAD_ENDPOINT.format(id),
                                headers={'Authorization': 'apikey token=' + self.get_api_key()})
        response.raise_for_status()
        with open(self.get_ontology_filename(id), 'wb') as file:
            file.write(response.content)

    def get_api_key(self):
        """ Get BioPortal API key

        Returns:
            :obj:`str`: key
        """
        return config_manager.get_config()['bioportal']['key']
