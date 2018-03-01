""" Downloads ontologies from BioPortal

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-23
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_source
import json
import kinetic_datanator.config.core
import os
import pickle
import pronto
import re
import requests
import sys


class BioPortal(data_source.DataSource):
    """ Loads ontologies from BioPortal

    Attributes:
        cache_dir (:obj:`str`): directory to store local copies of ontologies
        BIOPORTAL_ENDPOINT (:obj:`str`): URL pattern to download ontologies
        DEFAULT_CACHE_DIR (:obj:`str`): default directory to store local copies of ontologies
    """

    BIOPORTAL_ENDPOINT = 'http://data.bioontology.org'
    CCO_DOWNLOAD_URL = 'http://www.bio.ntnu.no/ontology/CCO/cco.obo'
    DEFAULT_CACHE_DIR = os.path.join(os.path.dirname(__file__), 'cache')

    def __init__(self, cache_dir=DEFAULT_CACHE_DIR):
        """
        Args:
            cache_dir (:obj:`str`, optional): directory to store local copies of ontologies
        """
        if not os.path.isdir(cache_dir):
            os.makedirs(cache_dir)
        self.cache_dir = cache_dir

    def get_ontologies(self):
        """ Get list of ontologies

        Returns:
            :obj:`list`: list of ontologies
        """
        filename = self.get_ontologies_filename()
        if not os.path.isfile(filename):
            self.download_ontologies()

        with open(filename, 'r') as file:
            return json.load(file)

    def get_ontologies_filename(self):
        """ Get the local filename to store a list of the ontologies

        Returns:
            :obj:`str`: filename
        """
        return os.path.join(self.cache_dir, 'bio_portal.ontologies.json')

    def download_ontologies(self):
        """
        Args:
            :obj:`list` of :obj:`str`: list of ontologies
        """
        response = requests.get(self.BIOPORTAL_ENDPOINT + '/ontologies',
                                headers={'Authorization': 'apikey token=' + self.get_api_key()})
        response.raise_for_status()

        ontologies_all = response.json()
        ontologies_lite = {}
        for ontology in ontologies_all:
            ontologies_lite[ontology['acronym']] = ontology['name']

        with open(self.get_ontologies_filename(), 'w') as file:
            json.dump(ontologies_lite, file)

    def get_ontology(self, id):
        """ Load ontology and download the ontology from BioPortal if neccessary

        Args:
            id (:obj:`str`): identifier of the ontology in BioPortal

        Returns:
            :obj:`pronto.Ontology`: ontology
        """
        filename = self.get_ontology_filename(id) + '.py' + str(sys.version_info[0]) + '.pkl'
        if not os.path.isfile(filename):
            self.download_ontology(id)

        with open(filename, 'rb') as file:
            return pickle.load(file)

    def get_ontology_filename(self, id):
        """ Get the local filename to store a copy of an ontology

        Args:
            id (:obj:`str`): identifier of the ontology in BioPortal

        Returns:
            :obj:`str`: filename
        """
        return os.path.join(self.cache_dir, id)

    def download_ontology(self, id):
        """ Download an ontology from BioPortal

        Args:
            id (:obj:`str`): identifier of the ontology in BioPortal
        """
        base_filename = self.get_ontology_filename(id)
        # download file
        if os.path.isfile(base_filename + '.obo'):
            filename_onto = base_filename + '.obo'
        elif os.path.isfile(base_filename + '.owl'):
            filename_onto = base_filename + '.owl'
        else:
            if id == 'CCO':
                response = requests.get(self.CCO_DOWNLOAD_URL)  # the BioPortal download link doesn't work
                response.raise_for_status()
                extension = 'obo'
            else:
                response = requests.get('{}/ontologies/{}/download'.format(self.BIOPORTAL_ENDPOINT, id),
                                        headers={'Authorization': 'apikey token=' + self.get_api_key()})
                response.raise_for_status()
                extension = re.findall('filename=".*?\.(.*?)"', response.headers['Content-Disposition'])[0]

            filename_onto = base_filename + '.' + extension
            with open(filename_onto, 'wb') as file:
                file.write(response.content)

        # parse and pickle result
        filename_pkl = base_filename + '.py' + str(sys.version_info[0]) + '.pkl'
        if not os.path.isfile(filename_pkl):
            ontology = pronto.Ontology(filename_onto)
            with open(filename_pkl, 'wb') as file:
                pickle.dump(ontology, file)

    def get_api_key(self):
        """ Get BioPortal API key

        Returns:
            :obj:`str`: key
        """
        return kinetic_datanator.config.core.get_config()['kinetic_datanator']['bioportal']['key']
