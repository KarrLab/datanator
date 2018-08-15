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


class BioPortal(data_source.CachedDataSource):
    """ Loads ontologies from BioPortal

    Attributes:
        ontologies (:obj:`list`): list of filenames of ontologies

        BIOPORTAL_ENDPOINT (:obj:`str`): URL pattern to download ontologies
        CCO_DOWNLOAD_URL (:obj:`str`): URL to download CCO ontology
        DEFAULT_CACHE_DIR (:obj:`str`): default directory to store local copies of ontologies
    """

    BIOPORTAL_ENDPOINT = 'http://data.bioontology.org'
    CCO_DOWNLOAD_URL = 'http://www.bio.ntnu.no/ontology/CCO/cco.obo'
    DEFAULT_CACHE_DIR = os.path.join(os.path.dirname(__file__), 'cache')

    DEFAULT_ONTOLOGIES = (
        'BTO.obo',
        'CCO.obo',
        'CL.owl',
        'DOID.obo',
        'EFO.owl',
        'FMA.owl',
        'GO.obo',
        'PW.obo',
        'SBO.obo',
    )

    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backups=True, verbose=False, flask=False,
                 quilt_owner=None, quilt_package=None, ontologies=None):
        """
        Args:
            name (:obj:`str`, optional): name
            cache_dirname (:obj:`str`, optional): directory to store the local copy of the data source
            clear_content (:obj:`bool`, optional): if :obj:`True`, clear the content of the sqlite local copy of the data source
            load_content (:obj:`bool`, optional): if :obj:`True`, load the content of the local sqlite database from the external source
            max_entries (:obj:`float`, optional): maximum number of entries to save locally
            commit_intermediate_results (:obj:`bool`, optional): if :obj:`True`, commit the changes throughout the loading
                process. This is particularly helpful for restarting this method when webservices go offline.
            download_backups (:obj:`bool`, optional): if :obj:`True`, load the local copy of the data source from the Karr Lab server
            verbose (:obj:`bool`, optional): if :obj:`True`, print status information to the standard output
            quilt_owner (:obj:`str`, optional): owner of Quilt package to save data
            quilt_package (:obj:`str`, optional): identifier of Quilt package to save data
            ontologies (:obj:`list`, optional): list of filenames of ontologies
        """
        if ontologies is None:
            ontologies = self.DEFAULT_ONTOLOGIES
        self.ontologies = ontologies

        if not name:
            name = self.__class__.__name__
        self.name = name

        if not cache_dirname:
            cache_dirname = data_source.CACHE_DIRNAME
        self.cache_dirname = cache_dirname

        # verbosity
        self.verbose = verbose

        # set Quilt configuration
        quilt_config = kinetic_datanator.config.get_config()['kinetic_datanator']['quilt']
        self.quilt_owner = quilt_owner or quilt_config['owner']
        self.quilt_package = quilt_package or quilt_config['package']

        """ load content as necessary """
        has_all_files = True
        for path in self.get_paths_to_backup():
            if not os.path.isfile(os.path.join(self.cache_dirname, path)):
                has_all_files = False
                break

        if has_all_files:
            if clear_content:
                self.clear_content()
            if load_content:
                self.load_content()
        elif download_backups:
            self.download_backups()
            if load_content:
                self.load_content()
        else:
            if clear_content:
                self.clear_content()
            if load_content:
                self.load_content()

    def get_engine(self):
        """ Get an engine for the sqlite database. If the database doesn't exist, initialize its structure.

        Returns:
            :obj:`sqlalchemy.engine.Engine`: database engine
        """
        pass

    def get_session(self):
        """ Get a session for the sqlite database

        Returns:
            :obj:`sqlalchemy.orm.session.Session`: database session
        """
        pass

    def get_paths_to_backup(self, download=False):
        """ Get a list of the files to backup/unpack

        Args:
            download (:obj:`bool`, optional): if :obj:`True`, prepare the files for uploading

        Returns:
            :obj:`list` of :obj:`str`: list of paths to backup
        """
        paths = ['bio_portal.ontologies.json']
        for filename in self.ontologies:
            paths.append(filename)
            name, ext = os.path.splitext(filename)
            paths.append(name + '.py3.pkl')
        return paths

    def clear_content(self):
        """ Clear the content of the sqlite database (i.e. drop and recreate all tables). """
        for path in self.get_paths_to_backup():
            if os.path.isfile(os.path.join(self.cache_dirname, path)):
                os.remove(os.path.join(self.cache_dirname, path))

    def load_content(self):
        """ Load the content of the local copy of the data source """
        if self.verbose:
            print('Downloading ontologies ...')
        self.download_ontologies()

        for filename in self.ontologies:
            if self.verbose:
                print('  Dowloading {} ...'.format(filename))
            name, ext = os.path.splitext(filename)
            self.download_ontology(name)

    def get_ontologies(self):
        """ Get list of ontologies

        Returns:
            :obj:`list`: list of ontologies
        """
        self.download_ontologies()

        filename = self.get_ontologies_filename()
        with open(filename, 'r') as file:
            return json.load(file)

    def get_ontologies_filename(self):
        """ Get the local filename to store a list of the ontologies

        Returns:
            :obj:`str`: filename
        """
        return os.path.join(self.cache_dirname, 'bio_portal.ontologies.json')

    def download_ontologies(self):
        """
        Args:
            :obj:`list` of :obj:`str`: list of ontologies
        """
        filename = self.get_ontologies_filename()
        if os.path.isfile(filename):
            return

        response = requests.get(self.BIOPORTAL_ENDPOINT + '/ontologies',
                                headers={'Authorization': 'apikey token=' + self.get_api_key()})
        response.raise_for_status()

        ontologies_all = response.json()
        ontologies_acronyms = {}
        for ontology in ontologies_all:
            ontologies_acronyms[ontology['acronym']] = ontology['name']

        with open(filename, 'w') as file:
            json.dump(ontologies_acronyms, file)

    def get_ontology(self, id):
        """ Load ontology and download the ontology from BioPortal if neccessary

        Args:
            id (:obj:`str`): identifier of the ontology in BioPortal

        Returns:
            :obj:`pronto.Ontology`: ontology
        """
        self.download_ontology(id)

        filename = self.get_ontology_filename(id) + '.py' + str(sys.version_info[0]) + '.pkl'
        with open(filename, 'rb') as file:
            return pickle.load(file)

    def get_ontology_filename(self, id):
        """ Get the local filename to store a copy of an ontology

        Args:
            id (:obj:`str`): identifier of the ontology in BioPortal

        Returns:
            :obj:`str`: filename
        """
        return os.path.join(self.cache_dirname, id)

    def download_ontology(self, id):
        """ Download an ontology from BioPortal

        Args:
            id (:obj:`str`): identifier of the ontology in BioPortal
        """
        base_filename = self.get_ontology_filename(id)
        filename_pkl = base_filename + '.py' + str(sys.version_info[0]) + '.pkl'
        if os.path.isfile(filename_pkl):
            return

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
