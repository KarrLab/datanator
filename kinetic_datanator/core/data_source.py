"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-08
:Copyright: 2017, Karr Lab
:License: MIT
"""

import abc
import kinetic_datanator.config
import os
import requests
import requests_cache
import shutil
import six
import sqlalchemy
import sqlalchemy.orm
import sys
import tarfile
import tempfile
import wc_utils.quilt

CACHE_DIRNAME = os.path.join(os.path.dirname(__file__), '..', 'data_source', 'cache')
# :obj:`str`: default path for the sqlite database


class DataSource(six.with_metaclass(abc.ABCMeta, object)):
    """ Represents an external data source

    Attributes:
        name (:obj:`str`): name
    """

    def __init__(self, name=None):
        """
        Args:
            name (:obj:`str`, optional): name
        """
        if not name:
            name = self.__class__.__name__
        self.name = name


class CachedDataSource(DataSource):
    """ Represents an external data source that is cached locally in a sqlite database

    Attributes:
        filename (:obj:`str`): path to sqlite copy of the data source
        cache_dirname (:obj:`str`): directory to store the local copy of the data source
        engine (:obj:`sqlalchemy.engine.Engine`): sqlalchemy engine
        session (:obj:`sqlalchemy.orm.session.Session`): sqlalchemy session
        max_entries (:obj:`float`): maximum number of entries to save locally
        commit_intermediate_results (:obj:`bool`): if :obj:`True`, commit the changes throughout the loading
            process. This is particularly helpful for restarting this method when webservices go offline.
        verbose (:obj:`bool`): if :obj:`True`, print status information to the standard output
        quilt_owner (:obj:`str`): owner of Quilt package to save data
        quilt_package (:obj:`str`): identifier of Quilt package to save data

        base_model (:obj:`Base`): base ORM model for the sqlite databse
    """

    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backups=True, verbose=False, flask=False,
                 quilt_owner=None, quilt_package=None):
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
        """

        super(CachedDataSource, self).__init__(name=name)

        """ Set settings """
        # name
        if not cache_dirname:
            cache_dirname = CACHE_DIRNAME
        self.cache_dirname = cache_dirname
        self.filename = os.path.join(cache_dirname, self.name + '.sqlite')
        # loading
        self.max_entries = max_entries

        # committing
        self.commit_intermediate_results = commit_intermediate_results

        # verbosity
        self.verbose = verbose

        # flaskosity
        self.flask = flask

        # set Quilt configuration
        quilt_config = kinetic_datanator.config.get_config()['kinetic_datanator']['quilt']
        self.quilt_owner = quilt_owner or quilt_config['owner']
        self.quilt_package = quilt_package or quilt_config['package']

        """ Create SQLAlchemy session and load content if necessary """
        if os.path.isfile(self.filename):
            self.engine = self.get_engine()
            if clear_content:
                self.clear_content()
            self.session = self.get_session()
            if load_content:
                self.load_content()
        elif download_backups:
            self.download_backups()
            self.engine = self.get_engine()
            self.session = self.get_session()
            if load_content:
                self.load_content()
        else:
            self.engine = self.get_engine()
            if clear_content:
                self.clear_content()
            self.session = self.get_session()
            if load_content:
                self.load_content()

    def get_engine(self):
        """ Get an engine for the sqlite database. If the database doesn't exist, initialize its structure.

        Returns:
            :obj:`sqlalchemy.engine.Engine`: database engine
        """
        if not os.path.isdir(os.path.dirname(self.filename)):
            os.makedirs(os.path.dirname(self.filename))

        if self.flask:
            self.app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///' + self.filename
            engine = self.base_model.engine
            if not os.path.isfile(self.filename):
                self.base_model.metadata.create_all(engine)
        else:
            engine = sqlalchemy.create_engine('sqlite:///' + self.filename)
            if not os.path.isfile(self.filename):
                self.base_model.metadata.create_all(engine)

        return engine

    def clear_content(self):
        """ Clear the content of the sqlite database (i.e. drop and recreate all tables). """
        if self.flask:
            self.base_model.drop_all()
            self.base_model.create_all()
        else:
            self.base_model.metadata.drop_all(self.engine)
            self.base_model.metadata.create_all(self.engine)

    def get_session(self):
        """ Get a session for the sqlite database

        Returns:
            :obj:`sqlalchemy.orm.session.Session`: database session
        """
        if self.flask:
            return self.base_model.session
        else:
            return sqlalchemy.orm.sessionmaker(bind=self.engine)()

    def upload_backups(self):
        """ Backup the local sqlite database to Quilt """

        # create temporary directory to checkout package
        tmp_dirname = tempfile.mkdtemp()

        # install and export package
        manager = wc_utils.quilt.QuiltManager(tmp_dirname, self.quilt_package, owner=self.quilt_owner)
        manager.download()

        # copy new files to package
        paths = self.get_paths_to_backup()
        for path in paths:
            if os.path.isfile(os.path.join(self.cache_dirname, path)):
                shutil.copyfile(os.path.join(self.cache_dirname, path), os.path.join(tmp_dirname, path))
            else:
                shutil.copytree(os.path.join(self.cache_dirname, path), os.path.join(tmp_dirname, path))

        # build and push package
        manager.upload()

        # cleanup temporary directory
        shutil.rmtree(tmp_dirname)

    def download_backups(self):
        """ Download the local sqlite database from Quilt """

        # create temporary directory to checkout package
        tmp_dirname = tempfile.mkdtemp()

        # install and export package
        manager = wc_utils.quilt.QuiltManager(tmp_dirname, self.quilt_package, owner=self.quilt_owner)
        manager.download()

        # copy requested files from package
        paths = self.get_paths_to_backup(download=True)
        for path in paths:
            if os.path.isfile(os.path.join(tmp_dirname, path)):
                shutil.copyfile(os.path.join(tmp_dirname, path), os.path.join(self.cache_dirname, path))
            else:
                shutil.copytree(os.path.join(tmp_dirname, path), os.path.join(self.cache_dirname, path))

        # cleanup temporary directory
        shutil.rmtree(tmp_dirname)

    def get_paths_to_backup(self, download=False):
        """ Get a list of the files to backup/unpack

        Args:
            download (:obj:`bool`, optional): if :obj:`True`, prepare the files for uploading

        Returns:
            :obj:`list` of :obj:`str`: list of paths to backup
        """
        paths = []
        paths.append(self.name + '.sqlite')
        return paths

    @abc.abstractmethod
    def load_content(self):
        """ Load the content of the local copy of the data source """
        pass

    def get_or_create_object(self, cls, **kwargs):
        """ Get the SQLAlchemy object of type :obj:`cls` with attribute/value pairs specified by `**kwargs`. If
        an object with these attribute/value pairs does not exist, create an object with these attribute/value pairs
        and add it to the SQLAlchemy session.

        Args:
            cls (:obj:`class`): child class of :obj:`base_model`
            **kwargs (:obj:`dict`, optional): attribute-value pairs of desired SQLAlchemy object of type :obj:`cls`

        Returns:
            :obj:`base_model`: SQLAlchemy object of type :obj:`cls`
        """
        q = self.session.query(cls).filter_by(**kwargs)
        self.session.flush()
        if q.count():
            return q.first()
        else:
            obj = cls(**kwargs)
            self.session.add(obj)
            return obj


class FtpDataSource(CachedDataSource):
    """ An external data source which can be obtained via a FTP interface

    Attributes:
        ENDPOINT_DOMAINS (:obj:`dict` of :obj:`str`, :obj:`str`): dictionary of domains to retry
    """
    ENDPOINT_DOMAINS = {}


class HttpDataSource(CachedDataSource):
    """ An external data source which can be obtained via a HTTP interface

    Attributes:
        requests_cache_filename (:obj:`str`): path to cache HTTP requests
        requests_session (:obj:`requests_cache.core.CachedSession`): cache-enabled HTTP request session
        ENDPOINT_DOMAINS (:obj:`dict` of :obj:`str`, :obj:`str`): dictionary of domains to retry
        MAX_HTTP_RETRIES (:obj:`int`): maximum number of times to retry each HTTP request
    """

    ENDPOINT_DOMAINS = {}
    MAX_HTTP_RETRIES = 5

    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backups=True, verbose=False,
                 clear_requests_cache=False, download_request_backup=False,
                 quilt_owner=None, quilt_package=None,
                 flask=False):
        """
        Args:
            name (:obj:`str`, optional): name
            cache_dirname (:obj:`str`, optional): directory to store the local copy of the data source and the HTTP requests cache
            clear_content (:obj:`bool`, optional): if :obj:`True`, clear the content of the sqlite local copy of the data source
            load_content (:obj:`bool`, optional): if :obj:`True`, load the content of the local sqlite database from the external source
            max_entries (:obj:`float`, optional): maximum number of entries to save locally
            commit_intermediate_results (:obj:`bool`, optional): if :obj:`True`, commit the changes throughout the loading
                process. This is particularly helpful for restarting this method when webservices go offline.
            download_backups (:obj:`bool`, optional): if :obj:`True`, load the local copy of the data source from the Karr Lab server
            verbose (:obj:`bool`, optional): if :obj:`True`, print status information to the standard output
            clear_requests_cache (:obj:`bool`, optional): if :obj:`True`, clear the HTTP requests cache
            download_request_backup (:obj:`bool`, optional): if :obj:`True`, download the request backup
            quilt_owner (:obj:`str`, optional): owner of Quilt package to save data
            quilt_package (:obj:`str`, optional): identifier of Quilt package to save data
            flask (:obj:`bool`, optional)
        """

        """ CachedDataSource settings """
        if not name:
            name = self.__class__.__name__
        if not cache_dirname:
            cache_dirname = CACHE_DIRNAME

        """ Request settings """
        # todo (enhancement): avoid python version-specific requests cache; this currently is necessary because request_cache uses
        # pickle which is not backwards compatible
        self.requests_cache_filename = os.path.join(cache_dirname, name + '.requests.py{}.sqlite'.format(sys.version_info[0]))
        self.requests_session = self.get_requests_session()

        if clear_requests_cache:
            self.clear_requests_cache()

        self.download_request_backup = download_request_backup

        """ Call superclass constructor which will optionally load content """
        super(HttpDataSource, self).__init__(name=name, cache_dirname=cache_dirname,
                                             clear_content=clear_content, load_content=load_content, max_entries=max_entries,
                                             commit_intermediate_results=commit_intermediate_results,
                                             download_backups=download_backups, verbose=verbose,
                                             quilt_owner=quilt_owner, quilt_package=quilt_package,
                                             flask=flask)

    def get_requests_session(self):
        """ Setup an cache-enabled HTTP request session

        Returns:
            :obj:`requests_cache.core.CachedSession`: cached-enable session
        """
        if not os.path.isdir(os.path.dirname(self.requests_cache_filename)):
            os.makedirs(os.path.dirname(self.requests_cache_filename))
        name, _, _ = self.requests_cache_filename.rpartition('.')

        # create caching session
        session = requests_cache.core.CachedSession(name, backend='sqlite', expire_after=None)

        # setup retrying
        for endpoint_domain in self.ENDPOINT_DOMAINS.values():
            session.mount(endpoint_domain, requests.adapters.HTTPAdapter(max_retries=self.MAX_HTTP_RETRIES))

        return session

    def clear_requests_cache(self):
        """ Clear the cache-enabled HTTP request session """
        self.requests_session.cache.clear()

    def get_paths_to_backup(self, download=False):
        """ Get a list of the files to backup/unpack

        Args:
            download (:obj:`bool`, optional): if :obj:`True`, prepare the files for uploading

        Returns:
            :obj:`list` of :obj:`str`: paths to backup
        """
        paths = super(HttpDataSource, self).get_paths_to_backup(download=download)
        if not download or self.download_request_backup:
            paths.append(self.name + '.requests.py{}.sqlite'.format(sys.version_info[0]))
        return paths


class WebserviceDataSource(DataSource):
    """ A data source that is a webservice

    Attributes:
        requests_session (:obj:`requests.Session`): cache-enabled HTTP request session
        ENDPOINT_DOMAINS (:obj:`dict` of :obj:`str`, :obj:`str`): dictionary of domains to retry
        MAX_HTTP_RETRIES (:obj:`int`): maximum number of times to retry each HTTP request
    """

    ENDPOINT_DOMAINS = {}
    MAX_HTTP_RETRIES = 5

    def __init__(self):
        self.requests_session = requests.Session()
        for endpoint_domain in self.ENDPOINT_DOMAINS.values():
            self.requests_session.mount(endpoint_domain, requests.adapters.HTTPAdapter(max_retries=self.MAX_HTTP_RETRIES))


class DataSourceWarning(UserWarning):
    """ Data source warning """
    pass
