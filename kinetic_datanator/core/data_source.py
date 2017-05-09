"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-08
:Copyright: 2017, Karr Lab
:License: MIT
"""

from wc_utils import backup
import abc
import os
import requests
import requests_cache
import six
import sqlalchemy
import sqlalchemy.orm
import sys

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
        backup_server_token (:obj:`str`): authentification token to upload/download copies of the data source to/from the Karr Lab server
        engine (:obj:`sqlalchemy.engine.Engine`): sqlalchemy engine
        session (:obj:`sqlalchemy.orm.session.Session`): sqlalchemy session
        max_entries (:obj:`float`): maximum number of entries to save locally
        verbose (:obj:`bool`): if :obj:`True`, print status information to the standard output

        base_model (:obj:`Base`): base ORM model for the sqlite databse
    """

    def __init__(self, name=None, cache_dirname=CACHE_DIRNAME, clear_content=False, load_content=False, max_entries=float('inf'),
                 download_backup=True, verbose=False):
        """
        Args:
            name (:obj:`str`, optional): name
            cache_dirname (:obj:`str`, optional): directory to store the local copy of the data source
            clear_content (:obj:`bool`, optional): if :obj:`True`, clear the content of the sqlite local copy of the data source
            load_content (:obj:`bool`, optional): if :obj:`True`, load the content of the local sqlite database from the external source
            max_entries (:obj:`float`, optional): maximum number of entries to save locally
            download_backup (:obj:`bool`, optional): if :obj:`True`, load the local copy of the data source from the Karr Lab server
            verbose (:obj:`bool`, optional): if :obj:`True`, print status information to the standard output
        """

        super(CachedDataSource, self).__init__(name=name)

        """ Set settings """
        # name
        self.cache_dirname = cache_dirname
        self.filename = os.path.join(cache_dirname, self.name + '.sqlite')

        # loading
        self.max_entries = max_entries

        # backup settings
        self.backup_server_token = os.getenv('CODE_SERVER_TOKEN')

        # verbosity
        self.verbose = verbose

        """ Create SQLAlchemy session and load content if necessary """
        if os.path.isfile(self.filename):
            self.engine = self.get_engine()
            if clear_content:
                self.clear_content()
            self.session = self.get_session()
            if load_content:
                self.load_content()
        elif download_backup and self.backup_server_token:
            self.download_backup()
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

        engine = sqlalchemy.create_engine('sqlite:///' + self.filename)

        if not os.path.isfile(self.filename):
            self.base_model.metadata.create_all(engine)

        return engine

    def clear_content(self):
        """ Clear the content of the sqlite database (i.e. drop and recreate all tables). """
        self.base_model.metadata.drop_all(self.engine)
        self.base_model.metadata.create_all(self.engine)

    def get_session(self):
        """ Get a session for the sqlite database

        Returns:
            :obj:`sqlalchemy.orm.session.Session`: database session
        """
        return sqlalchemy.orm.sessionmaker(bind=self.engine)()

    def upload_backup(self):
        """ Backup the local sqlite database to the Karr Lab server """
        backup.BackupManager(self.filename, arcname=self.name + '.sqlite', token=self.backup_server_token) \
            .create() \
            .upload() \
            .cleanup()

    def download_backup(self):
        """ Download the local sqlite database from the Karr Lab server """
        backup.BackupManager(self.filename, arcname=self.name + '.sqlite', token=self.backup_server_token) \
            .download() \
            .extract() \
            .cleanup()

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
        if q.count():
            return q.first()
        else:
            obj = cls(**kwargs)
            self.session.add(obj)
            return obj


class HttpDataSource(CachedDataSource):
    """ An external data source which can be obtained via a HTTP interface 

    Attributes:
        requests_cache_filename (:obj:`str`): path to cache HTTP requests
        requests_session (:obj:`requests_cache.core.CachedSession`): cache-enabled HTTP request session
    """

    def __init__(self, name=None, cache_dirname=CACHE_DIRNAME, clear_content=False, load_content=False, max_entries=float('inf'),
                 download_backup=True, verbose=False,
                 clear_requests_cache=False):
        """
        Args:
            name (:obj:`str`, optional): name
            cache_dirname (:obj:`str`, optional): directory to store the local copy of the data source and the HTTP requests cache
            clear_content (:obj:`bool`, optional): if :obj:`True`, clear the content of the sqlite local copy of the data source
            load_content (:obj:`bool`, optional): if :obj:`True`, load the content of the local sqlite database from the external source
            max_entries (:obj:`float`, optional): maximum number of entries to save locally
            download_backup (:obj:`bool`, optional): if :obj:`True`, load the local copy of the data source from the Karr Lab server
            verbose (:obj:`bool`, optional): if :obj:`True`, print status information to the standard output
            clear_requests_cache (:obj:`bool`, optional): if :obj:`True`, clear the HTTP requests cache
        """

        """ CachedDataSource settings """
        if not name:
            name = self.__class__.__name__

        """ Request settings """
        self.requests_cache_filename = os.path.join(cache_dirname, name + '.requests.py{}.sqlite'.format(sys.version_info[0]))
        self.requests_session = self.get_requests_session()

        if clear_requests_cache:
            self.clear_requests_cache()

        self.disable_warnings()

        """ Call superclass constructor which will optionally load content """
        super(HttpDataSource, self).__init__(name=name, cache_dirname=cache_dirname,
                                             clear_content=clear_content, load_content=load_content, max_entries=max_entries,
                                             download_backup=download_backup, verbose=verbose)

    def get_requests_session(self):
        """ Setup an cache-enabled HTTP request session

        Returns:
            :obj:`requests_cache.core.CachedSession`: cached-enable session
        """
        if not os.path.isdir(os.path.dirname(self.requests_cache_filename)):
            os.makedirs(os.path.dirname(self.requests_cache_filename))
        name, _, _ = self.requests_cache_filename.rpartition('.')
        session = requests_cache.core.CachedSession(name, backend='sqlite', expire_after=None)
        return session

    def clear_requests_cache(self):
        """ Clear the cache-enabled HTTP request session """
        self.requests_session.cache.clear()

    def disable_warnings(self):
        """ Disable insecure HTTP request warnings """
        requests.packages.urllib3.disable_warnings(requests.packages.urllib3.exceptions.InsecureRequestWarning)


class WebserviceDataSource(DataSource):
    """ A data source that is a webservice """
    pass
