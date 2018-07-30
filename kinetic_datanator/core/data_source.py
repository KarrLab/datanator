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
from sqlalchemy_utils.functions import database_exists, create_database
import sys
import tarfile
import psycopg2

CACHE_DIRNAME = os.path.join(os.path.dirname(__file__), '..', 'data_source', 'cache')
# :obj:`str`: default path for the sqlite database


class DataSource(six.with_metaclass(abc.ABCMeta, object)):
    """ Represents an external data source

    Attributes:
        name (:obj:`str`): name
    """

    def __init__(self, name=None, verbose=False):
        """
        Args:
            name (:obj:`str`, optional): name
        """
        if not name:
            name = self.__class__.__name__
        self.name = name
        self.verbose = verbose

    def vprint(self, str):
        if self.verbose:
            print(str)

class PostgresDataSource(DataSource):
    """ Represents a Postgres database

    Need to have the data source be able to dump and restore from a dump

    """

    def __init__(self, cache_dirname=None,name=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, restore_backup=True, verbose=False):

        super(PostgresDataSource, self).__init__(name=name)

        self.app = create_app()
        # max entries
        self.max_entries = max_entries
        # committing
        self.commit_intermediate_results = commit_intermediate_results

        self.engine = self.get_engine()
        if clear_content:
            self.clear_content()
        self.session = self.get_session()
        if load_content:
            self.load_content()

    def get_engine(self):
        """ Get an engine for the postgres database. If the database doesn't exist, initialize its structure.

        Returns:
            :obj:`sqlalchemy.engine.Engine`: database engine
        """

        engine = self.base_model.engine
        if not database_exists(engine.url):
            create_database(engine.url)
            self.base_model.metadata.create_all(engine)

        return engine

    def clear_content(self):
        """ Clear the content of the sqlite database (i.e. drop and recreate all tables). """
        self.base_model.drop_all()
        self.base_model.create_all()

    def get_session(self):
        """ Get a session for the sqlite database

        Returns:
            :obj:`sqlalchemy.orm.session.Session`: database session
        """

        return self.base_model.session


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
        print(q)
        self.session.flush()
        if q.count():
            return q.first()
        else:
            obj = cls(**kwargs)
            self.session.add(obj)
            return obj



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

        base_model (:obj:`Base`): base ORM model for the sqlite databse
    """

    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backups=True, verbose=False):
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
            self.app.config['WHOOSH_BASE'] = self.cache_dirname + '/whoosh_cache/'
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
        """ Backup the local sqlite database to the Karr Lab server """
        for a_backup in self.get_backups(set_metadata=True):
            backup.BackupManager() \
                .create(a_backup) \
                .upload(a_backup) \
                .cleanup(a_backup)

    def download_backups(self):
        """ Download the local sqlite database from the Karr Lab server """
        for a_backup in self.get_backups(download=True):
            backup_manager = backup.BackupManager()
            backup_manager \
                .download(a_backup) \
                .extract(a_backup) \
                .cleanup(a_backup)

    def get_backups(self, download=False, set_metadata=False):
        """ Get a list of the files to backup/unpack

        Args:
            download (:obj:`bool`, optional): if :obj:`True`, prepare the files for uploading
            set_metadata (:obj:`bool`, optional): if :obj:`True`, set the metadata of the backup files

        Returns:
            :obj:`list` of :obj:`backup.Backup`: backups
        """
        list_backups = []
        a_backup = backup.Backup()
        path = backup.BackupPath(self.filename, self.name + '.sqlite')
        a_backup.paths.append(path)
        a_backup.local_filename = os.path.join(os.path.dirname(self.filename), path.arc_path + '.tar.gz')
        a_backup.remote_filename = path.arc_path + '.tar.gz'

        if set_metadata:
            a_backup.set_username_ip_date()
            a_backup.set_package(os.path.join(os.path.dirname(__file__), '..', '..'))

        if self.flask:
            whoosh_backup = backup.Backup()
            path = backup.BackupPath(self.cache_dirname + '/whoosh_cache/', 'whoosh_cache')
            whoosh_backup.paths.append(path)
            whoosh_backup.local_filename = path.arc_path
            whoosh_backup.remote_filename = path.arc_path + '.tar.gz'
            if set_metadata:
                whoosh_backup.set_username_ip_date()
                whoosh_backup.set_package(os.path.join(os.path.dirname(__file__), '..', '..'))
            list_backups.append(whoosh_backup)

        list_backups.append(a_backup)

        return list_backups

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
                 clear_requests_cache=False, download_request_backup=False, flask=False):
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
                                             download_backups=download_backups, verbose=verbose, flask=flask)

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

    def get_backups(self, download=False, set_metadata=False):
        """ Get a list of the files to backup/unpack

        Args:
            download (:obj:`bool`, optional): if :obj:`True`, prepare the files for uploading
            set_metadata (:obj:`bool`, optional): if :obj:`True`, set the metadata of the backup files

        Returns:
            :obj:`list` of :obj:`backup.Backup`: backups
        """
        backups = super(HttpDataSource, self).get_backups(download=download, set_metadata=set_metadata)
        if download and not self.download_request_backup:
            return backups

        requests_cache_basename = self.name + '.requests.py{}.sqlite'.format(sys.version_info[0])
        requests_cache_filename = os.path.join(os.path.dirname(self.requests_cache_filename), requests_cache_basename)
        a_backup = backup.Backup()
        path = backup.BackupPath(requests_cache_filename, requests_cache_basename)
        a_backup.paths.append(path)
        a_backup.local_filename = os.path.join(os.path.dirname(self.filename), path.arc_path + '.tar.gz')
        a_backup.remote_filename = path.arc_path + '.tar.gz'
        if set_metadata:
            a_backup.set_username_ip_date()
            a_backup.set_package(os.path.join(os.path.dirname(__file__), '..', '..'))

        backups.append(a_backup)
        return backups


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
