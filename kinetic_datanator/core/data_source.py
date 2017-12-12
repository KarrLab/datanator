"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-08
:Copyright: 2017, Karr Lab
:License: MIT
"""

from wc_utils import backup
import abc
import kinetic_datanator.config.paths
import os
import requests
import requests_cache
import six
import sqlalchemy
import sqlalchemy.orm
import sys
import tarfile
import wc_utils.config.core



config_manager = wc_utils.config.core.ConfigManager(kinetic_datanator.config.paths.core)
# :obj:`wc_utils.config.core.ConfigManager`: configuration manager

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
        backup_server_hostname (:obj:`str`): hostname for server to upload/download copies of the data source to/from the Karr Lab server
        backup_server_username (:obj:`str`): username for server to upload/download copies of the data source to/from the Karr Lab server
        backup_server_password (:obj:`str`): password for server to upload/download copies of the data source to/from the Karr Lab server
        backup_server_remote_dirname (:obj:`str`): remote directory on server to upload/download copies of the data source to/from the Karr Lab server
        engine (:obj:`sqlalchemy.engine.Engine`): sqlalchemy engine
        session (:obj:`sqlalchemy.orm.session.Session`): sqlalchemy session
        max_entries (:obj:`float`): maximum number of entries to save locally
        commit_intermediate_results (:obj:`bool`): if :obj:`True`, commit the changes throughout the loading
            process. This is particularly helpful for restarting this method when webservices go offline.
        verbose (:obj:`bool`): if :obj:`True`, print status information to the standard output

        base_model (:obj:`Base`): base ORM model for the sqlite databse
    """

    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backup=True, verbose=False, flask=False):
        """
        Args:
            name (:obj:`str`, optional): name
            cache_dirname (:obj:`str`, optional): directory to store the local copy of the data source
            clear_content (:obj:`bool`, optional): if :obj:`True`, clear the content of the sqlite local copy of the data source
            load_content (:obj:`bool`, optional): if :obj:`True`, load the content of the local sqlite database from the external source
            max_entries (:obj:`float`, optional): maximum number of entries to save locally
            commit_intermediate_results (:obj:`bool`, optional): if :obj:`True`, commit the changes throughout the loading
                process. This is particularly helpful for restarting this method when webservices go offline.
            download_backup (:obj:`bool`, optional): if :obj:`True`, load the local copy of the data source from the Karr Lab server
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

        # backup settings
        self.backup_server_hostname = config_manager.get_config()['wc_utils']['backup']['hostname']
        self.backup_server_username = config_manager.get_config()['wc_utils']['backup']['username']
        self.backup_server_password = config_manager.get_config()['wc_utils']['backup']['password']
        self.backup_server_remote_dirname = config_manager.get_config()['wc_utils']['backup']['remote_dirname']

        # verbosity
        self.verbose = verbose

        #flaskosity
        self.flask = flask

        """ Create SQLAlchemy session and load content if necessary """
        if os.path.isfile(self.filename):
            self.engine = self.get_engine()
            if clear_content:
                self.clear_content()
            self.session = self.get_session()
            if load_content:
                self.load_content()
        elif download_backup and self.backup_server_hostname:
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

        if self.flask:
            self.app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///' + self.filename
            self.app.config['WHOOSH_BASE'] = self.cache_dirname +'/whoosh_cache/'
            engine = self.base_model.create_all()
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

    def upload_backup(self):
        """ Backup the local sqlite database to the Karr Lab server """
        for file in self.get_backup_files(set_metadata=True):
            backup.BackupManager(
                os.path.join(os.path.dirname(self.filename), file.arcname + '.tar.gz'),
                file.arcname + '.tar.gz',
                hostname=self.backup_server_hostname, username=self.backup_server_username,
                password=self.backup_server_password, remote_dirname=self.backup_server_remote_dirname) \
                .create([file]) \
                .upload() \
                .cleanup()

    def download_backup(self):
        """ Download the local sqlite database from the Karr Lab server """
        for file in self.get_backup_files(download=True):
            backup_manager = backup.BackupManager(
                os.path.join(os.path.dirname(self.filename), file.arcname + '.tar.gz'),
                file.arcname + '.tar.gz',
                hostname=self.backup_server_hostname, username=self.backup_server_username,
                password=self.backup_server_password, remote_dirname=self.backup_server_remote_dirname)
            backup_manager.download()

        for file in self.get_backup_files(download=True, get_metadata=True):
            backup_manager = backup.BackupManager(
                os.path.join(os.path.dirname(self.filename), file.arcname + '.tar.gz'),
                file.arcname + '.tar.gz',
                hostname=self.backup_server_hostname, username=self.backup_server_username,
                password=self.backup_server_password, remote_dirname=self.backup_server_remote_dirname)
            backup_manager.extract([file])
            backup_manager.cleanup()

    def get_backup_files(self, download=False, get_metadata=False, set_metadata=False):
        """ Get a list of the files to backup/unpack

        Args:
            download (:obj:`bool`, optional): if :obj:`True`, prepare the files for uploading
            get_metadata (:obj:`bool`, optional): if :obj:`True`, get the metadata of the backup files
            set_metadata (:obj:`bool`, optional): if :obj:`True`, set the metadata of the backup files

        Returns:
            :obj:`list` of :obj:`backup.BackupFile`: list of files to backup/unpack
        """
        files = []

        file = backup.BackupFile(self.filename, self.name + '.sqlite')
        if set_metadata:
            file.set_created_modified_time()
            file.set_username_ip()
            file.set_program_version_from_repo(os.path.join(os.path.dirname(__file__), '..', '..'))
        files.append(file)

        return files

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
        if self.flask:
            q = cls.query.filter_by(**kwargs)
            if q.count():
                return q.first()
            else:
                obj = cls(**kwargs)
                self.session.add(obj)
                self.session.commit()
                return obj
        else:
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
        ENDPOINT_DOMAINS (:obj:`dict` of :obj:`str`, :obj:`str`): dictionary of domains to retry
        MAX_HTTP_RETRIES (:obj:`int`): maximum number of times to retry each HTTP request
    """

    ENDPOINT_DOMAINS = {}
    MAX_HTTP_RETRIES = 5

    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backup=True, verbose=False,
                 clear_requests_cache=False, download_request_backup=False, flask = False):
        """
        Args:
            name (:obj:`str`, optional): name
            cache_dirname (:obj:`str`, optional): directory to store the local copy of the data source and the HTTP requests cache
            clear_content (:obj:`bool`, optional): if :obj:`True`, clear the content of the sqlite local copy of the data source
            load_content (:obj:`bool`, optional): if :obj:`True`, load the content of the local sqlite database from the external source
            max_entries (:obj:`float`, optional): maximum number of entries to save locally
            commit_intermediate_results (:obj:`bool`, optional): if :obj:`True`, commit the changes throughout the loading
                process. This is particularly helpful for restarting this method when webservices go offline.
            download_backup (:obj:`bool`, optional): if :obj:`True`, load the local copy of the data source from the Karr Lab server
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
                                             download_backup=download_backup, verbose=verbose, flask = flask)

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

    def get_backup_files(self, download=False, get_metadata=False, set_metadata=False):
        """ Get a list of the files to backup/unpack

        Args:
            download (:obj:`bool`, optional): if :obj:`True`, prepare the files for uploading
            get_metadata (:obj:`bool`, optional): if :obj:`True`, get the metadata of the backup files
            set_metadata (:obj:`bool`, optional): if :obj:`True`, set the metadata of the backup files

        Returns:
            :obj:`list` of :obj:`backup.BackupFile`: list of files to backup/unpack
        """
        files = super(HttpDataSource, self).get_backup_files(download=download, get_metadata=get_metadata, set_metadata=set_metadata)
        if download and not self.download_request_backup:
            return files

        if get_metadata:
            tar_file = tarfile.open(self.filename + '.tar.gz', "r:gz")

        requests_cache_basename = self.name + '.requests.py{}.sqlite'.format(sys.version_info[0])
        requests_cache_filename = os.path.join(os.path.dirname(self.requests_cache_filename), requests_cache_basename)
        file = backup.BackupFile(requests_cache_filename, requests_cache_basename)
        if get_metadata:
            try:
                tar_file.getmember(file.arcname)
            except KeyError:
                pass
        elif set_metadata and os.path.isfile(file.filename):
            file.set_created_modified_time()
            file.set_username_ip()
            file.set_program_version_from_repo(os.path.join(os.path.dirname(__file__), '..', '..'))
        files.append(file)

        if get_metadata:
            tar_file.close()

        return files


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
