""" Downloads and parses the ArrayExpress database

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-08-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

import datetime
import dateutil.parser
import pkg_resources
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
from kinetic_datanator.core import data_source


Base = sqlalchemy.ext.declarative.declarative_base()
# :obj:`Base`: base model for local sqlite database

sample_characteristic = sqlalchemy.Table(
    'sample_characteristic', Base.metadata,
    sqlalchemy.Column('sample__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('sample._id'), index=True),
    sqlalchemy.Column('characteristic__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('characteristic._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Sample:Characteristic many-to-many association table

sample_variable = sqlalchemy.Table(
    'sample_variable', Base.metadata,
    sqlalchemy.Column('sample__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('sample._id'), index=True),
    sqlalchemy.Column('variable__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('variable._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Sample:Variable many-to-many association table

experiment_organism = sqlalchemy.Table(
    'experiment_organism', Base.metadata,
    sqlalchemy.Column('experiment__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment._id'), index=True),
    sqlalchemy.Column('organism__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('organism._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Experiment:Organism many-to-many association table

experiment_experiment_type = sqlalchemy.Table(
    'experiment_experiment_type', Base.metadata,
    sqlalchemy.Column('experiment__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment._id'), index=True),
    sqlalchemy.Column('experiment_type__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment_type._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Experiment:ExperimentType many-to-many association table

experiment_experiment_design = sqlalchemy.Table(
    'experiment_experiment_design', Base.metadata,
    sqlalchemy.Column('experiment__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment._id'), index=True),
    sqlalchemy.Column('experiment_design__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment_design._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Experiment:ExperimentDesign many-to-many association table

experiment_data_format = sqlalchemy.Table(
    'experiment_data_format', Base.metadata,
    sqlalchemy.Column('experiment__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment._id'), index=True),
    sqlalchemy.Column('data_format__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('data_format._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Experiment:DataFormat many-to-many association table

extract_sample = sqlalchemy.Table(
    'extract_sample', Base.metadata,
    sqlalchemy.Column('extract__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('extract._id'), index=True),
    sqlalchemy.Column('sample__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('sample._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Extract:Sample many-to-many association table

experiment_protocol = sqlalchemy.Table(
    'experiment_protocol', Base.metadata,
    sqlalchemy.Column('experiment__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('experiment._id'), index=True),
    sqlalchemy.Column('protocol_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('protocol._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Experiment:Protocol many-to-many association table


class Characteristic(Base):
    """ Represents an experimental characteristic

    Attributes:
        _id (:obj:`int`): unique id
        category (:obj:`str`): name of the characteristic (e.g. organism)
        value (:obj:`str`): value of characteristic (e.g. Mus musculus)
        samples (:obj:`list` of :obj:`Sample`): samples
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    category = sqlalchemy.Column(sqlalchemy.String())
    value = sqlalchemy.Column(sqlalchemy.String())

    sqlalchemy.schema.UniqueConstraint(category, value)

    __tablename__ = 'characteristic'


class Variable(Base):
    """ Represents an experimental variable

    Attributes:
        _id (:obj:`int`): unique id
        name (:obj:`str`): name of the variable (e.g. genotype)
        value (:obj:`str`): value of variable (e.g control)
        unit (:obj:`str`): units of value (e.g control). This field is not always filled. 
        samples (:obj:`list` of :obj:`Sample`): samples
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String())
    value = sqlalchemy.Column(sqlalchemy.String())
    unit = sqlalchemy.Column(sqlalchemy.String())

    sqlalchemy.schema.UniqueConstraint(name, value, unit)

    __tablename__ = 'variable'


class Sample(Base):
    """ Represents an observed concentration

    Attributes:
        _id (:obj:`int`): unique id
        experiment_id: (:obj:`int`): the accesion number of the experiment
        experiment (:obj:`Experiment`): experiment that the sample belongs to
        index (:obj:`int`): index of the sample within the experiment
        name (:obj:`str`): name of the source of the sample
        assay (:obj:`str`): assay
        characteristics (:obj:`list` of :obj:`Characteristic'): characteristics
        variables (:obj:`list` of :obj:`Variable'): variables
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    experiment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('experiment._id'), index=True)
    experiment = sqlalchemy.orm.relationship('Experiment', backref=sqlalchemy.orm.backref('samples'), foreign_keys=[experiment_id])
    index = sqlalchemy.Column(sqlalchemy.Integer())
    name = sqlalchemy.Column(sqlalchemy.String())
    assay = sqlalchemy.Column(sqlalchemy.String())
    characteristics = sqlalchemy.orm.relationship(
        'Characteristic', secondary=sample_characteristic, backref=sqlalchemy.orm.backref('samples'))
    variables = sqlalchemy.orm.relationship('Variable', secondary=sample_variable, backref=sqlalchemy.orm.backref('samples'))

    sqlalchemy.schema.UniqueConstraint(experiment_id, index)

    __tablename__ = 'sample'


class Experiment(Base):
    """ Represents an experiment

    Attributes:
        _id (:obj:`int`): unique id
        id (:obj:`str`): unique string identifier assigned by ArrayExpress
        name (:obj:`str`): name
        name_2 (:obj:`str`): second name
        description (:obj:`str`): description
        organisms (:obj:`list` of :obj:`Organism`): list of organisms
        types (:obj:`list` of :obj:`ExperimentType`): list of experiment types
        designs (:obj:`list` of :obj:`ExperimentDesign`): list of experimental designs
        submission_date (:obj:`datetime.date`): submission date
        release_date (:obj:`datetime.date`): release date
        data_formats (:obj:`list` of :obj:`DataFormat`): list of data formats
        samples (:obj:`list` of :obj:`Sample`): list of samples
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    id = sqlalchemy.Column(sqlalchemy.String(), unique=True, index=True)
    name = sqlalchemy.Column(sqlalchemy.String())
    name_2 = sqlalchemy.Column(sqlalchemy.String())
    description = sqlalchemy.Column(sqlalchemy.String())
    organisms = sqlalchemy.orm.relationship('Organism', secondary=experiment_organism, backref=sqlalchemy.orm.backref('experiments'))
    types = sqlalchemy.orm.relationship('ExperimentType', secondary=experiment_experiment_type,
                                        backref=sqlalchemy.orm.backref('experiments'))
    designs = sqlalchemy.orm.relationship('ExperimentDesign', secondary=experiment_experiment_design,
                                          backref=sqlalchemy.orm.backref('experiments'))
    submission_date = sqlalchemy.Column(sqlalchemy.Date())
    release_date = sqlalchemy.Column(sqlalchemy.Date())
    data_formats = sqlalchemy.orm.relationship('DataFormat', secondary=experiment_data_format,
                                               backref=sqlalchemy.orm.backref('experiments'))

    __tablename__ = 'experiment'


class ExperimentDesign(Base):
    """ Represents and experimental design

    Attributes:
        _id (:obj:`int`): unique id
        name (:obj:`str`): name
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String(), unique=True)

    __tablename__ = 'experiment_design'


class ExperimentType(Base):
    """ Represents a type of experiment

    Attributes:
        _id (:obj:`int`): unique id
        name (:obj:`str`): name
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String(), unique=True)

    __tablename__ = 'experiment_type'


class DataFormat(Base):
    """ Represents a data format

    Attributes:
        _id (:obj:`int`): unique id    
        name (:obj:`str`): name
        experiments (:obj:`list` of :obj:`Experiment`): list of experiments
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String(), unique=True)

    __tablename__ = 'data_format'


class Organism(Base):
    """ Represents an organism

    Attributes:
        _id (:obj:`int`): unique id
        name (:obj:`str`): name
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String(), unique=True)

    __tablename__ = 'organism'


class Extract(Base):
    """ Represents an extract of a sample

    Attributes:
        _id (:obj:`int`): unique id
        name (:obj:`str`): name
        samples (:obj:`list` of :obj:`Sample`): list of samples
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String(), unique=True)
    samples = sqlalchemy.orm.relationship('Sample', secondary=extract_sample, backref=sqlalchemy.orm.backref('extracts'))

    __tablename__ = 'extract'


class Protocol(Base):
    """ Represents a protocol for an experiment

    Attributes:
        _id (:obj:`int`): unique id
        protocol_accession (:obj:`str`): array express identifier for protocol
        protocol_type (:obj:`list` of :obj:`Sample`): the type of exerpimental protocol (e.g. normalization, extraction, etc.)
        text (:obj:`str`): description the protocol
        experiments (:obj:`list` of :obj:`Experiment`): list of experiments that performed this protocol
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    protocol_accession = sqlalchemy.Column(sqlalchemy.String())
    protocol_type = sqlalchemy.Column(sqlalchemy.String())
    text = sqlalchemy.Column(sqlalchemy.String())
    experiments = sqlalchemy.orm.relationship('Experiment', secondary=experiment_protocol, backref=sqlalchemy.orm.backref('protocols'))
    
    __tablename__ = 'protocol'


class ArrayExpress(data_source.HttpDataSource):
    """ A local sqlite copy of the ArrayExpress database

    Attributes:
        EXCLUDED_DATASET_IDS (:obj:`list` of :obj:`str`): list of IDs of datasets to exclude
    """

    base_model = Base

    ENDPOINT_DOMAINS = {
        'array_express': 'https://www.ebi.ac.uk/arrayexpress/json/v3/experiments',
    }

    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backup=True, verbose=False,
                 clear_requests_cache=False, download_request_backup=False):
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
        super(ArrayExpress, self).__init__(name=name, cache_dirname=cache_dirname, clear_content=clear_content,
                                           load_content=load_content, max_entries=max_entries,
                                           commit_intermediate_results=commit_intermediate_results,
                                           download_backup=download_backup, verbose=verbose,
                                           clear_requests_cache=clear_requests_cache, download_request_backup=download_request_backup)

        with open(pkg_resources.resource_filename('kinetic_datanator', 'data_source/array_express_excluded_dataset_ids.txt'), 'r') as file:
            self.EXCLUDED_DATASET_IDS = [line.rstrip() for line in file]

    def load_content(self, start_year=2001, end_year=None):
        """
        Downloads all medatata from array exrpess on their samples and experiments. The metadata
        is saved as the text file. Within the text files, the data is stored as a JSON object. 

        Args:
            start_year (:obj:`int`, optional): the first year to retrieve experiments for
            end_year (:obj:`int`, optional): the last year to retrieve experiments for
        """
        if not end_year:
            end_year = datetime.datetime.now().year

        session = self.session

        # download and parse experiment ids
        if self.verbose:
            print('Loading metadata for experiments ...')

        self.load_experiment_metadata(start_year, end_year)

        if self.verbose:
            print('  done.')

        # download and parse experiments
        n_experiments = session.query(Experiment).count()

        if self.verbose:
            print('Loading samples and protocols for experiments ...')

        for i_experiment, experiment in enumerate(session.query(Experiment).all()):
            if self.verbose and i_experiment % 500 == 0:
                print('  Loading samples and protocols for experiment {} of {}'.format(i_experiment + 1, n_experiments))

            self.load_experiment_samples(experiment)
            self.load_experiment_protocols(experiment)

        if self.verbose:
            print('  done.')

        self.session.commit()

    def load_experiment_metadata(self, start_year=2001, end_year=None):
        """ Get a list of accession identifiers for the experiments from the year :obj:`start_year` to year :obj:`end_year`

        Args:
            start_year (:obj:`int`, optional): the first year to retrieve experiment acession ids for
            end_year (:obj:`int`, optional): the last year to retrieve experiment acession ids for

        Returns:
            :obj:`list` of :obj:`str`: list of experiment accession identifiers
        """
        if not end_year:
            end_year = datetime.datetime.now().year

        db_session = self.session

        for year in range(start_year, end_year + 1):
            response = self.requests_session.get(self.ENDPOINT_DOMAINS['array_express'] + '?date=[{}-01-01+{}-12-31]'.format(year, year))
            response.raise_for_status()
            for expt_json in response.json()['experiments']['experiment']:

                id = expt_json['accession']
                if id in self.EXCLUDED_DATASET_IDS:
                    continue

                experiment = self.get_or_create_object(Experiment, id=id)

                if isinstance(expt_json['name'], list):
                    experiment.name = expt_json['name'][0]
                    experiment.name_2 = expt_json['name'][1]
                else:
                    experiment.name = expt_json['name']

                if 'organism' in expt_json:
                    for organism_name in expt_json['organism']:
                        experiment.organisms.append(self.get_or_create_object(Organism, name=organism_name))

                if 'description' in expt_json:
                    experiment.description = expt_json['description'][0]['text']

                if 'experimenttype' in expt_json:
                    entries = expt_json['experimenttype']
                    for entry in entries:
                        experiment.types.append(self.get_or_create_object(ExperimentType, name=entry))

                if 'experimentdesign' in expt_json:
                    entries = expt_json['experimentdesign']
                    for entry in entries:
                        experiment.designs.append(self.get_or_create_object(ExperimentDesign, name=entry))

                if 'bioassaydatagroup' in expt_json:
                    for entry in expt_json['bioassaydatagroup']:
                        experiment.data_formats.append(self.get_or_create_object(DataFormat, name=entry['dataformat']))

                if 'submissiondate' in expt_json:
                    experiment.submission_date = dateutil.parser.parse(expt_json['submissiondate']).date()

                if 'releasedate' in expt_json:
                    experiment.release_date = dateutil.parser.parse(expt_json['releasedate']).date()

                db_session.add(experiment)

    def load_experiment_samples(self, experiment):
        """ Load the samples for an experiment

        Args:
            experiment (:obj:`Experiment`): experiment
        """
        response = self.requests_session.get(self.ENDPOINT_DOMAINS['array_express'] + "/{}/samples".format(experiment.id))
        response.raise_for_status()
        json = response.json()

        if 'experiment' not in json:
            return
        experiment_json = json['experiment']

        if 'sample' not in experiment_json:
            return
        samples = experiment_json['sample']

        session = self.session

        if not isinstance(samples, list):
            samples = [samples]

        for i_sample, sample in enumerate(samples):
            self.load_experiment_sample(experiment, sample, i_sample)

    def load_experiment_sample(self, experiment, sample_json, index):
        """ Load the samples for an experiment

        Args:
            experiment (:obj:`Experiment`): experiment
            sample_json (:obj:`dict`): sample
            index (:obj:`int`): index of the sample within the experiment
        """
        sample = Sample()

        experiment.samples.append(sample)
        sample.index = index

        if 'source' in sample_json:
            source = sample_json['source']
            if isinstance(source, list):
                sample.name = source[0]['name']
            else:
                sample.name = source['name']

        if 'extract' in sample_json:
            extracts = sample_json['extract']
            if not isinstance(extracts, list):
                extracts = [extracts]
            for extract in extracts:
                sample.extracts.append(self.get_or_create_object(Extract, name=extract['name']))

        if 'assay' in sample_json:
            sample.assay = sample_json['assay']['name']

        if 'characteristic' in sample_json:
            characteristics = sample_json['characteristic']
            if not isinstance(characteristics, list):
                characteristics = [characteristics]
            for characteristic in characteristics:
                sample.characteristics.append(self.get_or_create_object(
                    Characteristic, category=characteristic['category'], value=characteristic['value']))

        if 'variable' in sample_json:
            variables = sample_json['variable']
            if not isinstance(variables, list):
                variables = [variables]
            for variable in variables:
                unit = None
                if "unit" in variable:
                    unit = variable['unit']
                sample.variables.append(self.get_or_create_object(
                    Variable, name=variable['name'], value=variable['value'], unit=unit))

    
    def load_experiment_protocols(self, experiment):

        """ Load the protocols for an experiment

        Args:
            experiment (:obj:`Experiment`): experiment
        """
        response = self.requests_session.get(self.ENDPOINT_DOMAINS['array_express'] + "/{}/protocols".format(experiment.id))
        response.raise_for_status()
        json = response.json()

        session = self.session
        # todo: implement

        if 'protocols' not in json:
            return
        protocol_json = json['protocols']

        if 'protocol' not in protocol_json:
            return
        protocols = protocol_json['protocol']

        session = self.session

        if not isinstance(protocols, list):
            protocols = [protocols]

        for protocol in protocols:
            self.load_experiment_protocol(experiment, protocol)


    def load_experiment_protocol(self, experiment, protocol_json):
        
        """ Load the protocols for an experiment

        Args:
            experiment (:obj:`Experiment`): experiment
            protocol_json (:obj:`dict`): sample
        """

        db_session = self.session
        
        protocol = Protocol()

        if 'accession' not in protocol_json:
            return
        protocol_accession = protocol_json['accession']
        protocol = self.get_or_create_object(Protocol, protocol_accession=protocol_accession)

        if 'type' in protocol_json:
            protocol.protocol_type = protocol_json['type']
        if 'text' in protocol_json:
            if isinstance(protocol_json['text'], str):
                protocol.text = protocol_json['text']
            if isinstance(protocol_json['text'], list):
                details = ""
                for item in protocol_json['text']:
                    if isinstance(item, unicode):
                        if item[-1:] == ',':
                            item = item[:-1]
                        details = details + item + "\n"
                if details:
                    details = details[:-1]
                protocol.text = details
        
        protocol.experiments.append(experiment)

    def get_or_create_object(self, cls, **kwargs):
        """ Get the first instance of :obj:`cls` that has the property-values pairs described by kwargs, or create an instance of :obj:`cls`
        if there is no instance with the property-values pairs described by kwargs

        Args:
            cls (:obj:`class`): type of object to find or create
            **kwargs: values of the properties of the object

        Returns:
            :obj:`Base`: instance of :obj:`cls` hat has the property-values pairs described by kwargs
        """
        q = self.session.query(cls).filter_by(**kwargs)
        if self.session.query(q.exists()).scalar():
            return q.first()

        obj = cls(**kwargs)
        self.session.add(obj)
        return obj
