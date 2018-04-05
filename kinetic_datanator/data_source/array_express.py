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
from kinetic_datanator.data_source.array_express_tools import ensembl_tools
import requests
import time
from ete3 import NCBITaxa
import os


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

sample_url = sqlalchemy.Table(
    'sample_url', Base.metadata,
    sqlalchemy.Column('sample__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('sample._id'), index=True),
    sqlalchemy.Column('url_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('url._id'), index=True),
)

sample_ensemblInfo = sqlalchemy.Table(
    'sample_ensemblInfo', Base.metadata,
    sqlalchemy.Column('sample__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('sample._id'), index=True),
    sqlalchemy.Column('ensemblInfo_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('ensemblInfo._id'), index=True),
)


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


class Url(Base):
    """ Represents a url
    Attributes:
        _id (:obj:`int`): unique id
        url (:obj:`str`): the text of the url
        samples (:obj:`list` of :obj:`Sample`): samples
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    url = sqlalchemy.Column(sqlalchemy.String())
    sqlalchemy.schema.UniqueConstraint(url)

    __tablename__ = 'url'


class EnsemblInfo(Base):
    """ Represents a url
    Attributes:
        _id (:obj:`int`): unique id
        organism_strain (:obj:`str`): the particular strain that relates to the ensembl reference genome (e.g. escherichia_coli_k12)
        url (:obj:`str`): the download url for the CDNA file from ensembl
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    organism_strain = sqlalchemy.Column(sqlalchemy.String(), index=True)
    url = sqlalchemy.Column(sqlalchemy.String())
    sqlalchemy.schema.UniqueConstraint(url)

    __tablename__ = 'ensemblInfo'


class Sample(Base):
    """ Represents an observed concentration
    Attributes:
        _id (:obj:`int`): unique id
        experiment_id: (:obj:`int`): the id of the experiment the samaple belongs in
        experiment (:obj:`Experiment`): experiment that the sample belongs to
        index (:obj:`int`): index of the sample within the experiment
        name (:obj:`str`): name of the source of the sample (this is used to identify the sample in arraya express)
        assay (:obj:`str`): name of the assay
        ensembl_organism_strain (:obj:`str`): the particular strain that relates to the ensembl reference genome (e.g. escherichia_coli_k12)
        characteristics (:obj:`list` of :obj:`Characteristic'): characteristics
        variables (:obj:`list` of :obj:`Variable'): variablesassay (:obj:`str`): name of the assay
        fastq_urls (:obj:`list` of :obj:`Url'): variablesassay (:obj:`str`): name of the assay
        read_type (:obj:`str`): the nature of the FASTQ file reads. Either 'single', 'multiple', or 'parallel'
        ensembl_info (:obj:`list` of :obj:`Variable'): informtation about the ensembl reference genome
        full_strain_specificity (:obj:`bool`): whether or not ensembl reference genome matches the full strain specifity recoreded in array express
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    experiment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('experiment.id'), index=True)
    experiment = sqlalchemy.orm.relationship('Experiment', backref=sqlalchemy.orm.backref('samples'), foreign_keys=[experiment_id])
    index = sqlalchemy.Column(sqlalchemy.Integer())
    name = sqlalchemy.Column(sqlalchemy.String())
    assay = sqlalchemy.Column(sqlalchemy.String())
    ensembl_organism_strain = sqlalchemy.Column(sqlalchemy.String())
    characteristics = sqlalchemy.orm.relationship(
        'Characteristic', secondary=sample_characteristic, backref=sqlalchemy.orm.backref('samples'))
    variables = sqlalchemy.orm.relationship('Variable', secondary=sample_variable, backref=sqlalchemy.orm.backref('samples'))
    fastq_urls = sqlalchemy.orm.relationship('Url', secondary=sample_url, backref=sqlalchemy.orm.backref('samples'))
    read_type = sqlalchemy.Column(sqlalchemy.String())
    ensembl_info = sqlalchemy.orm.relationship('EnsemblInfo', secondary=sample_ensemblInfo, backref=sqlalchemy.orm.backref('samples'))
    full_strain_specificity = sqlalchemy.Column(sqlalchemy.Boolean())
    sqlalchemy.schema.UniqueConstraint(experiment_id, name)
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
        read_type (:obj:`str`): type of FASTQ files (in an RNA-Seq experiment)
        has_fastq_files (:obj:`bool`): whether this experiment has FASTQ files or not
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
    read_type = sqlalchemy.Column(sqlalchemy.String())
    has_fastq_files = sqlalchemy.Column(sqlalchemy.Boolean())

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
        bio_assay_data_cubes (:obj:`int`): number of dimensions to the data
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String())
    bio_assay_data_cubes = sqlalchemy.Column(sqlalchemy.Integer())

    sqlalchemy.schema.UniqueConstraint(name, bio_assay_data_cubes)
    __tablename__ = 'data_format'


class Organism(Base):
    """ Represents an organism
    Attributes:
        _id (:obj:`int`): unique id
        name (:obj:`str`): name
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String(), unique=True)
    ncbi_id = sqlalchemy.Column(sqlalchemy.Integer())

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
        performer (:obj:`str`): name of the person who did the experiment
        hardware (:obj:`str`): hardware (usually detection instruments) used in protocol
        software (:obj:`str`): software (usually for analyzing and normalizing the data)
        experiments (:obj:`list` of :obj:`Experiment`): list of experiments that performed this protocol
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    protocol_accession = sqlalchemy.Column(sqlalchemy.String())
    protocol_type = sqlalchemy.Column(sqlalchemy.String())
    text = sqlalchemy.Column(sqlalchemy.String())
    performer = sqlalchemy.Column(sqlalchemy.String())
    hardware = sqlalchemy.Column(sqlalchemy.String())
    software = sqlalchemy.Column(sqlalchemy.String())
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
                 commit_intermediate_results=False, download_backups=True, verbose=False,
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
        #with open(pkg_resources.resource_filename('kinetic_datanator', 'data_source/array_express_excluded_dataset_ids.txt'), 'r') as file:
        #    self.EXCLUDED_DATASET_IDS = [line.rstrip() for line in file]
        super(ArrayExpress, self).__init__(name=name, cache_dirname=cache_dirname, clear_content=clear_content,
                                           load_content=load_content, max_entries=max_entries,
                                           commit_intermediate_results=commit_intermediate_results,
                                           download_backups=download_backups, verbose=verbose,
                                           clear_requests_cache=clear_requests_cache, download_request_backup=download_request_backup)
        
        #with open(pkg_resources.resource_filename('kinetic_datanator', 'data_source/array_express_excluded_dataset_ids.txt'), 'r') as file:
        #    self.EXCLUDED_DATASET_IDS = [line.rstrip() for line in file]

    def load_content(self, test_url=""):
        """
        Downloads all medatata from array exrpess on their samples and experiments. The metadata
        is saved as the text file. Within the text files, the data is stored as a JSON object.
        Args:
            start_year (:obj:`int`, optional): the first year to retrieve experiments for
            end_year (:obj:`int`, optional): the last year to retrieve experiments for
        """

        session = self.session

        # download and parse experiment ids
        if self.verbose:
            print('Loading metadata for experiments ...')

        self.load_experiment_metadata(test_url)
        self.session.commit()

        if self.verbose:
            print('  done.')

        # download and parse experiments
        n_experiments = session.query(Experiment).count()

        if self.verbose:
            print('Loading samples and protocols for experiments ...')

       # download and parse experiments
        if self.verbose:
            print('Loading samples and protocols for experiments ...')

        total_experiments = session.query(Experiment).all()
        for i_experiment, experiment in enumerate(total_experiments):
            if self.verbose and i_experiment % 500 == 0:
                print('  Loading samples and protocols for experiment {} of {}'.format(i_experiment + 1, len(total_experiments)))
            self.load_experiment_samples(experiment)
            self.load_experiment_protocols(experiment)

        if self.verbose:
            print('  done.')
        self.session.commit()

    def load_experiment_metadata(self, test_url=""):
        """ Get a list of accession identifiers for the experiments from the year :obj:`start_year` to year :obj:`end_year`
        Args:
            start_year (:obj:`int`, optional): the first year to retrieve experiment acession ids for
            end_year (:obj:`int`, optional): the last year to retrieve experiment acession ids for
        Returns:
            :obj:`list` of :obj:`str`: list of experiment accession identifiers
        """

        db_session = self.session
        #?date=[{}-01-01+{}-12-31]
        if not test_url:
            response = self.requests_session.get(
                self.ENDPOINT_DOMAINS['array_express'] + '?exptype="RNA-seq of coding RNA from single cells"+OR+"RNA-seq of coding RNA"')
        elif test_url:
            response = self.requests_session.get(test_url)
        response.raise_for_status()
        for expt_json in response.json()['experiments']['experiment']:

            id = expt_json['accession']
            #if id in self.EXCLUDED_DATASET_IDS:
            #    continue

            experiment = self.get_or_create_object(Experiment, id=id)

            if isinstance(expt_json['name'], list):
                experiment.name = expt_json['name'][0]
                experiment.name_2 = expt_json['name'][1]
            else:
                experiment.name = expt_json['name']

            taxon_exceptions = {}
            exceptions_file = open("{}/array_express_tools/taxon_exceptions.txt".format(os.path.dirname(os.path.realpath(__file__))))
            for line in exceptions_file.readlines()[1:]:
                split = line.split(" -- ")
                taxon_exceptions[split[0]] = split[1][:-1]
            if 'organism' in expt_json:
                for organism_name in expt_json['organism']:
                    if organism_name[-1:] == " ":
                        organism_name = organism_name[:-1]
                    if organism_name not in taxon_exceptions:
                        experiment.organisms.append(self.get_or_create_object(Organism, name=organism_name, ncbi_id = NCBITaxa().get_name_translator([organism_name])[organism_name][0]))
                    else:
                        organisms = taxon_exceptions[organism_name].split(", ")
                        for org in organisms:
                            experiment.organisms.append(self.get_or_create_object(Organism, name=org, ncbi_id = NCBITaxa().get_name_translator([org])[org][0]))


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
                    experiment.data_formats.append(self.get_or_create_object(
                        DataFormat, name=entry['dataformat'], bio_assay_data_cubes=entry['bioassaydatacubes']))

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
        #while time.clock>1
        try:
            url = self.ENDPOINT_DOMAINS['array_express'] + "/{}/samples".format(experiment.id)
            response = self.requests_session.get(url)
            response.raise_for_status()
        except requests.HTTPError as resp:
            print(str(resp))
            if str(resp).startswith("500 Server Error: Internal Server Error for url:"):
                return 
            else:
                raise

        json = response.json()

        if 'experiment' not in json:
            return
        experiment_json = json['experiment']
        if 'sample' not in experiment_json:
            return
        samples = experiment_json['sample']

        paired_end = False
        list_of_source_names = []
        for sample in samples:
            if "source" in sample:
                list_of_source_names.append(sample['source']['name'])
        if ((len(samples) == len(set(list(list_of_source_names)))*2)  # if its paried end, then there two sample entries to every one source name
                # if its paired end, then each source name must appear twice
                and (list_of_source_names.count(list_of_source_names[0]) == 2)
                and (list_of_source_names.count(list_of_source_names[len(list_of_source_names)-1:len(list_of_source_names)][0]) == 2)
            ):
            paired_end = True

        if paired_end:
            experiment.read_type = "paired"
        else:
            experiment.read_type = "single"
        session = self.session

        if not isinstance(samples, list):
            samples = [samples]
        experiment.has_fastq_files = False

        sample_indeces = {}
        for name in set(list(list_of_source_names)):
            sample_indeces[name] = [i for i, j in enumerate(list_of_source_names) if j == name]
        i = 0
        for key, value in sample_indeces.items():
            new_sample = self.load_experiment_sample(experiment, samples[value[0]], i)
            new_sample.read_type = experiment.read_type
            i = i+1
            for num in value:
                if 'scan' in samples[num]:
                    if 'comment' in samples[num]['scan']:
                        if isinstance(samples[num]['scan']['comment'], list):
                            for comment in samples[num]['scan']['comment']:
                                if comment['name'] == "FASTQ_URI":
                                    if comment['value'] not in [u.url for u in new_sample.fastq_urls]:
                                        new_sample.fastq_urls.append(self.get_or_create_object(Url, url=comment['value']))
                                        experiment.has_fastq_files = True

                        elif isinstance(samples[num]['scan']['comment'], dict):
                            if sample['scan']['comment']['name'] == "FASTQ_URI":
                                if comment['value'] not in [u.url for u in new_sample.fastq_urls]:
                                    new_sample.fastq_urls.append(self.get_or_create_object(Url, url=comment['value']))
                                    experiment.has_fastq_files = True


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
                sample.name = source[0]['name'].encode('ascii', 'ignore')
            else:
                sample.name = source['name'].encode('ascii', 'ignore')

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
        
        try:
            strain_info = ensembl_tools.get_strain_info(sample)
        except LookupError:
            strain_info = None
        if strain_info:
            ensembl_info = self.session.query(EnsemblInfo).filter_by(organism_strain = strain_info.organism_strain).first()
            if ensembl_info:
                sample.ensembl_info.append(ensembl_info)
            else:
                try:
                    ftp_url = ensembl_tools.get_ftp_url(strain_info.download_url)
                except LookupError:
                    ftp_url = None
                if ftp_url != None:
                    sample.ensembl_info.append(
                        EnsemblInfo(organism_strain=strain_info.organism_strain, url=ftp_url))
                    sample.full_strain_specificity = strain_info.full_strain_specificity
            
        return sample

    def load_experiment_protocols(self, experiment):
        """ Load the protocols for an experiment
        Args:
            experiment (:obj:`Experiment`): experiment
        """
        response = self.requests_session.get(self.ENDPOINT_DOMAINS['array_express'] + "/{}/protocols".format(experiment.id))
        json = response.json()
        session = self.session
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
            if not(isinstance(protocol_json['text'], list) or isinstance(protocol_json['text'], dict)):
                protocol.text = protocol_json['text']
            if isinstance(protocol_json['text'], list):
                details = ""
                for item in protocol_json['text']:
                    if not(isinstance(item, list) or isinstance(item, dict)):
                        if item[-1:] == ',':
                            item = item[:-1]
                        details = details + item + "\n"
                if details:
                    details = details[:-1]
                protocol.text = details
        if 'performer' in protocol_json:
            protocol.performer = protocol_json['performer']
        if 'hardware' in protocol_json:
            protocol.hardware = protocol_json['hardware']
        if 'software' in protocol_json:
            protocol.software = protocol_json['software']
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
