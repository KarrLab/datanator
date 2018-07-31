"""
This code is a common schema for all the kinetic_datanator modules

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-31
:Copyright: 2017, Karr Lab
:License: MIT
"""
import os
import time
import re
import flask_whooshalchemy
import sqlalchemy.ext.declarative
from sqlalchemy import Column, BigInteger, Integer, Float, String, Text, ForeignKey, Boolean, Table,  Numeric, or_
from sqlalchemy.orm import relationship, backref, sessionmaker
from kinetic_datanator.util.build_util import timeit
from kinetic_datanator.config import config
from kinetic_datanator.core import data_source, models
from kinetic_datanator.app import create_app, db
from kinetic_datanator.util.constants import *
from kinetic_datanator.data_source import corum, pax, jaspar, jaspar, ecmdb, sabio_rk, intact, uniprot, array_express
from ete3 import NCBITaxa
from sqlalchemy.sql import func


class CommonSchema(data_source.PostgresDataSource):
    """
    A Local SQLlite copy of the aggregation of data_source modules

    """
    app = create_app()
    base_model = db
    text_indicies = [cls for name, cls in models.__dict__.items() if hasattr(cls,'__searchable__')]

    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, restore_backup=True, verbose=False, load_entire_small_DBs=False, test=True):
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
            verbose (:obj:`bool`, optional): if :obj:`True`, self.vprint status information to the standard output
            load_entire_small_DBs (:obj:`bool`, optional): Loads all entire databases that fall under 50 MB
            flask (:obj:`bool`, optional): Designates whether the database is defined as a Flask models
            test (:obj:`bool`, optional): Designates whether tests are being completed for brevity of tests
        """

        super(CommonSchema, self).__init__(name=name, cache_dirname=cache_dirname, clear_content=clear_content,
                                                load_content=False, max_entries=max_entries,
                                                commit_intermediate_results=commit_intermediate_results,
                                                restore_backup=restore_backup, verbose=verbose)

        self.load_entire_small_DBs = load_entire_small_DBs
        self.test = test

        if load_content:
            self.get_or_create_object(models.Progress, database_name='Pax', amount_loaded=1)
            self.get_or_create_object(models.Progress, database_name='Sabio', amount_loaded=1)
            self.get_or_create_object(models.Progress, database_name='Array Express', amount_loaded=1)
            self.get_or_create_object(models.Progress, database_name='IntAct', amount_loaded=0)
            self.session.commit()
            self.load_small_db_switch = True
            self.load_content()

    def load_content(self):
        """
        A wrapper for loading all the databases into common ORM database

        """
        self.vprint(''' \n
            ===================================
            |                                 |
            |                                 |
            |    Starting Datanator Build     |
            |                                 |
            |                                 |
            ===================================

            ''')
        t0 = time.time()

        observation = models.Observation()
        observation.physical_entity = models.PhysicalEntity()
        self.entity = observation.physical_entity
        observation.physical_property = models.PhysicalProperty()
        self.property = observation.physical_property

        # Chunk Larger DBs
        self.build_pax()
        self.build_intact_interactions()
        self.build_sabio()
        self.build_array_express()

        # Add complete smaller DBs
        if self.load_small_db_switch:
            self.build_intact_complexes()
            self.build_corum()
            self.build_jaspar()
            self.build_ecmdb()

        # Add missing subunit information
        self.build_uniprot()

        # Add missing Taxon information
        self.build_ncbi()


        self.vprint(''' \n
            =============================================
            |                                           |
            |             Finished Build                |
            |   Total time taken for build: %s secs     |
            |                                           |
            =============================================
            ''' % (str(round(time.time() - t0, 1))))

    @timeit
    def build_intact_interactions(self):
        """
        Collects IntAct.sqlite file and integrates interaction data into the common ORM

        """

        t0 = time.time()
        intactdb = intact.IntAct(cache_dirname=self.cache_dirname)
        batch = INTACT_INTERACTION_TEST_BATCH if self.test else INTACT_INTERACTION_BUILD_BATCH
        intact_progress = self.session.query(
            models.Progress).filter_by(database_name='IntAct').first()
        load_count = intact_progress.amount_loaded

        if self.max_entries == float('inf'):
            interactiondb = intactdb.session.query(
                intact.ProteinInteraction).all()
        else:
            interactiondb = intactdb.session.query(intact.ProteinInteraction).filter(intact.ProteinInteraction.index.in_
                (range(load_count, load_count + batch))).all()


        batch_list = []
        for i in interactiondb:

            metadata = self.get_or_create_object(models.Metadata, name='protein_interaction_' + str(i.index))
            metadata.method.append(self.get_or_create_object(models.Method, name=i.method))
            metadata.resource.append(self.get_or_create_object(models.Resource, namespace='pubmed', _id=i.publication))
            metadata.resource.append(self.get_or_create_object(models.Resource, namespace='paper', _id=i.publication_author))

            for c in dir(i):
                if getattr(i, c) == None and '__' not in c:
                    setattr(i, c, '')

            for type, protein, gene in [(i.type_a, i.protein_a, i.gene_a), (i.type_b, i.protein_b, i.gene_b)]:
                if type == 'protein':
                    self.get_or_create_object(models.ProteinSubunit, uniprot_id=protein, gene_name=gene)

            batch_list.append(models.ProteinInteraction(name=i.protein_a + " + " + i.protein_b, type='protein protein interaction', protein_a=i.protein_a,
                protein_b=i.protein_b,gene_a=i.gene_a, gene_b=i.gene_b, loc_a=i.feature_a, loc_b=i.feature_b, type_a=i.type_a, type_b=i.type_b,
                stoich_a=i.stoich_a, stoich_b=i.stoich_b, confidence=i.confidence,interaction_type=i.interaction_type,
                role_a = i.role_a, role_b= i.role_b, _metadata=metadata))

        self.session.add_all(batch_list)

        intact_progress.amount_loaded = load_count + batch

        self.vprint('Comitting')
        self.session.commit()


    @timeit
    def build_pax(self):
        """
        Collects Pax.sqlite file and integrates its data into the common ORM

        Total datasets: 493

        """

        paxdb = pax.Pax(cache_dirname=self.cache_dirname, verbose=self.verbose)
        pax_ses = paxdb.session
        batch = PAX_TEST_BATCH if self.test else PAX_BUILD_BATCH

        pax_progress = self.session.query(models.Progress).filter_by(database_name='Pax').first()

        load_count = pax_progress.amount_loaded

        if self.max_entries == float('inf'):
            pax_dataset = pax_ses.query(pax.Dataset).all()
        else:
            pax_dataset = pax_ses.query(pax.Dataset).filter(pax.Dataset.id.in_
                                                            (range(load_count, load_count + batch))).all()

        for dataset in pax_dataset:
            metadata = self.get_or_create_object(
                models.Metadata, name=dataset.file_name)
            metadata.taxon.append(self.get_or_create_object(
                models.Taxon, ncbi_id=dataset.taxon_ncbi_id))
            metadata.resource.append(self.get_or_create_object(
                models.Resource, namespace='url', _id=dataset.publication))
            self.property.abundance_dataset = self.get_or_create_object(models.AbundanceDataSet, type='Protein Abundance Dataset',
                                                                        name=dataset.file_name, file_name=dataset.file_name, score=dataset.score, weight=dataset.weight,
                                                                        coverage=dataset.coverage, _metadata=metadata)
            abundance = pax_ses.query(pax.Observation).filter_by(
                dataset_id=dataset.id).all()
            uni = [str(d.protein.uniprot_id) for d in abundance]

            self.session.bulk_insert_mappings(models.AbundanceData,
                                              [
                                                  dict(abundance=data.abundance, pax_load=data.dataset_id,
                                                       uniprot_id=data.protein.uniprot_id)
                                                  for data in abundance
                                              ])

            self.session.bulk_insert_mappings(models.ProteinSubunit,
                                              [
                                                  dict(uniprot_id=data.protein.uniprot_id, type='Protein Subunit',
                                                       pax_load=data.dataset_id) for data in abundance
                                              ])

            self.session.commit()

            for subunit in self.session.query(models.ProteinSubunit).filter_by(pax_load=dataset.id).all():
                subunit._metadata = metadata

            for rows in self.session.query(models.AbundanceData).filter_by(pax_load=dataset.id).all():
                rows.subunit = self.session.query(models.ProteinSubunit).filter_by(
                    pax_load=dataset.id).filter_by(uniprot_id=rows.uniprot_id).first()
                rows.dataset = self.property.abundance_dataset

        pax_progress.amount_loaded = load_count + batch

        self.vprint('Comitting')
        self.session.commit()

    @timeit
    def build_array_express(self):

        ae = array_express.ArrayExpress(cache_dirname=self.cache_dirname)
        batch = ARRAY_EXPRESS_TEST_BATCH if self.test else ARRAY_EXPRESS_BUILD_BATCH
        array_progress = self.session.query(models.Progress).filter_by(
            database_name='Array Express').first()
        load_count = array_progress.amount_loaded

        if self.max_entries == float('inf'):
            experiments = ae.session.query(array_express.Experiment).all()
        else:
            experiments = ae.session.query(array_express.Experiment).filter(array_express.Experiment._id.in_
                                                                            (range(load_count, load_count + batch)))

        for exp in experiments:

            exp_metadata = self.get_or_create_object(
                models.ExperimentMetadata,
                name="RNA-Seq Experiment Info: " + exp.id,
                description=exp.description
            )

            exp_metadata.taxon = [
                self.get_or_create_object(
                    models.Taxon,
                    name=org.name,
                    ncbi_id=org.ncbi_id
                )
                for org in exp.organisms]

            exp_metadata.resource.append(self.get_or_create_object(
                models.Resource,
                namespace="ArrayExpress",
                _id=exp.id,
                release_date=str(exp.release_date)
            ))

            exp_metadata.method = [
                self.get_or_create_object(
                    models.Method,
                    name=protocol.protocol_type,
                    comments=protocol.text,
                    performer=protocol.performer,
                    hardware=protocol.hardware,
                    software=protocol.software,
                )
                for protocol in exp.protocols]

            exp_metadata.experiment_design = [
                self.get_or_create_object(
                    models.ExperimentDesign,
                    name=exp_des.name,
                )
                for exp_des in exp.designs]

            exp_metadata.experiment_type = [
                self.get_or_create_object(
                    models.ExperimentType,
                    name=exp_type.name,
                )
                for exp_type in exp.types]

            exp_metadata.data_format = [
                self.get_or_create_object(
                    models.DataFormat,
                    name=data_format.name,
                    bio_assay_data_cubes=data_format.bio_assay_data_cubes,
                )
                for data_format in exp.data_formats]

            flask_experiment = self.get_or_create_object(
                models.RNASeqExperiment,
                #name = "exp.id",
                #type = "Array Express Experiment",
                exp_name=exp.name,
                accession_number=exp.id,
                _experimentmetadata=exp_metadata
            )
            for sample in ae.session.query(array_express.Sample).filter_by(experiment_id=exp.id).all():
                # m_name = "RNA-Seq Sample Info: " + sample.experiment_id + "_" + sample.name #{}_{}".format(sample.experiment_id, s_name)
                metadata = self.get_or_create_object(
                    models.Metadata,
                    name="RNA-Seq Sample Info: " + sample.experiment_id + "_" + sample.name
                )

                metadata.characteristic = [
                    self.get_or_create_object(models.Characteristic,
                                              category=characteristic.category,
                                              value=characteristic.value)
                    for characteristic in sample.characteristics
                ]

                metadata.variable = [
                    self.get_or_create_object(models.Variable,
                                              category=variable.name,
                                              value=variable.value,
                                              units=variable.unit)
                    for variable in sample.variables
                ]

                flask_sample = self.get_or_create_object(
                    models.RNASeqDataSet,
                    name=sample.experiment_id + "_" + sample.name,
                    type="Array Express Sample",
                    experiment_accession_number=sample.experiment_id,
                    sample_name=sample.name,
                    assay=sample.assay,
                    ensembl_organism_strain=sample.ensembl_organism_strain,
                    read_type=sample.read_type,
                    full_strain_specificity=sample.full_strain_specificity,
                    _metadata=metadata
                )

                flask_sample.reference_genome = [
                    self.get_or_create_object(models.ReferenceGenome,
                                              namespace="Ensembl",
                                              organism_strain=ensembl_info.organism_strain,
                                              download_url=ensembl_info.url,
                                              )
                    for ensembl_info in sample.ensembl_info]

                flask_experiment.samples.append(flask_sample)

        array_progress.amount_loaded = load_count + batch


        self.vprint('Comitting')
        self.session.commit()

    @timeit
    def build_sabio(self):
        """
        Collects SabioRK.sqlite file and integrates its data into the common ORM

        """

        t0 = time.time()
        sabiodb = sabio_rk.SabioRk(
            cache_dirname=self.cache_dirname, verbose=self.verbose)
        sabio_ses = sabiodb.session
        batch = SABIO_TEST_BATCH if self.test else SABIO_BUILD_BATCH
        sabio_progress = self.session.query(
            models.Progress).filter_by(database_name='Sabio').first()
        load_count = sabio_progress.amount_loaded

        if self.max_entries == float('inf'):
            sabio_entry = sabio_ses.query(sabio_rk.Entry)
        else:
            sabio_entry = sabio_ses.query(sabio_rk.Entry).filter(sabio_rk.Entry._id.in_
                                                                 (range(load_count, load_count + batch)))

        for item in sabio_entry:
            metadata = self.get_or_create_object(models.Metadata, name='Kinetic Law ' + str(
                item.id)) if item._type == 'kinetic_law' else self.get_or_create_object(models.Metadata, name=item.name)
            metadata.synonym = [self.get_or_create_object(
                models.Synonym, name=synonyms.name) for synonyms in item.synonyms]
            metadata.taxon = [self.get_or_create_object(
                models.Taxon, ncbi_id=docs.id) for docs in item.cross_references if docs.namespace == 'taxonomy']
            uniprot = [
                docs.id for docs in item.cross_references if docs.namespace == 'uniprot']
            metadata.resource = [self.get_or_create_object(
                models.Resource, namespace=docs.namespace, _id=docs.id) for docs in item.cross_references]
            compartment = self.get_or_create_object(
                models.CellCompartment, name=item.name) if item._type == 'compartment' else None
            if compartment:
                continue

            if item._type == 'compound':
                structure = None
                for struct in item.structures:
                    structure = self.get_or_create_object(models.Structure, type='Structure', name=item.name,
                                                          _value_smiles=struct.value, _value_inchi=struct._value_inchi,
                                                          _structure_formula_connectivity=struct._value_inchi_formula_connectivity, _metadata=metadata)\
                        if struct.format == 'smiles' else None

                self.entity.compound = self.get_or_create_object(models.Compound, type='Compound',
                                                                 name=item.name, compound_name=item.name,
                                                                 _is_name_ambiguous=sabio_ses.query(
                                                                     sabio_rk.Compound).get(item._id)._is_name_ambiguous,
                                                                 structure=structure, _metadata=metadata)
                continue

            elif item._type == 'enzyme':
                complx = sabio_ses.query(sabio_rk.Enzyme).get(item._id)
                self.entity.protein_complex = self.get_or_create_object(models.ProteinComplex, type='Enzyme',
                                                                        name=item.name, complex_name=item.name,
                                                                        molecular_weight=item.molecular_weight, funcat_dsc='Enzyme', _metadata=metadata)
                continue

            elif item._type == 'enzyme_subunit':
                result = self.session.query(models.ProteinComplex).filter_by(
                    complex_name=item.enzyme.name).first()
                self.entity.protein_subunit = self.get_or_create_object(models.ProteinSubunit, type='Enzyme Subunit',
                                                                        name=item.name, subunit_name=item.name,
                                                                        uniprot_id=uniprot, coefficient=item.coefficient,
                                                                        molecular_weight=item.molecular_weight,
                                                                        proteincomplex=result, _metadata=metadata)
                continue

            elif item._type == 'kinetic_law':

                catalyst = self.session.query(models.ProteinComplex).filter_by(
                    complex_name=item.enzyme.name).first() if item.enzyme_id else None
                metadata.resource.extend([self.get_or_create_object(
                    models.Resource, namespace=resource.namespace, _id=resource.id) for resource in item.references])
                metadata.taxon.append(self.get_or_create_object(
                    models.Taxon, ncbi_id=item.taxon))
                metadata.cell_line.append(self.get_or_create_object(
                    models.CellLine, name=item.taxon_variant))
                metadata.conditions.append(self.get_or_create_object(
                    models.Conditions, temperature=item.temperature, ph=item.ph, media=item.media))
                self.property.kinetic_law = self.get_or_create_object(models.KineticLaw, type='Kinetic Law', enzyme=catalyst,
                                                                      enzyme_type=item.enzyme_type, tissue=item.tissue, mechanism=item.mechanism, equation=item.equation, _metadata=metadata)

                def common_schema_compound(sabio_object):
                    compound_name = sabio_object.name
                    return self.session.query(models.Compound).filter_by(compound_name=compound_name).first()

                def common_schema_compartment(sabio_object):
                    if sabio_object:
                        compartment_name = sabio_object.name
                        return self.session.query(models.CellCompartment).filter_by(name=compartment_name).first()
                    else:
                        return None

                reactants = [models.Reaction(compound=common_schema_compound(r.compound),
                                             compartment=common_schema_compartment(r.compartment), _is_reactant=1, rxn_type=r.type,
                                             kinetic_law=self.property.kinetic_law) for r in item.reactants if item.reactants]

                products = [models.Reaction(compound=common_schema_compound(p.compound),
                                            compartment=common_schema_compartment(p.compartment), _is_product=1, rxn_type=p.type,
                                            kinetic_law=self.property.kinetic_law) for p in item.products if item.products]

                modifier = [models.Reaction(compound=common_schema_compound(m.compound),
                                            compartment=common_schema_compartment(m.compartment), _is_modifier=1, rxn_type=m.type,
                                            kinetic_law=self.property.kinetic_law) for m in item.modifiers if item.products]

                for param in item.parameters:
                    parameter = models.Parameter(sabio_type=param.type, value=param.value, error=param.error,
                                                 units=param.units, observed_name=param.observed_name, kinetic_law=self.property.kinetic_law,
                                                 observed_sabio_type=param.observed_type, observed_value=param.observed_value,
                                                 observed_error=param.observed_error, observed_units=param.observed_units)
                    parameter.compound = common_schema_compound(
                        param.compound) if param.compound else None
                continue

        sabio_progress.amount_loaded = load_count + batch

        self.vprint('Comitting')
        self.session.commit()


    @timeit
    def build_corum(self):
        """
        Collects Corum.sqlite file and integrates its data into the common ORM

        """

        corumdb = corum.Corum(
            cache_dirname=self.cache_dirname, verbose=self.verbose)
        corum_ses = corumdb.session

        corum_complex = corum_ses.query(corum.Complex).all()
        corum_subunit = corum_ses.query(corum.Subunit).all()

        max_entries = self.max_entries

        if self.load_entire_small_DBs:
            max_entries = float('inf')

        entries = 0
        for complx in corum_complex:
            if entries < max_entries:
                entry = complx.observation
                metadata = self.get_or_create_object(
                    models.Metadata, name=complx.complex_name)
                metadata.taxon.append(self.get_or_create_object(
                    models.Taxon, ncbi_id=entry.taxon_ncbi_id))
                metadata.resource.append(self.get_or_create_object(
                    models.Resource, namespace='pubmed', _id=str(entry.pubmed_id)))
                metadata.method.append(self.get_or_create_object(
                    models.Method, name='purification', comments=entry.pur_method))
                metadata.cell_line.append(self.get_or_create_object(
                    models.CellLine, name=entry.cell_line))
                self.entity.protein_complex = self.get_or_create_object(models.ProteinComplex,
                                                                        type='Protein Complex', name=complx.complex_name, complex_name=complx.complex_name, go_id=complx.go_id,
                                                                        go_dsc=complx.go_dsc, funcat_id=complx.funcat_id, funcat_dsc=complx.funcat_dsc, su_cmt=complx.su_cmt,
                                                                        complex_cmt=complx.complex_cmt, disease_cmt=complx.disease_cmt,  _metadata=metadata)
                entries += 1

        entries = 0
        for subunit in corum_subunit:
            if entries < max_entries:
                complx = self.session.query(models.ProteinComplex).filter_by(
                    complex_name=subunit.complex.complex_name).first()
                entry = self.session.query(
                    models.Observation).get(complx.complex_id)
                self.entity.protein_subunit = self.get_or_create_object(models.ProteinSubunit,
                                                                        type='Protein Subunit', uniprot_id=subunit.su_uniprot,
                                                                        entrez_id=subunit.su_entrezs, name=subunit.protein_name, subunit_name=subunit.protein_name, gene_name=subunit.gene_name,
                                                                        gene_syn=subunit.gene_syn, proteincomplex=complx, _metadata=self.session.query(models.Metadata).get(entry._metadata_id))
                entries += 1


        self.vprint('Comitting')
        self.session.commit()

    @timeit
    def build_jaspar(self):
        """
        Collects Jaspar.sqlite file and integrates its data into the common ORM

        Total datasets: 2404

        """

        jaspardb = jaspar.Jaspar(
            cache_dirname=self.cache_dirname, verbose=self.verbose)

        jasp_ses = jaspardb.session

        self.entity = self.entity
        self.property = self.property

        def list_to_string(list_):
            if list_:
                ans = ''
                for word in list_:
                    ans += word
                return ans
            else:
                return ''

        matrix = jasp_ses.query(jaspar.Matrix).all()

        max_entries = self.max_entries

        if self.load_entire_small_DBs:
            max_entries = float('inf')

        entries = 0
        for entry in matrix:
            if entries <= max_entries:
                annotations = jasp_ses.query(
                    jaspar.Annotation).filter_by(ID=entry.ID)
                class_ = [c.VAL for c in annotations.filter_by(
                    TAG='class').all()]
                family_ = [f.VAL for f in annotations.filter_by(
                    TAG='family').all()]
                pubmed = [p.VAL for p in annotations.filter_by(
                    TAG='medline').all()]
                type_ = [t.VAL for t in annotations.filter_by(
                    TAG='type').all()]
                species = [s.TAX_ID for s in jasp_ses.query(
                    jaspar.Species).filter_by(ID=entry.ID).all()]
                protein = [p.ACC for p in jasp_ses.query(
                    jaspar.Protein).filter_by(ID=entry.ID).all()]

                metadata = self.get_or_create_object(
                    models.Metadata, name=entry.NAME + ' Binding Motif')
                metadata.method.append(self.get_or_create_object(
                    models.Method, name=list_to_string(type_)))
                metadata.resource = [self.get_or_create_object(
                    models.Resource, namespace='pubmed', _id=ref) for ref in pubmed]
                metadata.taxon = [self.get_or_create_object(
                    models.Taxon, ncbi_id=int(tax)) for tax in species if tax != '-']
                if '::' in entry.NAME:
                    self.entity.protein_complex = self.get_or_create_object(models.ProteinComplex, type='Transcription Factor Complex',
                                                                            name=entry.NAME, complex_name=entry.NAME, complex_cmt='transcription factor', class_name=list_to_string(class_),
                                                                            family_name=list_to_string(family_), _metadata=metadata)
                    self.property.dna_binding_dataset = self.get_or_create_object(models.DNABindingDataset, type='DNA Binding Dataset',
                                                                                  name=entry.NAME, version=entry.VERSION, tf=self.entity.protein_complex, _metadata=metadata)
                else:
                    prot = protein[0] if protein else None
                    self.entity.protein_subunit = self.get_or_create_object(models.ProteinSubunit, uniprot_id=prot,
                                                                            type='Transcription Factor Subunit', name=entry.NAME, subunit_name=entry.NAME, gene_name=entry.NAME,
                                                                            class_name=list_to_string(class_), family_name=list_to_string(family_), _metadata=metadata)
                    self.property.dna_binding_dataset = self.get_or_create_object(models.DNABindingDataset, type='DNA Binding Dataset',
                                                                                  name=entry.NAME, version=entry.VERSION, subunit=self.entity.protein_subunit, _metadata=metadata)
                subquery = jasp_ses.query(jaspar.Data).filter_by(ID=entry.ID)
                for position in range(1, 1 + max(set([c.col for c in subquery.all()]))):
                    freq = subquery.filter_by(col=position)
                    self.get_or_create_object(models.DNABindingData, position=position, frequency_a=freq.filter_by(row='A').first().val,
                                              frequency_c=freq.filter_by(row='C').first(
                    ).val, frequency_t=freq.filter_by(row='T').first().val,
                        frequency_g=freq.filter_by(row='G').first().val, jaspar_id=entry.ID, dataset=self.property.dna_binding_dataset)
            entries += 1


        self.vprint('Comitting')
        self.session.commit()

    @timeit
    def build_ecmdb(self):
        """
        Collects ECMDB.sqlite file and integrates its data into the common ORM

        """

        ecmDB = ecmdb.Ecmdb(
            cache_dirname=self.cache_dirname, verbose=self.verbose)
        ecm_ses = ecmDB.session
        ecmdb_compound = ecm_ses.query(ecmdb.Compound).all()

        max_entries = self.max_entries

        if self.load_entire_small_DBs:
            max_entries = float('inf')

        entries = 0
        for compound in ecmdb_compound:
            if entries < max_entries:
                concentration = compound.concentrations
                ref = compound.cross_references
                compart = compound.compartments
                syn = compound.synonyms
                metadata = self.get_or_create_object(
                    models.Metadata, name=compound.name)
                tax = self.session.query(
                    models.Taxon).filter_by(ncbi_id=562).first()
                metadata.taxon.append(tax) if tax else metadata.taxon.append(
                    self.get_or_create_object(models.Taxon, ncbi_id=562))
                metadata.resource = [self.get_or_create_object(
                    models.Resource, namespace=docs.namespace, _id=docs.id) for docs in ref]
                metadata.cell_compartment = [self.get_or_create_object(
                    models.CellCompartment, name=areas.name) for areas in compart]
                metadata.synonym = [self.get_or_create_object(
                    models.Synonym, name=syns.name) for syns in syn]
                self.property.structure = self.get_or_create_object(models.Structure, type='Structure', name=compound.name,
                                                                    _value_inchi=compound.structure,
                                                                    _structure_formula_connectivity=compound._structure_formula_connectivity, _metadata=metadata)
                self.entity.compound = self.get_or_create_object(models.Compound, type='Compound', name=compound.name,
                                                                 compound_name=compound.name, description=compound.description, comment=compound.comment, structure=self.property.structure,
                                                                 _metadata=metadata)
                index = 0 if concentration else -1
                if index == -1:
                    continue
                for rows in concentration:
                    new_metadata = self.get_or_create_object(
                        models.Metadata, name=compound.name + ' Concentration ' + str(index))
                    new_metadata.taxon = metadata.taxon
                    new_metadata.cell_compartment = metadata.cell_compartment
                    new_metadata.synonym = metadata.synonym
                    new_metadata.cell_line.append(
                        self.get_or_create_object(models.CellLine, name=rows.strain))
                    new_metadata.conditions.append(self.get_or_create_object(models.Conditions, growth_status=rows.growth_status,
                                                                             media=rows.media, temperature=rows.temperature, growth_system=rows.growth_system))
                    new_metadata.resource = [self.get_or_create_object(
                        models.Resource, namespace=docs.namespace, _id=docs.id) for docs in rows.references]
                    self.property.concentration = self.get_or_create_object(models.Concentration, type='Concentration', name=compound.name + ' Concentration ' + str(index),
                                                                            value=rows.value, error=rows.error, units='uM', _metadata=new_metadata, compound=self.entity.compound)
                    index += 1
                entries += 1


        self.vprint('Comitting')
        self.session.commit()


    @timeit
    def build_intact_complexes(self):
        """
        Collects IntAct.sqlite file and integrates complex data into the common ORM

        """
        intactdb = intact.IntAct(cache_dirname=self.cache_dirname)
        intact_ses = intactdb.session

        def parse_string(regex, string):
            result = re.findall(regex, string)
            return '|'.join(result)

        max_entries = self.max_entries

        if self.load_entire_small_DBs:
            max_entries = float('inf')

        complexdb = intactdb.session.query(intact.ProteinComplex).all() if max_entries == float('inf') \
            else intactdb.session.query(intact.ProteinComplex).limit(self.max_entries).all()

        for row in complexdb:
            metadata = self.get_or_create_object(
                models.Metadata, name=row.name)
            metadata.taxon.append(self.get_or_create_object(
                models.Taxon, ncbi_id=row.ncbi))
            metadata.resource.append(self.get_or_create_object(
                models.Resource, namespace='ebi id', _id=row.identifier))

            go_dsc = parse_string(re.compile(".*?\((.*?)\)"), row.go_annot)
            go_id = parse_string(re.compile("GO:(.*?)\("), row.go_annot)
            subunits = re.findall(re.compile(
                ".*?\|(.*?)\("), '|' + row.subunits)
            self.entity.protein_complex = self.get_or_create_object(models.ProteinComplex, name=row.name, type='Protein Comlplex',
                                                                    complex_name=row.name, su_cmt=row.subunits, go_dsc=go_dsc, go_id=go_id,
                                                                    complex_cmt=row.desc, _metadata=metadata)
            for sub in subunits:
                if 'CHEBI' not in sub:
                    self.entity.protein_subunit = self.get_or_create_object(models.ProteinSubunit,
                                                                            type='Protein Subunit', uniprot_id=sub,
                                                                            proteincomplex=self.entity.protein_complex, _metadata=metadata)

        self.vprint('Comitting')
        self.session.commit()


    @timeit
    def build_uniprot(self):
        """
        Collects Uniprot.sqlite file and integrates data into existing ProteinSubunit table

        """

        unidb = uniprot.Uniprot(cache_dirname=self.cache_dirname)
        unidb_ses = unidb.session

        com_unis = self.session.query(models.ProteinSubunit).filter_by(
            uniprot_checked=None).all()

        for subunit in com_unis:
            info = unidb_ses.query(uniprot.UniprotData).filter_by(
                uniprot_id=subunit.uniprot_id).first()
            subunit.uniprot_checked = True
            if info:
                subunit.subunit_name = subunit.name = info.entry_name if not subunit.subunit_name else subunit.subunit_name
                subunit.entrez_id = info.entrez_id if not subunit.entrez_id else subunit.entrez_id
                subunit.gene_name = info.gene_name if not subunit.gene_name else subunit.gene_name
                subunit.canonical_sequence = info.canonical_sequence if not subunit.canonical_sequence else subunit.canonical_sequence
                subunit.length = info.length if not subunit.length else subunit.length
                subunit.mass = info.mass if not subunit.mass else subunit.mass


        self.vprint('Comitting')
        self.session.commit()

    def build_ncbi(self):
        """
        Uses NCBI package to integrate data into existing Taxon table

        """

        ncbi = NCBITaxa()

        ncbi_ids = []
        species = self.session.query(models.Taxon).all()
        for items in species:
            ncbi_ids.append(items.ncbi_id)

        species_dict = ncbi.get_taxid_translator(ncbi_ids)

        for tax in species:
            if tax.ncbi_id in species_dict.keys():
                tax.name = species_dict[tax.ncbi_id]


        self.vprint('Comitting')
        self.session.commit()
