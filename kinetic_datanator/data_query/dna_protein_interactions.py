"""
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-09-28
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model, data_query,  common_schema
from kinetic_datanator.app import models, flask_common_schema
from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import tempfile
import shutil
import csv


class ProteintoDNAInteractionQueryGenerator(data_query.CachedDataSourceQueryGenerator):
    """ Queries Proteins for their DNA """

    def __init__(self,
                 taxon=None, max_taxon_dist=None, taxon_dist_scale=None, include_variants=False,
                 temperature=37., temperature_std=1.,
                 ph=7.5, ph_std=0.3):
        """
        Args:
            taxon (:obj:`str`, optional): target taxon
            max_taxon_dist (:obj:`int`, optional): maximum taxonomic distance to include
            taxon_dist_scale (:obj:`float`, optional): The scale of the taxonomic distance scoring distribution.
                This determines how quickly the score falls to zero away from zero.
            include_variants (:obj:`bool`, optional): if :obj:`True`, also include observations from mutant taxa
            temperature (:obj:`float`, optional): desired temperature to search for
            temperature_std (:obj:`float`, optional): how much to penalize observations from other temperatures
            ph (:obj:`float`, optional): desired pH to search for
            ph_std (:obj:`float`, optional): how much to penalize observations from other pHs
        """
        super(ProteintoDNAInteractionQueryGenerator, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=common_schema.CommonSchema())

    def get_observed_values(self, protein):
        """ Find the DNA binding motif for a given protein

        Args:
            protein (:obj:`data_model.ProteinSpecie`): Protein to find data for

        Returns:
            :obj:`list` of :obj:`data_model.Observable`: list of observables

        """
        versions = self.get_DNA_by_protein(protein)

        index = 0
        observable = []
        for motif in versions:
            binding_matrix = []
            for position in motif.all():
                binding_matrix.append([position.frequency_a, position.frequency_c, position.frequency_g, position.frequency_t])
            binding_matrix = map(list, zip(*binding_matrix))
            self.cache_dirname = tempfile.mkdtemp()
            with open(self.cache_dirname+'/data.pfm', 'w') as pfm:
                for items in binding_matrix:
                    writer = csv.writer(pfm, delimiter = '\t')
                    writer.writerow(items)

            m = motifs.read(open(self.cache_dirname+'/data.pfm'), 'pfm')

            dna_specie = data_model.DnaSpecie(binding_matrix = m.counts, sequence = str(m.counts.consensus))
            interaction = data_model.Interaction(name = 'Transcription Factor DNA Binding Site')
            observable.append(data_model.Observable(specie = dna_specie, interaction = interaction))

            for position in motif.all():
                observable[index].specie.cross_references = data_model.Resource(namespace ='pubmed',\
                id = position.dataset._metadata.resource[0]._id),
                break

            shutil.rmtree(self.cache_dirname)
            index += 1

        return observable


    def get_DNA_by_protein(self, protein, select = common_schema.DNABindingData):
        """
        Args:
            specie (:obj:`data_model.ProteinSpecie`): species to find data for

        Returns:
            :obj:`list` of :obj:`sqlalchemy.orm.query.Query`: list of versions queries for DNA binding observed for similar proteins
        """

        q = []

        condition = common_schema.ProteinSubunit.uniprot_id == protein.uniprot_id

        versions = self.data_source.session.query(common_schema.DNABindingDataset).\
            join((common_schema.ProteinSubunit, common_schema.DNABindingDataset.subunit)).filter(condition).all()

        for version in range(1,len(versions)+1):
            q.append(self.data_source.session.query(select).join((common_schema.DNABindingDataset, select.dataset))\
                .filter(common_schema.DNABindingDataset.version == version)\
                .join((common_schema.ProteinSubunit, common_schema.DNABindingDataset.subunit)).filter(condition))

        return q



class DNAtoProteinInteractionQueryGenerator(data_query.CachedDataSourceQueryGenerator):
    """ Queries DNA reigons for relevant proteins """

    def __init__(self,
                 taxon=None, max_taxon_dist=None, taxon_dist_scale=None, include_variants=False,
                 temperature=37., temperature_std=1.,
                 ph=7.5, ph_std=0.3):
        """
        Args:
            taxon (:obj:`str`, optional): target taxon
            max_taxon_dist (:obj:`int`, optional): maximum taxonomic distance to include
            taxon_dist_scale (:obj:`float`, optional): The scale of the taxonomic distance scoring distribution.
                This determines how quickly the score falls to zero away from zero.
            include_variants (:obj:`bool`, optional): if :obj:`True`, also include observations from mutant taxa
            temperature (:obj:`float`, optional): desired temperature to search for
            temperature_std (:obj:`float`, optional): how much to penalize observations from other temperatures
            ph (:obj:`float`, optional): desired pH to search for
            ph_std (:obj:`float`, optional): how much to penalize observations from other pHs
        """
        super(DNAtoProteinInteractionQueryGenerator, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=common_schema.CommonSchema())


    def get_observed_values(self, DnaSpecie):
        """ Find observed concentrations for the metabolite or similar metabolites

        Args:
            specie (:obj:`data_model.DnaSpecie`): species to find data for

        Returns:
            :obj:`list` of :obj:`data_model.Observable`: list of observable
        """
        proteins = self.get_protein_by_DNA_sequence(DnaSpecie.sequence)

        observable = []
        for subunits in proteins:

            protein = data_model.ProteinSpecie(uniprot_id = subunits[0].uniprot_id,\
                entrez_id = subunits[0].entrez_id, gene_name = subunits[0].gene_name,\
                length = subunits[0].length, mass = subunits[0].mass)
            protein.cross_references = []
            for doc in subunits[0]._metadata.resource:
                protein.cross_references.append(data_model.Resource(namespace = 'pubmed', id = doc._id))

            interaction = data_model.Interaction(name = 'Transcription Factor DNA Binding Site', \
                position = subunits[1], score = subunits[2])

            observable.append(data_model.Observable(specie = protein, interaction = interaction))

        return observable


    def get_protein_by_DNA_sequence(self,sequence, select = common_schema.ProteinSubunit):
        """
        Args:
            sequence (:obj:`data_model.DnaSpecie.sequence`): sequence of DNA segment

        Returns:
            :obj:`list` of :obj:`tuple`: Returns the query for a protein, sequence position, and score

        """
        #TODO: Make more efficient. Add gene location filter

        ans = []

        size = len(sequence)
        all_matricies = self.data_source.session.query(common_schema.DNABindingDataset).all()
        all_sequences = []
        for matricies in all_matricies:
            if len(matricies.dna_binding_data) == size:
                binding_matrix = []
                for position in matricies.dna_binding_data:
                    binding_matrix.append([position.frequency_a, position.frequency_c, position.frequency_g, position.frequency_t])
                binding_matrix = map(list, zip(*binding_matrix))
                self.cache_dirname = tempfile.mkdtemp()
                with open(self.cache_dirname+'/data.pfm', 'w') as pfm:
                    for items in binding_matrix:
                        writer = csv.writer(pfm, delimiter = '\t')
                        writer.writerow(items)
                m = motifs.read(open(self.cache_dirname+'/data.pfm'), 'pfm')
                my_seq = Seq(sequence, IUPAC.unambiguous_dna)
                ##TODO: Include selective threshold
                # distribution = m.pssm.distribution().threshold_paster()
                # print distribution.threshold_paster()
                for position, score in m.pssm.search(my_seq, threshold = 2):
                    ans.append((matricies.subunit, position, score))

        return ans


class FlaskProteintoDNAInteractionQueryGenerator(data_query.CachedDataSourceQueryGenerator):
    """ Queries Proteins for their DNA """

    def __init__(self,
                 taxon=None, max_taxon_dist=None, taxon_dist_scale=None, include_variants=False,
                 temperature=37., temperature_std=1.,
                 ph=7.5, ph_std=0.3):
        """
        Args:
            taxon (:obj:`str`, optional): target taxon
            max_taxon_dist (:obj:`int`, optional): maximum taxonomic distance to include
            taxon_dist_scale (:obj:`float`, optional): The scale of the taxonomic distance scoring distribution.
                This determines how quickly the score falls to zero away from zero.
            include_variants (:obj:`bool`, optional): if :obj:`True`, also include observations from mutant taxa
            temperature (:obj:`float`, optional): desired temperature to search for
            temperature_std (:obj:`float`, optional): how much to penalize observations from other temperatures
            ph (:obj:`float`, optional): desired pH to search for
            ph_std (:obj:`float`, optional): how much to penalize observations from other pHs
        """
        super(FlaskProteintoDNAInteractionQueryGenerator, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=flask_common_schema.FlaskCommonSchema())

    def get_observed_values(self, protein):
        """ Find the DNA binding motif for a given protein

        Args:
            protein (:obj:`data_model.ProteinSpecie`): Protein to find data for

        Returns:
            :obj:`list` of :obj:`data_model.Observable`: list of observables

        """
        versions = self.get_DNA_by_protein(protein)

        index = 0
        observable = []
        for motif in versions:
            binding_matrix = []
            for position in motif.all():
                binding_matrix.append([position.frequency_a, position.frequency_c, position.frequency_g, position.frequency_t])
            binding_matrix = map(list, zip(*binding_matrix))
            self.cache_dirname = tempfile.mkdtemp()
            with open(self.cache_dirname+'/data.pfm', 'w') as pfm:
                for items in binding_matrix:
                    writer = csv.writer(pfm, delimiter = '\t')
                    writer.writerow(items)

            m = motifs.read(open(self.cache_dirname+'/data.pfm'), 'pfm')

            dna_specie = data_model.DnaSpecie(binding_matrix = m.counts, sequence = str(m.counts.consensus))
            interaction = data_model.Interaction(name = 'Transcription Factor DNA Binding Site')
            observable.append(data_model.Observable(specie = dna_specie, interaction = interaction))

            for position in motif.all():
                observable[index].specie.cross_references = data_model.Resource(namespace ='pubmed',\
                id = position.dataset._metadata.resource[0]._id),
                break

            shutil.rmtree(self.cache_dirname)
            index += 1

        return observable


    def get_DNA_by_protein(self, protein, select = models.DNABindingData):
        """
        Args:
            specie (:obj:`data_model.ProteinSpecie`): species to find data for

        Returns:
            :obj:`list` of :obj:`sqlalchemy.orm.query.Query`: list of versions queries for DNA binding observed for similar proteins
        """

        q = []

        condition = models.ProteinSubunit.uniprot_id == protein.uniprot_id

        versions = self.data_source.session.query(models.DNABindingDataset).\
            join((models.ProteinSubunit, models.DNABindingDataset.subunit)).filter(condition).all()

        for version in range(1,len(versions)+1):
            q.append(self.data_source.session.query(select).join((models.DNABindingDataset, select.dataset))\
                .filter(models.DNABindingDataset.version == version)\
                .join((models.ProteinSubunit, models.DNABindingDataset.subunit)).filter(condition))

        return q

class FlaskDNAtoProteinInteractionQueryGenerator(data_query.CachedDataSourceQueryGenerator):
    """ Queries DNA reigons for relevant proteins """

    def __init__(self,
                 taxon=None, max_taxon_dist=None, taxon_dist_scale=None, include_variants=False,
                 temperature=37., temperature_std=1.,
                 ph=7.5, ph_std=0.3):
        """
        Args:
            taxon (:obj:`str`, optional): target taxon
            max_taxon_dist (:obj:`int`, optional): maximum taxonomic distance to include
            taxon_dist_scale (:obj:`float`, optional): The scale of the taxonomic distance scoring distribution.
                This determines how quickly the score falls to zero away from zero.
            include_variants (:obj:`bool`, optional): if :obj:`True`, also include observations from mutant taxa
            temperature (:obj:`float`, optional): desired temperature to search for
            temperature_std (:obj:`float`, optional): how much to penalize observations from other temperatures
            ph (:obj:`float`, optional): desired pH to search for
            ph_std (:obj:`float`, optional): how much to penalize observations from other pHs
        """
        super(FlaskDNAtoProteinInteractionQueryGenerator, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=flask_common_schema.FlaskCommonSchema())


    def get_observed_values(self, DnaSpecie):
        """ Find observed concentrations for the metabolite or similar metabolites

        Args:
            specie (:obj:`data_model.DnaSpecie`): species to find data for

        Returns:
            :obj:`list` of :obj:`data_model.Observable`: list of observable
        """
        proteins = self.get_protein_by_DNA_sequence(DnaSpecie.sequence)

        observable = []
        for subunits in proteins:

            protein = data_model.ProteinSpecie(uniprot_id = subunits[0].uniprot_id,\
                entrez_id = subunits[0].entrez_id, gene_name = subunits[0].gene_name,\
                length = subunits[0].length, mass = subunits[0].mass)
            protein.cross_references = []
            for doc in subunits[0]._metadata.resource:
                protein.cross_references.append(data_model.Resource(namespace = 'pubmed', id = doc._id))

            interaction = data_model.Interaction(name = 'Transcription Factor DNA Binding Site', \
                position = subunits[1], score = subunits[2])

            observable.append(data_model.Observable(specie = protein, interaction = interaction))

        return observable


    def get_protein_by_DNA_sequence(self,sequence, select = models.ProteinSubunit):
        """
        Args:
            sequence (:obj:`data_model.DnaSpecie.sequence`): sequence of DNA segment

        Returns:
            :obj:`list` of :obj:`tuple`: Returns the query for a protein, sequence position, and score

        """
        #TODO: Make more efficient. Add gene location filter

        ans = []

        size = len(sequence)
        all_matricies = self.data_source.session.query(models.DNABindingDataset).all()
        all_sequences = []
        for matricies in all_matricies:
            if len(matricies.dna_binding_data) == size:
                binding_matrix = []
                for position in matricies.dna_binding_data:
                    binding_matrix.append([position.frequency_a, position.frequency_c, position.frequency_g, position.frequency_t])
                binding_matrix = map(list, zip(*binding_matrix))
                self.cache_dirname = tempfile.mkdtemp()
                with open(self.cache_dirname+'/data.pfm', 'w') as pfm:
                    for items in binding_matrix:
                        writer = csv.writer(pfm, delimiter = '\t')
                        writer.writerow(items)
                m = motifs.read(open(self.cache_dirname+'/data.pfm'), 'pfm')
                my_seq = Seq(sequence, IUPAC.unambiguous_dna)
                ##TODO: Include selective threshold
                # distribution = m.pssm.distribution().threshold_paster()
                # print distribution.threshold_paster()
                for position, score in m.pssm.search(my_seq, threshold = 2):
                    ans.append((matricies.subunit, position, score))

        return ans
