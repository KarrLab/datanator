# -*- coding: utf-8 -*-

"""
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-09-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.core import data_query
from kinetic_datanator.core import common_schema

class ProteinConcentrationsQueryGenerator(data_query.CachedDataSourceQueryGenerator):
    """ Finds relevant concentration observations for proteins """

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
        super(ProteinConcentrationsQueryGenerator, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=common_schema.CommonSchema())

    def get_observed_values(self, protein):
        """ Find the observed values for protein abundance

        Args:
            protein (:obj:`data_model.ProteinSpecie`): Protein to find data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observed values

        """
        pass


    def get_abundance_by_uniprot(self, uniprot, select = common_schema.AbundanceData):
        """ Find the abundance from uniprot

        Args:
            uniprot (:obj:`str`): protein id from Uniprot Database
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`common_schema.CommonSchema` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching abundance rows
        """
        q = self.data_source.session.query(select).join(common_schema.ProteinSubunit, common_schema.AbundanceData.subunit)
        condition = common_schema.ProteinSubunit.uniprot_id == uniprot
        return q.filter(condition)

    def get_abundance_by_gene_name(self, gene_name, select = common_schema.AbundanceData):
        """

        """
        ##TODO: Figure out a gene name from the string of gene_name in the common database. So if name = gen_name (in string)
        q = self.data_source.session.query(select).join(common_schema.ProteinSubunit, common_schema.AbundanceData.subunit)
        condition = common_schema.ProteinSubunit.gene_name == gene_name
        return q.filter(condition)


    def get_abundance_by_sequence(self, sequence, select = common_schema.AbundanceData):
        """

        """
        q = self.data_source.session.query(select).join(common_schema.ProteinSubunit, common_schema.AbundanceData.subunit)
        condition = common_schema.ProteinSubunit.canonical_sequence == sequence
        return q.filter(condition)

    def get_abundance_by_entrez(self, entrez_id, select = common_schema.AbundanceData):
        """

        """
        q = self.data_source.session.query(select).join(common_schema.ProteinSubunit, common_schema.AbundanceData.subunit)
        condition = common_schema.ProteinSubunit.entrez_id == entrez_id
        return q.filter(condition)

    def get_abundance_by_mass(self, mass, select = common_schema.AbundanceData):
        """

        """
        q = self.data_source.session.query(select).join(common_schema.ProteinSubunit, common_schema.AbundanceData.subunit)
        condition = common_schema.ProteinSubunit.mass == mass
        return q.filter(condition)

    def get_abundance_by_length(self, length, select = common_schema.AbundanceData):
        """

        """
        q = self.data_source.session.query(select).join(common_schema.ProteinSubunit, common_schema.AbundanceData.subunit)
        condition = common_schema.ProteinSubunit.length == length
        return q.filter(condition)
