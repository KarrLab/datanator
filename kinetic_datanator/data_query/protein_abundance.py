"""
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-09-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model, data_query, flask_common_schema, models

class ProteinAbundanceQuery(data_query.CachedDataSourceQueryGenerator):
    """ Finds relevant concentration observations for proteins """

    def __init__(self,
                 taxon=None, max_taxon_dist=None, taxon_dist_scale=None, include_variants=False,
                 temperature=37., temperature_std=1.,
                 ph=7.5, ph_std=0.3, cache_dirname = None):
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

        super(ProteinAbundanceQuery, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=flask_common_schema.FlaskCommonSchema(cache_dirname = cache_dirname))

    def get_observed_result(self, protein):
        """ Find the observed values for protein abundance

        Args:
            protein (:obj:`models.ProteinSubunit`): Protein Subunit to find data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observed values

        """
        abundances = self.get_abundance_by_uniprot(protein.uniprot_id)
        observed_vals = []

        for abundance in abundances:

            metadata = data_model.ObservedResultMetadata(
                genetics=data_model.Genetics(
                    taxon=abundance.dataset._metadata.taxon[0].name
                )
            ) if abundance.dataset._metadata.taxon else None

            observable = data_model.Observable(
                specie=data_model.ProteinSpecie(name=protein.subunit_name,
                                                uniprot_id=protein.uniprot_id, entrez_id=protein.entrez_id,
                                                gene_name=protein.gene_name, length=protein.length,
                                                mass=protein.mass, sequence = protein.canonical_sequence)
            )

            observable.specie.cross_references = [
                data_model.Resource(namespace='publication',
                                    id=abundance.dataset.file_name),
                data_model.Resource(
                    namespace='url', id=abundance.dataset._metadata.resource[0]._id)
            ]

            observed_vals.append(data_model.ObservedValue(
                metadata=metadata,
                observable=observable,
                value=abundance.abundance,
                error=0,
                units='PPM',
            ))

        return observed_vals

    def get_abundance_by_uniprot(self, uniprot, select=models.AbundanceData):
        """ Find the abundance from uniprot

        Args:
            uniprot (:obj:`str`): protein id from Uniprot Database
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`common_schema.CommonSchema` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching abundance rows
        """
        q = self.data_source.session.query(select).join(
            models.ProteinSubunit, models.AbundanceData.subunit)
        condition = models.ProteinSubunit.uniprot_id == uniprot
        return q.filter(condition)

    def get_abundance_by_gene_name(self, gene_name, select=models.AbundanceData):
        """ Find the abundance from gene_name

        Args:
            uniprot (:obj:`str`): protein id from Uniprot Database
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`common_schema.CommonSchema` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching abundance rows
        """

        # TODO: Figure out a gene name from the string of gene_name in the common database. So if name = gen_name (in string)
        q = self.data_source.session.query(select).join(
            models.ProteinSubunit, models.AbundanceData.subunit)
        condition = models.ProteinSubunit.gene_name == gene_name
        return q.filter(condition)

    def get_abundance_by_sequence(self, sequence, select=models.AbundanceData):
        """ Find the abundance from uniprot

        Args:
            uniprot (:obj:`str`): protein id from Uniprot Database
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`common_schema.CommonSchema` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching abundance rows
        """
        q = self.data_source.session.query(select).join(
            models.ProteinSubunit, models.AbundanceData.subunit)
        condition = models.ProteinSubunit.canonical_sequence == sequence
        return q.filter(condition)

    def get_abundance_by_entrez(self, entrez_id, select=models.AbundanceData):
        """ Find the abundance from uniprot

        Args:
            uniprot (:obj:`str`): protein id from Uniprot Database
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`common_schema.CommonSchema` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching abundance rows
        """
        q = self.data_source.session.query(select).join(
            models.ProteinSubunit, models.AbundanceData.subunit)
        condition = models.ProteinSubunit.entrez_id == entrez_id
        return q.filter(condition)

    def get_abundance_by_mass(self, mass, select=models.AbundanceData):
        """ Find the abundance from uniprot

        Args:
            uniprot (:obj:`str`): protein id from Uniprot Database
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`common_schema.CommonSchema` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching abundance rows
        """
        q = self.data_source.session.query(select).join(
            models.ProteinSubunit, models.AbundanceData.subunit)
        condition = models.ProteinSubunit.mass == mass
        return q.filter(condition)

    def get_abundance_by_length(self, length, select=models.AbundanceData):
        """ Find the abundance from uniprot

        Args:
            uniprot (:obj:`str`): protein id from Uniprot Database
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`common_schema.CommonSchema` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching abundance rows
        """
        q = self.data_source.session.query(select).join(
            models.ProteinSubunit, models.AbundanceData.subunit)
        condition = models.ProteinSubunit.length == length
        return q.filter(condition)
