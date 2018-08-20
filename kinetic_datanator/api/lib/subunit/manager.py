"""
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-20
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model, data_query, common_schema, models
from kinetic_datanator.api.lib.data_manager import BaseManager
from kinetic_datanator.util.constants import DATA_CACHE_DIR

class ProteinSubunitManager(BaseManager):
    """ Manages subunit information for API """

    def __init__(self, cache_dirname = DATA_CACHE_DIR):
        self.data_source = common_schema.CommonSchema(cache_dirname=cache_dirname)

    def get_subunit_by_id(self, id):
        return self.data_source.session.query(models.ProteinSubunit).get(id)

    def _search(self, value):
        return models.ProteinSubunit.query.search(value).all()

    def _port(self, protein):
        return data_model.ProteinSpecie(name=protein.subunit_name,
                                        uniprot_id=protein.uniprot_id, entrez_id=protein.entrez_id,
                                        gene_name=protein.gene_name, length=protein.length,
                                        mass=protein.mass, sequence = protein.canonical_sequence)


    def get_observed_abundances(self, protein):
        """ Find the observed values for protein abundance

        Args:
            protein (:obj:`models.ProteinSubunit`): Protein Subunit to find data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observed values

        """
        abundances = self.get_abundance_by_uniprot(protein.uniprot_id)
        observed_vals = []

        for abundance in abundances:
            metadata = self.metadata_dump(abundance.dataset)

            observable = data_model.Observable(
                specie=self._port(protein)
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
        """ Find the abundance from a uniprot id

        Args:
            uniprot (:obj:`str`): protein id from Uniprot Database

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
            gene_name (:obj:`str`): gene name for a given protein

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
            sequence (:obj:`str`): amino acid sequence for a given protein

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
            entrez_id (:obj:`str`): NCBI entrez id for a given protein

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
            mass (:obj:`int`): mass of a protein

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
            length (:obj:`str`): number of amino acids in a protein

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching abundance rows
        """
        q = self.data_source.session.query(select).join(
            models.ProteinSubunit, models.AbundanceData.subunit)
        condition = models.ProteinSubunit.length == length
        return q.filter(condition)


subunit_manager = ProteinSubunitManager()
