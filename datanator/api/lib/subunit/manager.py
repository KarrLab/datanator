"""
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-20
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.core import data_model, data_query, common_schema, models
from datanator.api.lib.data_manager import BaseManager
from datanator.util.constants import DATA_CACHE_DIR
from datanator.api.lib.complex.manager import complex_manager
from sqlalchemy import or_

class ProteinSubunitManager(BaseManager):
    """ Manages subunit information for API """

    def __init__(self, cache_dirname = DATA_CACHE_DIR):
        self.data_source = common_schema.CommonSchema(cache_dirname=cache_dirname)

    def get_subunit_by_id(self, id):
        return self.data_source.session.query(models.ProteinSubunit).get(id)

    def _search_simple(self, value):
        return models.ProteinSubunit.query.search(value, vector=models.ProteinSubunit.simple_search_vector).all()


    def _search_complex(self, value):
        return models.ProteinSubunit.query.search(value, vector=models.ProteinSubunit.complex_search_vector).all()


    def _port(self, protein):
        return data_model.ProteinSpecie(name=protein.subunit_name,
                                        uniprot_id=protein.uniprot_id, entrez_id=protein.entrez_id,
                                        gene_name=protein.gene_name, length=protein.length,
                                        mass=protein.mass, sequence = protein.canonical_sequence)

    def get_observable_interactions(self, protein_subunit):
        """ Get known protein interactions that were observed for a given subunit

        Args:
            protein_subunit (:obj:`models.ProteinSubunit`): subunit to find interactions for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedInteraction`: list of Observed Protein Interactions

        """
        interact = self.get_interaction_by_subunit(protein_subunit.uniprot_id).all()
        observed_interaction= []

        for item in interact:

            participants = []
            for participant, protein in [(item.type_a, item.protein_a), (item.type_b, item.protein_b)]:
                if participant == None:
                    specie = None

                if participant == 'protein':
                    prot = self.data_source.session.query(models.ProteinSubunit).filter_by(uniprot_id = protein).first()
                    if prot:
                        specie =  data_model.ProteinSpecie(name=prot.uniprot_id, uniprot_id =prot.uniprot_id,
                            entrez_id=prot.entrez_id, gene_name=prot.gene_name, length=prot.length,
                            mass=prot.mass, sequence=prot.canonical_sequence)
                elif participant == 'peptide':
                    specie = data_model.PolymerSpecie(name=protein, sequence=protein)

                participants.append(specie)

            metadata = self.metadata_dump(item)
            resource = [data_model.Resource(namespace=source.namespace, id=source._id) for source in item._metadata.resource]
            interaction = data_model.SpecieInteraction(specie_a=participants[0], specie_b=participants[1], stoichiometry_a = item.stoich_a, stoichiometry_b = item.stoich_b,
            loc_a=item.loc_a, loc_b=item.loc_b, cross_references=resource, name=item.name, confidence=item.confidence, type_a = item.type_a, type_b=item.type_b, interaction_type=item.interaction_type)
            observed_interaction.append(data_model.ObservedInteraction(interaction=interaction, metadata=metadata))


        return observed_interaction

    def get_observable_complex(self, protein_subunit):
        """ Get known protein complex that were observed for a given subunit

        Args:
            protein_subunit (:obj:`models.ProteinSubunit`): subunit to find complex for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedSpecie`: list of Protein Complexes
        """
        observed_specie = []
        _complex = self.get_known_complex_by_subunit(protein_subunit.uniprot_id)

        for complex in _complex:
            metadata= self.metadata_dump(complex)
            ported_complex = complex_manager._port(complex)

            observed_specie.append(data_model.ObservedSpecie(specie = ported_complex, metadata=metadata))

        return observed_specie

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


    def get_interaction_by_subunit(self, uniprot, select = models.ProteinInteraction):
        """ Get interactions that were observed for a given uniprot id

        Args:
            uniprot (:obj:`str`): uniprot id to search for

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for protein interactions that contain the uniprot_id
        """
        return self.data_source.session.query(select).filter(or_(select.protein_a == uniprot,
                select.protein_b == uniprot))

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


    def get_known_complex_by_subunit(self, uniprot, select = models.ProteinComplex):
        """ Get known complexes that were observed for a given uniprot id subunit

        Args:
            uniprot (:obj:`str`): uniprot id to search for

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for protein complexes that contain the uniprot_id
        """
        q = self.data_source.session.query(select).join(models.ProteinSubunit, select.protein_subunit)
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

    # def get_abundance_by_entrez(self, entrez_id, select=models.AbundanceData):
    #     """ Find the abundance from uniprot
    #
    #     Args:
    #         entrez_id (:obj:`str`): NCBI entrez id for a given protein
    #
    #     Returns:
    #         :obj:`sqlalchemy.orm.query.Query`: query for matching abundance rows
    #     """
    #     q = self.data_source.session.query(select).join(
    #         models.ProteinSubunit, models.AbundanceData.subunit)
    #     condition = models.ProteinSubunit.entrez_id == entrez_id
    #     return q.filter(condition)

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
