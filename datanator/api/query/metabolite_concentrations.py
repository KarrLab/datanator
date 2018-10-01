"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-05-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model, data_query, common_schema, models
from kinetic_datanator.util import molecule_util



class MetaboliteConcentrationQuery(data_query.CachedDataSourceQueryGenerator):
    """ Finds relevant concentration observations for metabolites """

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
        super(MetaboliteConcentrationQuery, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=common_schema.CommonSchema(cache_dirname= cache_dirname))

        # self.filters.append(data_query.SpecieStructuralSimilarityFilter())
        # self.filters.append(data_query.MolecularSimilarityFilter())

    def get_observed_result(self, metabolite):
        """ Find observed concentrations for the metabolite or similar metabolites

        Args:
            metabolite (:obj:`models.Metabolite`): metabolite to find data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observations
        """

        concentrations = self.get_concentration_by_structure(metabolite.structure._value_inchi, only_formula_and_connectivity=False).all()
        observed_values = []

        references = [data_model.Resource(namespace=item.namespace, id=item._id) for item in metabolite._metadata.resource]

        for c in concentrations:
            metadata = self.metadata_dump(c)

            observable = data_model.Observable(
                specie = data_model.Specie(name = metabolite.metabolite_name,
                cross_references = references , structure=metabolite.structure._value_inchi),
                compartment = data_model.Compartment(name = c._metadata.cell_compartment[0].name)
            )

            observed_values.append(data_model.ObservedValue(
                metadata = metadata,
                observable = observable,
                value = c.value,
                error = c.error,
                units = c.units
            ))

        return observed_values

    def get_concentration_by_structure(self, inchi, only_formula_and_connectivity=True, select = models.Concentration):
        """
        Args:
            inchi (:obj:`str`): inchi structure to find concentrations

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching concentration rows
        """

        q = self.data_source.session.query(select).join((models.Metabolite, select.metabolite)).\
            join((models.Structure, models.Metabolite.structure))

        if only_formula_and_connectivity:
            formula_and_connectivity = molecule_util.InchiMolecule(inchi).get_formula_and_connectivity()
            condition = models.Structure._structure_formula_connectivity == formula_and_connectivity
        else:
            condition = models.Structure._value_inchi == inchi

        return q.filter(condition)
