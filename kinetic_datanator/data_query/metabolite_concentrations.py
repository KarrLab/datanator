"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-05-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.core import data_query
from kinetic_datanator.core import common_schema
from kinetic_datanator.util import molecule_util
from kinetic_datanator.flask_datanator import flask_common_schema, models



class MetaboliteConcentrationsQueryGenerator(data_query.CachedDataSourceQueryGenerator):
    """ Finds relevant concentration observations for metabolites """

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
        super(MetaboliteConcentrationsQueryGenerator, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=common_schema.CommonSchema())

        self.filters.append(data_query.SpecieStructuralSimilarityFilter())
        # self.filters.append(data_query.MolecularSimilarityFilter())

    def get_observed_values(self, specie):
        """ Find observed concentrations for the metabolite or similar metabolites

        Args:
            specie (:obj:`data_model.Specie`): species to find data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observations
        """
        concentrations = self.get_concentration_by_structure(specie.structure).all()
        observed_vals = []

        for c in concentrations:

        #TODO: Figure out the multi metadata issue.

            observation = data_model.Observation(
                genetics = data_model.Genetics(
                    taxon = c._metadata.taxon[0].name,
                    variation = c._metadata.cell_line[0].name
                ),
                environment = data_model.Environment(
                    temperature = c._metadata.conditions[0].temperature,
                    ph = c._metadata.conditions[0].ph,
                    media = c._metadata.conditions[0].media,
                    growth_status = c._metadata.conditions[0].growth_status,
                    growth_system = c._metadata.conditions[0].growth_system,
                )
            )

            observable = data_model.Observable(
                specie = specie,
                compartment = data_model.Compartment(name = c._metadata.cell_compartment[0].name)
            )

            observed_vals.append(data_model.ObservedValue(
                observation = observation,
                observable = observable,
                value = c.value,
                error = c.error,
                #TODO: Confirm the units are the same
                units = 'uM'
            ))

        return observed_vals



    def get_concentration_by_structure(self, inchi, only_formula_and_connectivity=True, select = common_schema.Concentration):
        """

        Args:
            specie (:obj:`data_model.Specie`): species to find data for

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching concentrations

        """

        q = self.data_source.session.query(select).join((common_schema.Compound, select.compound)).\
            join((common_schema.Structure, common_schema.Compound.structure))

        if only_formula_and_connectivity:
            formula_and_connectivity = molecule_util.InchiMolecule(inchi).get_formula_and_connectivity()
            condition = common_schema.Structure._structure_formula_connectivity == formula_and_connectivity
        else:
            condition = common_schema.Structure._value_inchi == inchi
        return q.filter(condition)


class FlaskMetaboliteConcentrationsQueryGenerator(data_query.CachedDataSourceQueryGenerator):
    """ Finds relevant concentration observations for metabolites """

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
        super(FlaskMetaboliteConcentrationsQueryGenerator, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=flask_common_schema.FlaskCommonSchema())

        self.filters.append(data_query.SpecieStructuralSimilarityFilter())
        # self.filters.append(data_query.MolecularSimilarityFilter())

    def get_observed_values(self, compound):
        """ Find observed concentrations for the metabolite or similar metabolites

        Args:
            specie (:obj:`data_model.Specie`): species to find data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observations
        """
        concentrations = self.get_concentration_by_structure(compound.structure._value_inchi, only_formula_and_connectivity=False).all()
        observed_vals = []

        for c in concentrations:

        #TODO: Figure out the multi metadata issue.

            observation = data_model.Observation(
                genetics = data_model.Genetics(
                    taxon = c._metadata.taxon[0].name,
                    variation = c._metadata.cell_line[0].name
                ),
                environment = data_model.Environment(
                    temperature = c._metadata.conditions[0].temperature,
                    ph = c._metadata.conditions[0].ph,
                    media = c._metadata.conditions[0].media,
                    growth_status = c._metadata.conditions[0].growth_status,
                    growth_system = c._metadata.conditions[0].growth_system,
                )
            )

            observable = data_model.Observable(
                specie = data_model.Specie(name = compound.compound_name, structure=compound.structure._value_inchi),
                compartment = data_model.Compartment(name = c._metadata.cell_compartment[0].name)
            )

            observed_vals.append(data_model.ObservedValue(
                observation = observation,
                observable = observable,
                value = c.value,
                error = c.error,
                #TODO: Confirm the units are the same
                units = 'uM'
            ))

        return observed_vals



    def get_concentration_by_structure(self, inchi, only_formula_and_connectivity=True, select = models.Concentration):
        """

        Args:
            specie (:obj:`data_model.Specie`): species to find data for

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching concentrations

        """

        q = self.data_source.session.query(select).join((models.Compound, select.compound)).\
            join((models.Structure, models.Compound.structure))

        if only_formula_and_connectivity:
            formula_and_connectivity = molecule_util.InchiMolecule(inchi).get_formula_and_connectivity()
            condition = models.Structure._structure_formula_connectivity == formula_and_connectivity
        else:
            condition = models.Structure._value_inchi == inchi
        return q.filter(condition)
