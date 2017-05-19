"""
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.core import data_query
from kinetic_datanator.core import filter
from kinetic_datanator.data_source import sabio_rk
from kinetic_datanator.util import molecule_util
from wc_utils.util import string
import sqlalchemy
import sqlalchemy.orm


class ReactionKineticsQuery(data_query.CachedDataSourceQuery):
    """ Finds relevant kinetics observations for reactions

    1. Find kinetics observed for the reaction or similar reactions

      a. Find kinetics observed for the reaction

        i. Find the SABIO-RK compound(s) of each participant
        ii. Find the SABIO-RK reaction(s) which contain all of these SABIO-RK compounds
        iii. Find the SABIO-RK kinetic laws associated with these SABIO-RK reactions

      b. Find kinetics observed for similar reactions

        i. Find kinetics observed for the assigned EC number(s)
        ii. Find kinetics observed for EC number(s) predicted by tools such as E-zyme

    2. Filter out observations from disimilar genetic and environmental conditions and
       rank the remaing observations by their similarity to the desired genetic and environmental
       conditions
    3. Calculate a statistical representation of the relevant observations
    """

    def __init__(self,
                 reaction=None,
                 taxon=None, max_taxon_dist=None, taxon_dist_scale=None, include_variants=False,
                 temperature=37., temperature_std=1.,
                 ph=7.5, ph_std=0.3):
        """
        Args:
            reaction (:obj:`data_model.Reaction`, optional): reaction to find data for
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
        super(ReactionKineticsQuery, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=sabio_rk.SabioRk())

        if reaction:
            self.filters.append(filter.ReactionSimilarityFilter(reaction))

    def get_observed_values(self, reaction):
        """ Find observed kinetics for the reaction or similar reactions

        1. Find kinetics observed for the reaction

          a. Find the SABIO-RK compound(s) of each participant
          b. Find the SABIO-RK reaction(s) which contain all of these SABIO-RK compounds
          c. Find the SABIO-RK kinetic laws associated with these SABIO-RK reactions

        2. Find kinetics observed for similar reactions

          a. Find kinetics observed for the assigned EC number(s)
          b. Find kinetics observed for EC number(s) predicted by tools such as E-zyme

        Args:
            reaction (:obj:`data_model.Reaction`): reaction to find data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observed values
        """

        laws = self.get_kinetic_laws_by_reaction(reaction)
        observed_vals = []
        for law in laws.all():
            observation = data_model.Observation()
            for parameter in law.parameters:
                observable = data_model.Property(name=parameter.type, parent=reaction)

                if parameter.compound:
                    observable.parent.parent = data_model.Specie(
                        name=parameter.compound.name,
                        structure=parameter.compound.structures[0] if parameter.compound.structures else None,
                    )
                if parameter.compartment:
                    observable.parent.parent.parent = data_model.Compartment(
                        name=parameter.compartment.name,
                        parent=reaction,
                    )
  
                observed_vals.append(data_model.ObservedValue(
                    observation=observation,
                    observable=observable,
                    value=parameter.value,
                    error=None,
                    units=parameter.units,
                ))

        return observed_vals

    def get_kinetic_laws_by_reaction(self, reaction, select=sabio_rk.KineticLaw):
        """ Get kinetic laws that were observed for similar reactions (same participants or same EC class)

        Args:
            reaction (:obj:`data_model.Reaction`): reaction to find data for
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`sabio_rk.KineticLaw` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for kinetic laws observed for similar reactions
        """

        # by participants
        participants = reaction.participants
        q = self.get_kinetic_laws_by_participants(participants, select=select)
        if q.count():
            return q

        # by assigned EC numbers
        ec_numbers = [xr.id for xr in reaction.get_manual_ec_numbers()]
        q = self.get_kinetic_laws_by_ec_numbers(ec_numbers, select=select)
        if q.count():
            return q

        # by predicted EC numbers
        ec_numbers = [xr.id for xr in reaction.get_predicted_ec_numbers()]
        q = self.get_kinetic_laws_by_ec_numbers(ec_numbers, select=select)
        if q.count():
            return q

        # return empty list if no relevant observations were found
        return self.data_source.session.query(select).filter_by(id=-1)

    def get_kinetic_laws_by_participants(self, participants, only_formula_and_connectivity=True, include_water_hydrogen=False,
                                         select=sabio_rk.KineticLaw):
        """ Get kinetic laws which were observed for reactions with the participants :obj:`participants`

        Args:
            participants (:obj:`list` of :obj:`data_model.ReactionParticipant`): list of reaction participants
            only_formula_and_connectivity (:obj:`bool`, optional): if :obj:`True`, find reactions which contain species with the same
                InChI formula and connectivity layers
            include_water_hydrogen (:obj:`bool`, optional): if :obj:`True`, restrict reactions based on their water, hydroxide, and
                hydrogen participants
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional): 
                :obj:`sabio_rk.KineticLaw` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for kinetic laws that contain all of the participants
        """
        q_reactions = self.get_reactions_by_participants(
            participants,
            only_formula_and_connectivity=only_formula_and_connectivity,
            include_water_hydrogen=include_water_hydrogen,
            select=sabio_rk.Reaction)
        return self.data_source.session \
            .query(select) \
            .join((q_reactions.subquery(), sabio_rk.KineticLaw.reaction))

    def get_reactions_by_participants(self, participants, only_formula_and_connectivity=True, include_water_hydrogen=False,
                                      select=sabio_rk.Reaction):
        """ Get reactions with the participants :obj:`participants`

        Args:
            participants (:obj:`list` of :obj:`data_model.ReactionParticipant`): list of reaction participants
            only_formula_and_connectivity (:obj:`bool`, optional): if :obj:`True`, find reactions which contain species with the same
                InChI formula and connectivity layers
            include_water_hydrogen (:obj:`bool`, optional): if :obj:`True`, restrict reactions based on their water, hydroxide, and
                hydrogen participants
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`sabio_rk.Reaction` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for reactions that contain all of the participants
        """
        q_reactions = None
        for i_part, part in enumerate(participants):
            try:
                structure = part.specie.to_inchi(only_formula_and_connectivity=only_formula_and_connectivity)
            except ValueError:
                return self.data_source.session.query(select).filter(sabio_rk.Reaction.id == -1)

            if not include_water_hydrogen:
                if only_formula_and_connectivity:
                    formula_and_connectivity = structure
                else:
                    formula_and_connectivity = part.specie.to_inchi(only_formula_and_connectivity=True)
                if formula_and_connectivity in ['', 'H2O', 'H2O']:
                    continue

            if part.coefficient < 0:
                role = 'reactant'
            elif part.coefficient > 0:
                role = 'product'

            q_part = self.get_reactions_by_compound(
                structure, only_formula_and_connectivity=only_formula_and_connectivity, role=role,
                select=select)

            if not q_reactions:
                q_reactions = q_part
            else:
                q_reactions = q_reactions.intersect(q_part)

        return q_reactions

    def get_reactions_by_compound(self, structure, only_formula_and_connectivity=False, role='reactant', select=sabio_rk.Reaction):
        """ Get reactions that contain a structure in role :obj:`role`

        Args:
            structure (:obj:`str`): InChI structure or formula and connectivity layers to search for
            only_formula_and_connectivity (:obj:`bool`, optional): if :obj:`True`, find reactions which contain species with the same
                InChI formula and connectivity layers
            role (:obj:`str`, optional): reaction role (reactant, or product) to search for species
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`sabio_rk.Reaction` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for reactions that contain the structure in role :obj:`role`
        """
        if only_formula_and_connectivity:
            condition = sabio_rk.CompoundStructure._value_inchi_formula_connectivity == structure
        else:
            condition = sabio_rk.CompoundStructure._value_inchi == structure

        if role == 'reactant':
            participant_type = sabio_rk.Reaction.reactants
        else:
            participant_type = sabio_rk.Reaction.products

        return self.data_source.session.query(select) \
            .join((sabio_rk.ReactionParticipant, participant_type)) \
            .join((sabio_rk.Compound, sabio_rk.ReactionParticipant.compound)) \
            .join((sabio_rk.CompoundStructure, sabio_rk.Compound.structures)) \
            .filter(condition) \
            .distinct(sabio_rk.Reaction.id)

    def get_compounds_by_structure(self, inchi, only_formula_and_connectivity=True, select=sabio_rk.Compound):
        """ Get compounds with the same structure. Optionally, get compounds which only have
        the same core empirical formula and core atom connecticity (i.e. same InChI formula
        and connectivity layers).

        Args:
            inchi (:obj:`str`): molecule structure in InChI format
            only_formula_and_connectivity (:obj:`bool`, optional): if :obj:`True`, get compounds which only have
                the same core empirical formula and core atom connecticity. if :obj:`False`, get compounds with the 
                identical structure.
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`sabio_rk.Compound` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching compounds
        """
        q = self.data_source.session.query(select).join((sabio_rk.CompoundStructure, sabio_rk.Compound.structures))
        if only_formula_and_connectivity:
            formula_and_connectivity = molecule_util.InchiMolecule(inchi).get_formula_and_connectivity()
            condition = sabio_rk.CompoundStructure._value_inchi_formula_connectivity == formula_and_connectivity
        else:
            condition = sabio_rk.CompoundStructure._value_inchi == inchi
        return q.filter(condition)

    def get_kinetic_laws_by_ec_numbers(self, ec_numbers, match_levels=4, select=sabio_rk.KineticLaw):
        """ Get kinetic laws which have one of a list of EC numbers or, optionally,
        belong to one of a list of EC classes.

        Args:
            ec_numbers (:obj:`list` of :obj:`str`): EC numbers to search for
            match_levels (:obj:`int`): number of EC levels that the EC number must match
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`sabio_rk.KineticLaw` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching kinetic laws
        """
        # find kinetic laws with identical EC number
        q = self.data_source.session.query(select) \
            .join((sabio_rk.Resource, sabio_rk.KineticLaw.cross_references)) \
            .filter(sabio_rk.Resource.namespace == 'ec-code')

        if match_levels == 4:
            result = q.filter(sabio_rk.Resource.id.in_(ec_numbers))
        else:
            conditions = []
            for ec_number in ec_numbers:
                ec_class, _, _ = string.partition_nth(ec_number, '.', match_levels)
                conditions.append(sabio_rk.Resource.id.like('{}.%'.format(ec_class)))
            result = q.filter(sqlalchemy.or_(*conditions))

        return result
