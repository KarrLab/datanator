"""
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.core import data_query
from kinetic_datanator.data_source import sabio_rk
from kinetic_datanator.util import molecule_util
from wc_utils.util import string
import sqlalchemy
import sqlalchemy.orm


class ReactionKineticsQueryGenerator(data_query.CachedDataSourceQueryGenerator):
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
        super(ReactionKineticsQueryGenerator, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=sabio_rk.SabioRk())

        self.filters.append(data_query.ReactionSimilarityFilter())
        self.filters.append(data_query.ReactionParticipantFilter())

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

        q_law = self.get_kinetic_laws_by_reaction(reaction)
        observed_vals = []
        for law in q_law.all():
            sabiork_reaction_id = next(xr.id for xr in law.cross_references if xr.namespace == 'sabiork.reaction')

            reaction = data_model.Reaction(
                cross_references=[
                    data_model.Resource(namespace='sabiork.kineticrecord', id=str(law.id)),
                    data_model.Resource(namespace='sabiork.reaction', id=sabiork_reaction_id),
                ],
            )
            species = {}
            compartments = {}

            for reactant in law.reactants:
                part = data_model.ReactionParticipant(coefficient=-1)

                if reactant.compound.id not in species:
                    species[reactant.compound.id] = data_model.Specie(name=reactant.compound.name)
                part.specie = species[reactant.compound.id]

                if reactant.compound.structures:
                    part.specie.structure = reactant.compound.structures[0].value

                if reactant.compartment:
                    if reactant.compartment.name not in compartments:
                        compartments[reactant.compartment.name] = data_model.Compartment(name=reactant.compartment.name)
                    part.compartment = compartments[reactant.compartment.name]

                reaction.participants.append(part)

            for product in law.products:
                part = data_model.ReactionParticipant(coefficient=1)

                if product.compound.id not in species:
                    species[product.compound.id] = data_model.Specie(name=product.compound.name)
                part.specie = species[product.compound.id]

                if product.compound.structures:
                    part.specie.structure = product.compound.structures[0].value

                if product.compartment:
                    if product.compartment.name not in compartments:
                        compartments[product.compartment.name] = data_model.Compartment(name=product.compartment.name)
                    part.compartment = compartments[product.compartment.name]

                reaction.participants.append(part)

            for modifier in law.modifiers:
                part = data_model.ReactionParticipant(coefficient=0)

                if modifier.compound.id not in species:
                    species[modifier.compound.id] = data_model.Specie(name=modifier.compound.name)
                part.specie = species[modifier.compound.id]

                if modifier.compound.structures:
                    part.specie.structure = modifier.compound.structures[0].value

                if modifier.compartment:
                    if modifier.compartment.name not in compartments:
                        compartments[modifier.compartment.name] = data_model.Compartment(name=modifier.compartment.name)
                    part.compartment = compartments[modifier.compartment.name]

                reaction.participants.append(part)

            observation = data_model.Observation(
                genetics=data_model.Genetics(
                    taxon=law.taxon,
                    variation=law.taxon_variant,
                ),
                environment=data_model.Environment(
                    temperature=law.temperature,
                    ph=law.ph,
                    media=law.media,
                ),
            )
            for parameter in law.parameters:
                if parameter.value is None:
                    continue

                if parameter.type not in sabio_rk.Parameter.TYPES:
                    continue

                observable = data_model.Observable(
                    interaction=reaction,
                    property=sabio_rk.Parameter.TYPES[parameter.type],
                )

                if parameter.compound:
                    observable.specie = species[parameter.compound.id]
                    if parameter.compartment:
                        observable.compartment = data_model.Compartment(
                            id=parameter.compartment.name,
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
        """ Get kinetic laws with the participants :obj:`participants`

        Args:
            participants (:obj:`list` of :obj:`data_model.ReactionParticipant`): list of reaction participants
            only_formula_and_connectivity (:obj:`bool`, optional): if :obj:`True`, find kinetic laws which contain species with the same
                InChI formula and connectivity layers
            include_water_hydrogen (:obj:`bool`, optional): if :obj:`True`, restrict kinetic laws based on their water, hydroxide, and
                hydrogen participants
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`sabio_rk.KineticLaw` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for kinetic laws that contain all of the participants
        """
        q_laws = None
        for i_part, part in enumerate(participants):
            try:
                structure = part.specie.to_inchi(only_formula_and_connectivity=only_formula_and_connectivity)
            except ValueError:
                return self.data_source.session.query(select).filter(sabio_rk.KineticLaw.id == -1)

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

            q_part = self.get_kinetic_laws_by_compound(
                structure, only_formula_and_connectivity=only_formula_and_connectivity, role=role,
                select=select)

            if not q_laws:
                q_laws = q_part
            else:
                q_laws = q_laws.intersect(q_part)

        return q_laws

    def get_kinetic_laws_by_compound(self, structure, only_formula_and_connectivity=False, role='reactant', select=sabio_rk.KineticLaw):
        """ Get kinetic laws that contain a structure in role :obj:`role`

        Args:
            structure (:obj:`str`): InChI structure or formula and connectivity layers to search for
            only_formula_and_connectivity (:obj:`bool`, optional): if :obj:`True`, find kinetic laws which contain species with the same
                InChI formula and connectivity layers
            role (:obj:`str`, optional): role (reactant, or product) to search for species
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`sabio_rk.KineticLaw` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for kinetic laws that contain the structure in role :obj:`role`
        """

        if only_formula_and_connectivity:
            condition = sabio_rk.CompoundStructure._value_inchi_formula_connectivity == structure
        else:
            condition = sabio_rk.CompoundStructure._value_inchi == structure

        if role == 'reactant':
            participant_type = sabio_rk.ReactionParticipant.reactant_kinetic_law_id
        else:
            participant_type = sabio_rk.ReactionParticipant.product_kinetic_law_id

        session = self.data_source.session

        q_structure = session.query(sabio_rk.CompoundStructure) \
            .filter(condition) \
            .subquery()

        return session.query(select) \
            .join(sabio_rk.ReactionParticipant,
                  sabio_rk.KineticLaw._id == participant_type) \
            .join(sabio_rk.compound_compound_structure,
                  sabio_rk.ReactionParticipant.compound_id == sabio_rk.compound_compound_structure.c.get('compound__id')) \
            .join(q_structure,
                  sabio_rk.compound_compound_structure.c.get('compound_structure__id') == q_structure.c._id) \
            .filter(participant_type != None) \
            .distinct(sabio_rk.KineticLaw._id)

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
