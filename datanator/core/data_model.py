"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.util import molecule_util
import copy
import enum
import numpy
import obj_tables.core
import obj_tables.bio
import six.moves


class ConsensusMethod(enum.Enum):
    """ Represents the method by which a consensus was chosen """
    manual = 0
    mean = 1
    median = 2
    mode = 3
    weighted_mean = 4
    weighted_median = 5
    weighted_mode = 6


class Consensus(obj_tables.core.Model):
    """ Represents a consensus of one or more observed values of an attribute of a component of a model

    Attributes:
        observable (:obj:`Observable`): biological component that was estimated
        value (:obj:`float`): consensus value of the attribute of the model component
        error (:obj:`float`): uncertainty of the value of the attribute of the model component
        units (:obj:`str`): units of the value of the attribute of the model component
        evidence (:obj:`list` of :obj:`Evidence`): list of evidence which the consensus
            value is based on
        method (:obj:`ConsensusMethod`): method used to calculate the consensus value and error
        user (:obj:`str`): user who generated the consensus
        date (:obj:`datetime.datetime`): date and time when the consensus was generated
    """
    observable = obj_tables.core.ManyToOneAttribute('Observable', related_name='consensus')
    value = obj_tables.core.FloatAttribute()
    error = obj_tables.core.FloatAttribute()
    units = obj_tables.core.StringAttribute()
    evidence = obj_tables.core.ManyToManyAttribute('Evidence', related_name='consensus')
    method = obj_tables.core.EnumAttribute(ConsensusMethod)
    user = obj_tables.core.StringAttribute()
    date = obj_tables.core.DateTimeAttribute()


class Evidence(obj_tables.core.Model):
    """ Represents the observed values and their relevance which support a consensus

    Attributes:
        value (:obj:`ObservedValue`): observed value
        relevance (:obj:`float`): numeric score which indicates the relevance of the observed value to the
            consensus
    """
    value = obj_tables.core.ManyToOneAttribute('ObservedValue', related_name='evidence')
    relevance = obj_tables.core.FloatAttribute()


class ObservedResultMetadata(obj_tables.core.Model):
    """ Represents an observation (one or more observed values) about a biological system

    Attributes:
        genetics (:obj:`Genetics`): the taxon, and any genetic variation from the wildtype taxon, that the component was observed in
        environment (:obj:`Environment`): environment that the component was observed in
        values (:obj:`list` of :obj:`ObservedValue`): observed values
        reference (:obj:`Reference`): reference to the reference
    """
    genetics = obj_tables.core.ManyToOneAttribute('Genetics', related_name='observations')
    environment = obj_tables.core.ManyToOneAttribute('Environment', related_name='observations')
    cross_references = obj_tables.core.ManyToManyAttribute('Resource', related_name='observations')
    method = obj_tables.core.ManyToOneAttribute('Method', related_name='observations')
    synonym = obj_tables.core.ManyToManyAttribute('Synonym', related_name='observations')

class ObservedResult(obj_tables.core.Model):
    """ Represents a base dataset for a queried response

    Attributes:
        observation (:obj:`Observaton`): the collection of covariate observed values
        method (:obj:`str`): method that was used to make the observation
    """

    metadata = obj_tables.core.ManyToOneAttribute('ObservedResultMetadata', related_name='observed_result')


class ObservedInteraction(ObservedResult):
    """ Represents an observed interaction of a biological system

    Attributes:
        interaction (:obj:`Interaction`): observed interaction
    """

    interaction = obj_tables.core.ManyToOneAttribute('Interaction', related_name='observed_interaction')

class ObservedSpecie(ObservedResult):
    """ Represents an observed interaction of a biological system

    Attributes:
        specie (:obj:`Specie`): observed specie
    """

    specie = obj_tables.core.ManyToOneAttribute('Specie', related_name='observed_specie')

class ObservedValue(ObservedResult):
    """ Represents an observed value of a biological system

    Attributes:
        observable (:obj:`Observaton`): the observed interaction or specie for which the value corresponds
        value (:obj:`float`): observed value
        error (:obj:`float`): uncertainty of the observed value
        units (:obj:`units`): SI units of the observed value
    """
    observable = obj_tables.core.ManyToOneAttribute('Observable', related_name='observed_values')
    value = obj_tables.core.FloatAttribute()
    error = obj_tables.core.FloatAttribute()
    units = obj_tables.core.StringAttribute()


class Observable(obj_tables.core.Model):
    """ Represents an observable of a biological system

    Attributes:
        interaction (:obj:`Interaction`): observed interaction
        specie (:obj:`Specie`): observed species
        compartment (:obj:`Compartment`): compartment that the spcies/interaction was observed in
        property (:obj:`str`): property that was observed
    """
    interaction = obj_tables.core.ManyToOneAttribute('Interaction', related_name='observed_values')
    specie = obj_tables.core.ManyToOneAttribute('Specie', related_name='observed_values')
    compartment = obj_tables.core.ManyToOneAttribute('Compartment', related_name='observed_values')
    property = obj_tables.core.StringAttribute()


class EntityInteractionOrProperty(obj_tables.core.Model):
    """ Represents an observable of a biological system

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        cross_references (:obj:`list` of :obj:`Resource`): list of cross references to external resources
    """
    id = obj_tables.core.StringAttribute()
    name = obj_tables.core.StringAttribute()
    cross_references = obj_tables.core.ManyToManyAttribute('Resource', related_name='observables')


class Compartment(EntityInteractionOrProperty):
    """ Representes a compartment in a biological system """
    pass


class Specie(EntityInteractionOrProperty):
    """ Represents a molecular species in a biological system

    Attributes:
        structure (:obj:`str`): structure
    """
    structure = obj_tables.core.LongStringAttribute()

    def to_inchi(self, only_formula_and_connectivity=False):
        """ Get the structure in InChi format

        Args:
            only_formula_and_connectivity (:obj:`bool`): if :obj:`True`, return only the
                formula and connectivity layers

        Returns:
            :obj:`str`: structure in InChi format or just the formula and connectivity layers
                if :obj:`only_formula_and_connectivity` is :obj:`True`
        """
        inchi = molecule_util.Molecule(structure=self.structure).to_inchi()
        if only_formula_and_connectivity:
            return molecule_util.InchiMolecule(inchi).get_formula_and_connectivity()
        else:
            return inchi

    def to_mol(self):
        """ Get the structure in .mol format

        Returns:
            :obj:`str`: structure in .mol format
        """
        return molecule_util.Molecule(structure=self.structure).to_mol()

    def to_openbabel(self):
        """ Get the structure as a Open Babel molecule

        Returns:
            :obj:`openbabel.OBMol`: structure as a Open Babel molecule
        """
        return molecule_util.Molecule(structure=self.structure).to_openbabel()

    def to_pybel(self):
        """ Get the structure as a Pybel molecule

        Returns:
            :obj:`pybel.Molecule`: structure as a Pybel molecule
        """
        return molecule_util.Molecule(structure=self.structure).to_pybel()

    def to_smiles(self):
        """ Get the structure in canonical SMILES format

        Returns:
            :obj:`str`: structure in canonical SMILES format
        """
        return molecule_util.Molecule(structure=self.structure).to_smiles()

    def get_similarity(self, other, fingerprint_type='fp2'):
        """ Calculate the similarity with another species

        Args:
            other (:obj:`Specie`): a second species
            fingerprint_type (:obj:`str`, optional): fingerprint type to use to calculate similarity

        Returns:
            :obj:`float`: the similarity with the other molecule
        """
        self_mol = molecule_util.Molecule(structure=self.structure)
        other_mol = molecule_util.Molecule(structure=other.structure)
        return self_mol.get_similarity(other_mol, fingerprint_type=fingerprint_type)


class PolymerSpecie(Specie):
    """ Represents a polymer

    Attributes:
        sequence (:obj:`str`): sequence
    """
    sequence = obj_tables.core.LongStringAttribute()


class DnaSpecie(PolymerSpecie):
    """ Represents a DNA polymer

    Attributes:
        binding_matrix (:obj:`Bio.motifs.matrix.FrequencyPositionMatrix`): Binding motif
    """
    binding_matrix = obj_tables.bio.FrequencyPositionMatrixAttribute()


class RnaSpecie(PolymerSpecie):
    """ Represents a RNA polymer """
    pass


class ProteinSpecie(PolymerSpecie):
    """ Represents a protein polymer

    Attributes:
        uniprot_id (:obj:`str`): Uniprot Identifier
        entrez_id (:obj:`int`): Entrez Identifier
        gene_name (:obj:`str`): gene name from which protein stems
        length (:obj:`int`): Length of the amino acid sequence
        mass (:obj:`int`): Mass of the protein in KDa

    """
    uniprot_id  = obj_tables.core.StringAttribute()
    entrez_id = obj_tables.core.IntegerAttribute()
    gene_name = obj_tables.core.StringAttribute()
    length = obj_tables.core.IntegerAttribute()
    mass = obj_tables.core.IntegerAttribute()

class ProteinComplexSpecie(ProteinSpecie):
    """ Represents a protein interaction

        Attributes:

    """

    go_id = obj_tables.core.StringAttribute()
    go_dsc = obj_tables.core.StringAttribute()
    funcat_id = obj_tables.core.StringAttribute()
    funcat_dsc = obj_tables.core.StringAttribute()
    su_cmt = obj_tables.core.StringAttribute()
    complex_cmt = obj_tables.core.StringAttribute()
    disease_cmt = obj_tables.core.StringAttribute()
    class_name = obj_tables.core.StringAttribute()
    family_name = obj_tables.core.StringAttribute()
    molecular_weight = obj_tables.core.FloatAttribute()



class Interaction(EntityInteractionOrProperty):
    """ Represents an interaction

        Attributes:
            position (:obj:`int`): position at which interaction occurs
            score (:obj:`float`): ranking of the response of the interaction
    """
    #TODO: Assess difference between score and confidence
    position = obj_tables.core.IntegerAttribute()
    score = obj_tables.core.FloatAttribute(default=0)
    confidence = obj_tables.core.StringAttribute()

class SpecieInteraction(Interaction):
    """ Represents a protein interaction

    Attributes:
        specie_a (:obj:`str`):
        specie_b (:obj:`str`):
        stoichiometry_a (:obj:`int`):
        stoichiometry_b (:obj:`int`):
        loc_a (:obj:`str`):
        loc_b (:obj:`str`):
    """

    specie_a = obj_tables.core.OneToOneAttribute('Specie', related_name='specie_interaction')
    specie_b = obj_tables.core.OneToOneAttribute('Specie', related_name='specie_interaction')
    stoichiometry_a = obj_tables.core.IntegerAttribute()
    stoichiometry_b = obj_tables.core.IntegerAttribute()
    loc_a = obj_tables.core.StringAttribute()
    loc_b = obj_tables.core.StringAttribute()
    type_a = obj_tables.core.StringAttribute()
    type_b = obj_tables.core.StringAttribute()
    interaction_type = obj_tables.core.StringAttribute()



class Reaction(Interaction):
    """ Represents a reaction

    Attributes:
        participants (:obj:`list` of :obj:`ReactionParticipant`): list of participants
        reversible (:obj:`bool`): :obj:`True` if the reaction is reversible
    """
    kinetic_law_id = obj_tables.core.IntegerAttribute()
    participants = obj_tables.core.ManyToManyAttribute('ReactionParticipant', related_name='reactions')
    reversible = obj_tables.core.BooleanAttribute()

    def get_reactants(self):
        """ Get the reactants

        Returns:
            :obj:`list` of :obj:`ReactionParticipant`: list of reactants
        """
        return list(filter(lambda p: p.coefficient < 0, self.participants))

    def get_products(self):
        """ Get the products

        Returns:
            :obj:`list` of :obj:`ReactionParticipant`: list of products
        """
        return list(filter(lambda p: p.coefficient > 0, self.participants))

    def get_modifiers(self):
        """ Get the modifiers

        Returns:
            :obj:`list` of :obj:`ReactionParticipant`: list of modifiers
        """
        return list(filter(lambda p: p.coefficient == 0, self.participants))

    def get_ordered_participants(self, collapse_repeated=True):
        """ Get an ordered list of the participants

        Args:
            collapse_repeated (:obj:`bool`): if :obj:`True`, collapse any repeated participants

        Returns:
            :obj:`list` of :obj:`ReactionParticipant`: ordered list of reaction participants
        """

        # create copy of participants
        participants = copy.copy(list(self.participants))

        # collapse repeated participants
        if collapse_repeated:
            i_part = 0
            while i_part < len(participants):
                part = participants[i_part]
                if part.coefficient == 0:
                    participants.pop(i_part)
                    continue

                i_other_part = i_part + 1
                while i_other_part < len(participants):
                    other_part = participants[i_other_part]
                    if other_part.specie.structure == part.specie.structure and \
                            other_part.specie.id == part.specie.id and \
                            ((other_part.compartment is None and part.compartment is None) or
                                other_part.compartment.id == part.compartment.id) and \
                            numpy.sign(other_part.coefficient) == numpy.sign(part.coefficient):
                        part.coefficient += other_part.coefficient
                        participants.pop(i_other_part)
                    else:
                        i_other_part += 1

                i_part += 1

        # order participants
        participants.sort(key=lambda part: part.order if part.order is not None else 1e10)

        # return
        return participants

    def stringify(self):
        #TODO: Add the modifier

        participants = self.get_ordered_participants()
        result = ''
        for item in participants:
            if item.coefficient == -1 and len(result) < 1:
                result += str(item.specie.name)
            elif item.coefficient == -1 and len(result) > 1:
                result += ' + ' + str(item.specie.name)
            elif item.coefficient == 1 and '-->' not in result:
                result += ' --> ' + str(item.specie.name)
            elif item.coefficient == 1:
                result += ' + ' + str(item.specie.name)
            else:
                continue

        return result

    def get_reactant_product_pairs(self):
        """ Get list of pairs of similar reactants and products

        Note: This requires the modeler to have ordered the reactans and products by their similarity. The modeler is required to
        specify this pairing because it cannot easily be computed. In particular, we have tried to use Tanitomo similarity to
        predict reactant-product pairings, but this doesn't adequately capture reaction centers.

        Returns:
            :obj:`list` of :obj:`tuple` of obj:`ReactionParticipant`, :obj:`ReactionParticipant`: list of pairs of similar reactants and products
        """
        participants = self.get_ordered_participants()
        reactants = list(filter(lambda p: p.coefficient < 0, participants))
        products = list(filter(lambda p: p.coefficient > 0, participants))
        return list(six.moves.zip_longest(reactants, products))

    def get_ec_numbers(self):
        """ Get the EC numbers from the list of cross references

        Returns:
            :obj:`list` of :obj:`str`: list of EC numbers
        """
        return list(filter(lambda xr: xr.namespace == 'ec-code', self.cross_references))

    def get_manual_ec_numbers(self):
        """ Get the manually assigned EC numbers from the list of cross references

        Returns:
            :obj:`list` of :obj:`str`: list of EC manually assigned numbers
        """
        return list(filter(lambda xr: xr.namespace == 'ec-code' and xr.assignment_method ==
                           ResourceAssignmentMethod.manual, self.cross_references))

    def get_predicted_ec_numbers(self):
        """ Get the predicted EC numbers from the list of cross references

        Returns:
            :obj:`list` of :obj:`str`: list of predicted EC numbers
        """
        return list(filter(lambda xr: xr.namespace == 'ec-code' and xr.assignment_method ==
                           ResourceAssignmentMethod.predicted, self.cross_references))

    def get_ec_number(self):
        """ Get the most relevant EC number from the list of cross references

        * If the reaction has a single manually-assigned EC number, return that
        * If the reaction has multiple manually-assigned EC numbers, return an error
        * Otherwise, return the most relevant predicted EC number

        Returns:
            :obj:`str`: most relevant EC number
        """

        # return
        manual_ec_numbers = list(filter(lambda xr: xr.namespace == 'ec-code' and xr.assignment_method ==
                                        ResourceAssignmentMethod.manual, self.cross_references))
        if len(manual_ec_numbers) == 1:
            return manual_ec_numbers[0].id
        elif len(manual_ec_numbers) > 1:
            raise ValueError('Reaction {} has multiple EC numbers: {}'.format(self.id, ', '.join(xr.id for xr in manual_ec_numbers)))

        # find most relevant predicted EC number
        max_relevance = float('-inf')
        max_id = ''
        for xr in self.cross_references:
            if xr.assignment_method == ResourceAssignmentMethod.predicted and xr.relevance > max_relevance:
                max_relevance = xr.relevance
                max_id = xr.id

        return max_id


class ReactionParticipant(obj_tables.core.Model):
    """ Represents a participant in a reaction

    Attributes:
        specie (:obj:`Specie`): molecular species
        compartment (:obj:`Compartment`): compartment
        coefficient (:obj:`float`): coefficient
        order (:obj:`int`): order in the list of participants
    """
    specie = obj_tables.core.ManyToOneAttribute(Specie, related_name='reaction_participants')
    compartment = obj_tables.core.ManyToOneAttribute(Compartment, related_name='reaction_participants')
    coefficient = obj_tables.core.FloatAttribute()
    order = obj_tables.core.IntegerAttribute()


class Genetics(obj_tables.core.Model):
    """ Represents a taxon

    Attributes:
        taxon (:obj:`str`): taxon name
        variation (:obj:`str`): the genetic variation from the wildtype taxon

        observations (:obj:`list` of :obj:`Observation`): list of observations
    """
    taxon = obj_tables.core.StringAttribute()
    variation = obj_tables.core.StringAttribute()

    def is_wildtype(self):
        """ Determine if the taxon is the wildtype taxon

        Returns:
            obj:`bool`: `True` if the taxon doesn't have any genetic perturbation(s)
        """
        return not self.variation

    def is_variant(self):
        """ Determine if the taxon is the wildtype stain

        Returns:
            :obj:`bool`: `True` if the taxon has at least one genetic perturbation
        """
        return not self.is_wildtype()


class Environment(obj_tables.core.Model):
    """ Represents the environment (temperature, pH, media chemical composition) of an observation

    Attributes:
        temperature (:obj:`float`): temperature in Celcius
        ph (:obj:`float`): pH
        media (:obj:`str`): the chemical composition of the environment

        observations (:obj:`list` of :obj:`Observation`): list of observations
    """
    temperature = obj_tables.core.FloatAttribute()
    ph = obj_tables.core.FloatAttribute()
    media = obj_tables.core.LongStringAttribute()
    growth_status = obj_tables.core.LongStringAttribute()
    growth_system = obj_tables.core.LongStringAttribute()


class Method(obj_tables.core.Model):
    """ Represents a method used to generate an observation

    Attributes:
        name (:obj:`str`): name
        description (:obj:`str`): description

        observations (:obj:`list` of :obj:`Observation`): list of observations
    """
    name = obj_tables.core.StringAttribute()
    description = obj_tables.core.LongStringAttribute()
    performer = obj_tables.core.StringAttribute()
    hardware = obj_tables.core.StringAttribute()
    software = obj_tables.core.StringAttribute()

class Synonym(obj_tables.core.Model):
    """
    Represents a synonym of a given physical entity or property

    Attributes:
        name (:obj:`str`): Name of the Synonym

    """
    name = obj_tables.core.StringAttribute()


class ExperimentalMethod(Method):
    """ Represents a experimental method used to generate an observation """
    pass


class ComputationalMethod(Method):
    """ Represents a computational method used to generate an observation

    Attributes:
        version (:obj:`str`): version
        arguments (:obj:`str`): string representation of the arguments to the method
    """
    version = obj_tables.core.StringAttribute()
    arguments = obj_tables.core.LongStringAttribute()


class ResourceAssignmentMethod(enum.Enum):
    """ Represents the method used to assign a cross reference to an observable """
    manual = 0
    predicted = 1


class Resource(obj_tables.core.Model):
    """ Represents an object in an external resource

    Attributes:
        namespace (:obj:`str`): namespace of the identifier (e.g. pubmed)
        id (:obj:`str`): identifier within :obj:`namespace` (e.g. PMID)
        relevance (:obj:`float`): numerical indicator relevance of the external resource to the observable
        assignment_method (:obj:`ResourceAssignmentMethod`): method used to assign the cross reference to the observable
    """
    namespace = obj_tables.core.StringAttribute()
    id = obj_tables.core.StringAttribute()
    relevance = obj_tables.core.FloatAttribute()
    assignment_method = obj_tables.core.EnumAttribute(ResourceAssignmentMethod)


class Reference(obj_tables.core.Model):
    """ Represent a reference for an observation

    Attributes:
        title (:obj:`str`): title
        editor (:obj:`str`): editor
        author (:obj:`str`): author
        year (:obj:`int`): year
        publication (:obj:`str`): publication
        volume (:obj:`str`): volume
        number (:obj:`str`): number
        chapter (:obj:`str`): chapter
        pages (:obj:`str`): pages
        url (:obj:`str`): url

        observations (:obj:`list` of :obj:`Observation`): list of observations
    """
    title = obj_tables.core.StringAttribute()
    editor = obj_tables.core.StringAttribute()
    author = obj_tables.core.StringAttribute()
    year = obj_tables.core.IntegerAttribute()
    publication = obj_tables.core.StringAttribute()
    volume = obj_tables.core.StringAttribute()
    number = obj_tables.core.StringAttribute()
    chapter = obj_tables.core.StringAttribute()
    pages = obj_tables.core.StringAttribute()
    url = obj_tables.core.StringAttribute()
