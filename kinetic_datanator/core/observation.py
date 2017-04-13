"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from obj_model import core
import kinetic_datanator.core.model


class Observation(core.Model):
    """ Represents an observation about a biological system

    Attributes:
        component (:obj:`str`): biological component that was observed
        attribute (:obj:`str`): attribute of the biological component that was observed
        value (:obj:`float`): observed value
        error (:obj:`float`): uncertainty of the observed value
        units (:obj:`units`): SI units of the observed value
        taxon (:obj:`Taxon`): the taxon (and any genetic perturbations) that the component was observed in
        environment (:obj:`Environment`): environment that the component was observed in
        method (:obj:`Method`): method that was used to make the observation
        reference (:obj:`Reference`): reference to the reference

        parameters (:obj:`list` of :obj:`Parameter`): list of parameters
    """
    component = core.ManyToOneAttribute('Component', related_name='observations')
    attribute = core.StringAttribute()
    value = core.FloatAttribute()
    error = core.FloatAttribute()
    units = core.StringAttribute()
    taxon = core.ManyToOneAttribute('Taxon', related_name='observations')
    environment = core.ManyToOneAttribute('Environment', related_name='observations')
    method = core.ManyToOneAttribute('Method', related_name='observations')
    reference = core.ManyToOneAttribute('Reference', related_name='observations')

class Component(core.Model):
    """ Represents a component of a biological system

    Attributes:
        id (:obj:`str`): identifier
        namespace (:obj:`str`): namespace of the identifier
    """
    id = core.StringAttribute()
    namespace = core.StringAttribute()


class MoleculeComponent(Component):
    """ Represents a molecule component of a biological system

    Attributes:
        structure (:obj:`str`): structure
    """
    structure = core.LongStringAttribute()


class Taxon(core.Model):
    """ Represents a taxon

    Attributes:
        name (obj:`str`): name
        perturbations (:obj:`str`): the genetic perturbations to the wildtype taxon

        observations (:obj:`list` of :obj:`Observation`): list of observations
    """
    name = core.StringAttribute()
    perturbations = core.StringAttribute()

    def is_wildtype(self):
        """ Determine if the taxon is the wildtype taxon

        Returns:
            obj:`bool`: `True` if the taxon doesn't have any genetic perturbation(s)
        """
        return not self.perturbations

    def is_mutant(self):
        """ Determine if the taxon is the wildtype stain

        Returns:
            :obj:`bool`: `True` if the taxon has at least one genetic perturbation
        """
        return not self.is_wildtype()


class Environment(core.Model):
    """ Represents the environment (temperature, pH, media chemical composition) of an observation

    Attributes:
        temperature (:obj:`float`): temperature in Celcius
        ph (:obj:`float`): pH
        media (:obj:`str`): the chemical composition of the environment

        observations (:obj:`list` of :obj:`Observation`): list of observations
    """
    temperature = core.FloatAttribute()
    ph = core.FloatAttribute()
    media = core.LongStringAttribute()


class Method(core.Model):
    """ Represents a method used to generate an observation

    Attributes:
        name (:obj:`str`): name

        observations (:obj:`list` of :obj:`Observation`): list of observations
    """
    name = core.StringAttribute()


class ExperimentalMethod(Method):
    """ Represents a experimental method used to generate an observation

    Attributes:
        description (:obj:`str`): description
    """
    description = core.LongStringAttribute()


class ComputationalMethod(Method):
    """ Represents a computational method used to generate an observation

    Attributes:
        version (:obj:`str`): version
        arguments (:obj:`str`): string representation of the arguments to the method
    """
    version = core.StringAttribute()
    arguments = core.LongStringAttribute()


class Reference(core.Model):
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
        uri (:obj:`str`): uri

        observations (:obj:`list` of :obj:`Observation`): list of observations
    """
    title = core.StringAttribute()
    editor = core.StringAttribute()
    author = core.StringAttribute()
    year = core.IntegerAttribute()
    publication = core.StringAttribute()
    volume = core.StringAttribute()
    number = core.StringAttribute()
    chapter = core.StringAttribute()
    pages = core.StringAttribute()
    uri = core.StringAttribute()
