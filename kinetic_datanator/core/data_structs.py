"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from enum import Enum
from kinetic_datanator.core import observation


class Consensus(object):
    """ Represents a consensus of one or more observations of an attribute of a component of a model

    Attributes:
        component (:obj:`str`): biological component that was observed
        attribute (:obj:`str`): attribute of the biological component that was observed
        value (:obj:`float`): consensus value of the attribute of the model component
        error (:obj:`float`): uncertainty of the value of the attribute of the model component
        units (:obj:`str`): units of the value of the attribute of the model component
        method (:obj:`str`): method used to calculate the consensus value and error
        observations (:obj:`list` of :obj:`observation.Observation`)
    """

    def __init__(self, component=None, attribute=None, value=None, error=None, units=None, method=None, observations=None):
        """
        Args:
            component (:obj:`str`, optional): biological component that was observed
            attribute (:obj:`str`, optional): attribute of the biological component that was observed
            value (:obj:`float`, optional): consensus value of the attribute of the model component
            error (:obj:`float`, optional): uncertainty of the value of the attribute of the model component
            units (:obj:`str`, optional): units of the value of the attribute of the model component
            method (:obj:`str`, optional): method used to calculate the consensus value and error
            observations (:obj:`list` of :obj:`observation.Observation`, optional): list of observations which the consensus
                value is based on
        """
        self.component = component
        self.attribute = attribute
        self.value = value
        self.error = error
        self.units = units
        self.method = method
        self.observations = observations or []


class Compound(object):

    def __init__(self, id, inchi_smiles="", sabioNames=[]):
        self.id = id
        self.inchi_smiles = inchi_smiles
        self.sabioNames = sabioNames


class CrossReferenceAssignmentMethod(Enum):
    """ Represents the method by which a cross reference was assigned """
    manual = 0
    predicted = 1


class CrossReference(object):
    """ Represent a cross reference

    Attributes:
        namespace (:obj:`str`): namespace
        id (:obj:`str`): identifier within the namespace
        relevance (:obj:`float`): numeric score of the relevance of the cross reference
        assignment_method (:obj:`CrossReferenceAssignmentMethod`): method by which the cross reference was assigned
    """

    def __init__(self, namespace='', id='', relevance=float('nan'), assignment_method=None):
        """
        Args:
            namespace (:obj:`str`, optional): namespace
            id (:obj:`str`, optional): identifier within the namespace
            relevance (:obj:`float`, optional): numeric score of the relevance of the cross reference
            assignment_method (:obj:`CrossReferenceAssignmentMethod`, optional): method by which the cross reference was assigned
        """
        self.namespace = namespace
        self.id = id
        self.relevance = relevance
        self.assignment_method = assignment_method
