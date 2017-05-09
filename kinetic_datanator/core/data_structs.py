from enum import Enum


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
