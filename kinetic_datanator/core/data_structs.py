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
        source (:obj:`str`): source
        id (:obj:`str`): identifier within the source
        relevance (:obj:`float`): numeric score of the relevance of the cross reference
        assignment_method (:obj:`CrossReferenceAssignmentMethod`): method by which the cross reference was assigned
    """

    def __init__(self, source='', id='', relevance=float('nan'), assignment_method=None):
        """
        Args:
            source (:obj:`str`, optional): source
            id (:obj:`str`, optional): identifier within the source
            relevance (:obj:`float`, optional): numeric score of the relevance of the cross reference
            assignment_method (:obj:`CrossReferenceAssignmentMethod`, optional): method by which the cross reference was assigned
        """
        self.source = source
        self.id = id
        self.relevance = relevance
        self.assignment_method = assignment_method
