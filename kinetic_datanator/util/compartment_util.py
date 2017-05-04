""" Utilities for dealing with compartments

:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2017-04-13
:Copyright: 2017, Karr Lab
:License: MIT
"""


class Compartment(object):
    """ Represents a compartment

    Attributes:
        id (:obj:`str`): identifier
        name (:obj:`str`): name
        cross_references (:obj:`list` of :obj:`CrossReference`): list of cross references
    """

    def __init__(self, id='', name='', cross_references=None):
        """
        Args:
            id (:obj:`str`, optional): identifier
            name (:obj:`str`, optional): name
            cross_references (:obj:`list` of :obj:`CrossReference`, optional): list of cross references
        """
        self.id = id
        self.name = name
        self.cross_references = cross_references = []
