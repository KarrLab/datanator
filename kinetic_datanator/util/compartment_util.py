""" Utilities for dealing with compartments

:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2017-04-13
:Copyright: 2017, Karr Lab
:License: MIT
"""


class Compartment(object):
    """ Represents a compartment

    Attributes:
        name (:obj:`str`): name
    """

    def __init__(self, name):
        """
        Args:
            name (:obj:`str`): name
        """
        self.name = name
