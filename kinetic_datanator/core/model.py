"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from obj_model import core
from kinetic_datanator.core import observation


class Parameter(core.Model):
    """ Represents the consensus of a set of observations

    Attributes:
        component (:obj:`str`): model component
        attribute (:obj:`str`): attribute of the model component
        value (:obj:`float`): value of the attribute of the model component
        error (:obj:`float`): uncertainty of the value of the attribute of the model component
        units (:obj:`str`): SI units fo the value
        evidence (:obj:`list` of :obj:`Observation`): list of observations
        consensus_method (obj:`str`): method used to determine the consensus (e.g. mean, median, mode) value from the individual observations
    """
    component = core.StringAttribute()
    attribute = core.StringAttribute()
    value = core.FloatAttribute()
    error = core.FloatAttribute()
    units = core.StringAttribute()
    evidence = core.ManyToManyAttribute('kinetic_datanator.core.observation.Observation', related_name='parameters')
    consensus_method = core.StringAttribute()
