"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-05-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model, common_schema, models
from kinetic_datanator.api.lib.data_manager import BaseManager
from kinetic_datanator.util import molecule_util

class MetaboliteManager(BaseManager):

    """ Manages metabolite information for API """

    def __init__(self, cache_dirname = None):
        super(MetaboliteManager, self).__init__(cache_dirname=cache_dirname)

    def get_compound_by_id(self, id):
        return self.data_source.query(models.Compound).get(id)

    def get_observed_concentrations(self, compound):
        """ Find observed concentrations for the metabolite or similar metabolites

        Args:
            compound (:obj:`models.Compound`): compound to find data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observations
        """

        concentrations = self.get_concentration_by_structure(compound.structure._value_inchi, only_formula_and_connectivity=False)
        observed_values = []

        for c in concentrations:
            metadata = self.metadata_dump(c)

            observable = data_model.Observable(
                specie = self._port(compound),
                compartment = data_model.Compartment(name = c._metadata.cell_compartment[0].name)
            )

            observed_values.append(data_model.ObservedValue(
                metadata = metadata,
                observable = observable,
                value = c.value,
                error = c.error,
                units = c.units
            ))

        return observed_values

    def get_concentration_by_structure(self, inchi, only_formula_and_connectivity=True, select = models.Concentration):
        """
        Args:
            inchi (:obj:`str`): inchi structure to find concentrations

        Returns:
            :obj:`list`: List of models.Concentration Objects
        """

        q = self.data_source.session.query(select).join((models.Compound, select.compound)).\
            join((models.Structure, models.Compound.structure))

        if only_formula_and_connectivity:
            formula_and_connectivity = molecule_util.InchiMolecule(inchi).get_formula_and_connectivity()
            condition = models.Structure._structure_formula_connectivity == formula_and_connectivity
        else:
            condition = models.Structure._value_inchi == inchi

        return q.filter(condition).all()

    def _search(self, value):
        """
        Args:
            value (:obj:`str`): string to search table for

        Returns:
            :obj:`list` of `models.Compound`: List of found models.Compound objects
        """

        return models.Compound.query.search(value).all()


    def _port(self, compound):
        structure = compound.structure._value_inchi if compound.structure else None
        references = [data_model.Resource(namespace=item.namespace, id=item._id) for item in compound._metadata.resource]
        return data_model.Specie(id=compound.compound_id, name=compound.compound_name, structure = structure, cross_references = references)


metabolite_manager = MetaboliteManager()
