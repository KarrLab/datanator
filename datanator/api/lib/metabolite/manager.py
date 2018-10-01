"""
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-05-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.core import data_model, common_schema, models
from datanator.api.lib.data_manager import BaseManager
from datanator.util import molecule_util
from datanator.util.constants import DATA_CACHE_DIR

class MetaboliteManager(BaseManager):

    """ Manages metabolite information for API """

    def __init__(self, cache_dirname = DATA_CACHE_DIR):
        self.data_source = common_schema.CommonSchema(cache_dirname=cache_dirname)
        # super(MetaboliteManager, self).__init__(cache_dirname=cache_dirname)

    def get_metabolite_by_id(self, id):
        return self.data_source.session.query(models.Metabolite).get(id)

    def get_observed_concentrations(self, metabolite):
        """ Find observed concentrations for the metabolite or similar metabolites

        Args:
            metabolite (:obj:`models.Metabolite`): metabolite to find data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedValue`: list of relevant observations
        """

        concentrations = self.get_concentration_by_structure(metabolite.structure._value_inchi, only_formula_and_connectivity=False)
        observed_values = []

        for c in concentrations:
            metadata = self.metadata_dump(c)

            observable = data_model.Observable(
                specie = self._port(metabolite),
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

        q = self.data_source.session.query(select).join((models.Metabolite, select.metabolite)).\
            join((models.Structure, models.Metabolite.structure))

        if only_formula_and_connectivity:
            formula_and_connectivity = molecule_util.InchiMolecule(inchi).get_formula_and_connectivity()
            condition = models.Structure._structure_formula_connectivity == formula_and_connectivity
        else:
            condition = models.Structure._value_inchi == inchi

        return q.filter(condition).all()


    def get_metabolite_by_structure(self, inchi, only_formula_and_connectivity=False, select=models.Metabolite):
        """ Get metabolites with the same structure. Optionally, get metabolites which only have
        the same core empirical formula and core atom connecticity (i.e. same InChI formula
        and connectivity layers).

        Args:
            inchi (:obj:`str`): molecule structure in InChI format
            only_formula_and_connectivity (:obj:`bool`, optional): if :obj:`True`, get metabolites which only have
                the same core empirical formula and core atom connecticity. if :obj:`False`, get metabolites with the
                identical structure.

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for matching metabolites
        """
        q = self.data_source.session.query(select).join((models.Structure, models.Metabolite.structure))
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
            :obj:`list` of `models.Metabolite`: List of found models.Metabolite objects
        """

        return models.Metabolite.query.search(value).all()


    def _port(self, metabolite):
        structure = metabolite.structure._value_inchi if metabolite.structure else None
        references = [data_model.Resource(namespace=item.namespace, id=item._id) for item in metabolite._metadata.resource]
        return data_model.Specie(id=metabolite.metabolite_id, name=metabolite.metabolite_name, structure = structure, cross_references = references)


metabolite_manager = MetaboliteManager()
