"""
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-20
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model, data_query, common_schema, models
from kinetic_datanator.api.lib.data_manager import BaseManager
from kinetic_datanator.util.constants import DATA_CACHE_DIR

class ProteinComplexManager(BaseManager):
    """ Manages protein complex information for API """

    def __init__(self, cache_dirname = DATA_CACHE_DIR):
        self.data_source = common_schema.CommonSchema(cache_dirname=cache_dirname)

    def get_complex_by_id(self, id):
        return self.data_source.session.query(models.ProteinComplex).get(id)

    def _search(self, value):
        return models.ProteinComplex.query.search(value).all()

    def _port(self, complex):
        resource = data_model.Resource(namespace = complex._metadata.resource[0].namespace,
            id = complex._metadata.resource[0]._id)
        return data_model.ProteinComplexSpecie(name = complex.complex_name,
        go_id = complex.go_id, go_dsc = complex. go_dsc, funcat_id = complex.funcat_id,
        funcat_dsc = complex.funcat_dsc, su_cmt = complex.su_cmt, complex_cmt = complex.complex_cmt,
        disease_cmt = complex.disease_cmt, class_name = complex.class_name,
        family_name = complex.family_name, molecular_weight= complex.molecular_weight,
        cross_references = [resource])

    def get_observable_complex(self, protein_subunit):
        """ Get known protein complex that were observed for a given subunit

        Args:
            protein_subunit (:obj:`models.ProteinSubunit`): subunit to find complex for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedSpecie`: list of Protein Complexes
        """
        observed_specie = []
        _complex = self.get_known_complex_by_subunit(protein_subunit.uniprot_id)

        for complex in _complex:
            metadata= self.metadata_dump(item)
            self._port(complex)

            observed_specie.append(data_model.ObservedSpecie(specie = complex, metadata=metadata))

        return observed_specie


complex_manager = ProteinComplexManager()
