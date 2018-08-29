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
        resource = [data_model.Resource(namespace = resource.namespace,
            id = resource._id) for resource in complex._metadata.resource]
        return data_model.ProteinComplexSpecie(name = complex.complex_name,
        go_id = complex.go_id, go_dsc = complex. go_dsc, funcat_id = complex.funcat_id,
        funcat_dsc = complex.funcat_dsc, su_cmt = complex.su_cmt, complex_cmt = complex.complex_cmt,
        disease_cmt = complex.disease_cmt, class_name = complex.class_name,
        family_name = complex.family_name, molecular_weight= complex.molecular_weight,
        cross_references = resource)

    def get_observable_subunits(self, protein_complex):
        """ Get known protein subunit that were observed for a given protein complex

        Args:
            protein_subunit (:obj:`models.ProteinComplex`): complex to find subunits for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedSpecie`: list of Protein Subunits
        """

        subunits = self.get_subunits_by_known_complex(protein_complex.complex_name).all()
        observed_specie = []
        for item in subunits:
            metadata = self.metadata_dump(item)
            resource = data_model.Resource(namespace = item._metadata.resource[0].namespace,
                id = item._metadata.resource[0]._id)
            specie = data_model.ProteinSpecie(name = item.subunit_name, uniprot_id = item.uniprot_id,
                sequence = item.canonical_sequence, entrez_id = item.entrez_id,
                gene_name = item.gene_name, length = item.length, mass= item.mass,
                cross_references = [resource])

            observed_specie.append(data_model.ObservedSpecie(specie = specie, metadata=metadata))

        return observed_specie

    def get_known_complex_by_subunit(self, uniprot, select = models.ProteinComplex):
        """ Get known complexes that were observed for a given uniprot id subunit

        Args:
            uniprot (:obj:`str`): uniprot id to search for

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for protein complexes that contain the uniprot_id
        """
        q = self.data_source.session.query(select).join(models.ProteinSubunit, select.protein_subunit)
        condition = models.ProteinSubunit.uniprot_id == uniprot
        return q.filter(condition)


complex_manager = ProteinComplexManager()
