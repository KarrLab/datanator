from kinetic_datanator.core import data_model
from kinetic_datanator.core import data_query
from kinetic_datanator.core import common_schema

class ProteintoProteinInteractionQueryGenerator(data_query.CachedDataSourceQueryGenerator):
    """ Queries Proteins to find Interactions with other Proteins """

    def __init__(self,
                 taxon=None, max_taxon_dist=None, taxon_dist_scale=None, include_variants=False,
                 temperature=37., temperature_std=1.,
                 ph=7.5, ph_std=0.3):
        """
        Args:
            taxon (:obj:`str`, optional): target taxon
            max_taxon_dist (:obj:`int`, optional): maximum taxonomic distance to include
            taxon_dist_scale (:obj:`float`, optional): The scale of the taxonomic distance scoring distribution.
                This determines how quickly the score falls to zero away from zero.
            include_variants (:obj:`bool`, optional): if :obj:`True`, also include observations from mutant taxa
            temperature (:obj:`float`, optional): desired temperature to search for
            temperature_std (:obj:`float`, optional): how much to penalize observations from other temperatures
            ph (:obj:`float`, optional): desired pH to search for
            ph_std (:obj:`float`, optional): how much to penalize observations from other pHs
        """
        super(ProteintoProteinInteractionQueryGenerator, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=common_schema.CommonSchema())

    def get_observed_values(self):
        pass

    def get_observed_interactions_and_complex(self,protein):
        pass
    # outputs datamodel interaction and complex

    def get_interaction_by_subunit(self, uniprot, select = common_schema.ProteinInteractions):
        # q = self.data_source.session.query(select).join(common_schema.ProteinInteractions, common_schema.AbundanceData.subunit)
        # condition = common_schema.ProteinSubunit.gene_name == gene_name
        # return q.filter(condition)
        pass

    def get_known_complex_by_subunit(self, uniprot, select = common_schema.ProteinComplex):
        q = self.data_source.session.query(select).join(common_schema.ProteinSubunit, select.protein_subunit)
        condition = common_schema.ProteinSubunit.uniprot_id == uniprot
        return q.filter(condition)


    def get_subunits_by_known_complex(self, complex, select = common_schema.ProteinSubunit):
        q = self.data_source.session.query(select).join(common_schema.ProteinSubunit, select.protein_subunit)
        condition = common_schema.ProteinSubunit.uniprot_id == uniprot
