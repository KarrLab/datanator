from kinetic_datanator.core import data_model
from kinetic_datanator.core import data_query
from kinetic_datanator.core import common_schema
from sqlalchemy import or_
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

    def get_observable_interactions_and_complex(self,protein):

        interact = self.get_interaction_by_subunit(protein.uniprot_id).all()

        interaction = []
        index = 1
        for item in interact:
            resource = data_model.Resource(namespace = item._metadata.resource[0].namespace,
                id = item._metadata.resource[0]._id)
            interaction.append(data_model.ProteinInteraction(participant_a = item.participant_a,
            participant_b = item.participant_b, interaction_id = item.interaction,
            stoichiometry_a = item.stoich_a, stoichiometry_b = item.stoich_b,
            site_a=item.site_a, site_b=item.site_b, cross_references = [resource],
            name = 'Interaction '+ str(index)))
            index += 1

        complex_ = self.get_known_complex_by_subunit(protein.uniprot_id)

        plex = []
        for item in complex_:
            resource = data_model.Resource(namespace = item._metadata.resource[0].namespace,
                id = item._metadata.resource[0]._id)
            plex.append(data_model.KnownProteinComplex(name = item.complex_name,
            go_id = item.go_id, go_dsc = item. go_dsc, funcat_id = item.funcat_id,
            funcat_dsc = item.funcat_dsc, su_cmt = item.su_cmt, complex_cmt = item.complex_cmt,
            disease_cmt = item.disease_cmt, class_name = item.class_name,
            family_name = item.family_name, molecular_weight= item.molecular_weight,
            cross_references = [resource]))


        ##FIXME: Put into encompassing structure

        return interaction, plex



    def get_interaction_by_subunit(self, uniprot, select = common_schema.ProteinInteractions):
        q = self.data_source.session.query(select).filter(or_(select.participant_a == 'uniprotkb:'+uniprot,
            select.participant_b == 'uniprotkb:'+uniprot))
        return q


    def get_known_complex_by_subunit(self, uniprot, select = common_schema.ProteinComplex):
        q = self.data_source.session.query(select).join(common_schema.ProteinSubunit, select.protein_subunit)
        condition = common_schema.ProteinSubunit.uniprot_id == uniprot
        return q.filter(condition)


    def get_subunits_by_known_complex(self, complex_name, select = common_schema.ProteinSubunit):
        q = self.data_source.session.query(select).join(common_schema.ProteinComplex, select.proteincomplex)
        condition = common_schema.ProteinComplex.complex_name == complex_name
        return q.filter(condition)
