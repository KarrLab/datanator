from kinetic_datanator.core import data_model, data_query, flask_common_schema, models
from sqlalchemy import or_

class ProteinInteractionandComplexQueryGenerator(data_query.CachedDataSourceQueryGenerator):
    """ Queries Proteins to find Interactions with other Proteins """

    def __init__(self,
                 taxon=None, max_taxon_dist=None, taxon_dist_scale=None, include_variants=False,
                 temperature=37., temperature_std=1.,
                 ph=7.5, ph_std=0.3, cache_dirname = None):
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
        super(ProteinInteractionandComplexQueryGenerator, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=flask_common_schema.FlaskCommonSchema(cache_dirname=cache_dirname))

    def get_observed_values(self):
        pass

    def get_observable_interactions(self, protein_subunit):
        """ Get known protein interactions that were observed for a given subunit

        Args:
            protein_subunit (:obj:`models.ProteinSubunit`): subunit to find interactions for

        Returns:
            :obj:`list` of :obj:`data_model.ProteinInteraction`: list of Protein Interactions
        """

        interact = self.get_interaction_by_subunit(protein_subunit.uniprot_id).all()

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

        return interaction

    def get_observable_complex(self, protein_subunit):
        """ Get known protein complex that were observed for a given subunit

        Args:
            protein_subunit (:obj:`models.ProteinSubunit`): subunit to find complex for

        Returns:
            :obj:`list` of :obj:`data_model.KnownProteinComplex`: list of Protein Complexes
        """

        complex_ = self.get_known_complex_by_subunit(protein_subunit.uniprot_id)
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

        return plex

    def get_observable_subunits(self, protein_complex):
        """ Get known protein subunit that were observed for a given protein complex

        Args:
            protein_subunit (:obj:`models.ProteinComplex`): complex to find subunits for

        Returns:
            :obj:`list` of :obj:`data_model.ProteinSpecie`: list of Protein Subunits
        """

        subunits = self.get_subunits_by_known_complex(protein_complex.complex_name).all()
        ans = []
        for item in subunits:
            resource = data_model.Resource(namespace = item._metadata.resource[0].namespace,
                id = item._metadata.resource[0]._id)
            ans.append(data_model.ProteinSpecie(name = item.subunit_name, uniprot_id = item.uniprot_id,
                sequence = item.canonical_sequence, entrez_id = item.entrez_id,
                gene_name = item.gene_name, length = item.length, mass= item.mass,
                cross_references = [resource]))

        return ans


    def get_interaction_by_subunit(self, uniprot, select = models.ProteinInteractions):
        """ Get interactions that were observed for a given uniprot id

        Args:
            uniprot (:obj:`str`): uniprot id to search for
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`models.ProteinInteractions` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for protein interactions that contain the uniprot_id
        """
        q = self.data_source.session.query(select).filter(or_(select.participant_a == 'uniprotkb:'+uniprot,
            select.participant_b == 'uniprotkb:'+uniprot))
        return q


    def get_known_complex_by_subunit(self, uniprot, select = models.ProteinComplex):
        """ Get known complexes that were observed for a given uniprot id subunit

        Args:
            uniprot (:obj:`str`): uniprot id to search for
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`models.ProteinComplex` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for protein complexes that contain the uniprot_id
        """
        q = self.data_source.session.query(select).join(models.ProteinSubunit, select.protein_subunit)
        condition = models.ProteinSubunit.uniprot_id == uniprot
        return q.filter(condition)


    def get_subunits_by_known_complex(self, complex_name, select = models.ProteinSubunit):
        """ Get known protein subunits that were observed for a given complex

        Args:
            complex_name (:obj:`str`): complex to find subunits for
            select (:obj:`sqlalchemy.ext.declarative.api.DeclarativeMeta` or :obj:`sqlalchemy.orm.attributes.InstrumentedAttribute`, optional):
                :obj:`models.ProteinSubunit` or one of its columns

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for protein subunits that are within the given protein complex
        """
        q = self.data_source.session.query(select).join(models.ProteinComplex, select.proteincomplex)
        condition = models.ProteinComplex.complex_name == complex_name
        return q.filter(condition)
