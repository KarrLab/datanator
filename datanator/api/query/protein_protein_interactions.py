"""
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-05-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from datanator.core import data_model, data_query, common_schema, models
from sqlalchemy import or_

class ProteinInteractionandComplexQuery(data_query.CachedDataSourceQueryGenerator):
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
        super(ProteinInteractionandComplexQuery, self).__init__(
            taxon=taxon, max_taxon_dist=max_taxon_dist, taxon_dist_scale=taxon_dist_scale, include_variants=include_variants,
            temperature=temperature, temperature_std=temperature_std,
            ph=ph, ph_std=ph_std,
            data_source=common_schema.CommonSchema(cache_dirname=cache_dirname))


    def get_observed_result(self, component):
        """ Get known observed interactions, complexes, or subunits for a given component

        Args:
            component (:obj:`models.ProteinSubunit`) or (:obj:`models.ProteinComplex`): subunit or complex to find interactions for

        Returns:
            :obj:`list` of :obj:`data_model.SpecieInteraction`: list of Protein Interactions
        """

        if component.__class__.__name__ == 'ProteinSubunit':
            return self.get_observable_interactions(component), self.get_observable_complex(component)
        elif component.__class__.__name__ == 'ProteinComplex':
            return self.get_observable_subunits(component)

    def get_observable_interactions(self, protein_subunit):
        """ Get known protein interactions that were observed for a given subunit

        Args:
            protein_subunit (:obj:`models.ProteinSubunit`): subunit to find interactions for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedInteraction`: list of Observed Protein Interactions

        """
        interact = self.get_interaction_by_subunit(protein_subunit.uniprot_id).all()
        observed_interaction= []

        for item in interact:

            participants = []
            for participant, protein in [(item.type_a, item.protein_a), (item.type_b, item.protein_b)]:
                if participant == None:
                    specie = None

                if participant == 'protein':
                    prot = self.data_source.session.query(models.ProteinSubunit).filter_by(uniprot_id = protein).first()
                    if prot:
                        specie =  data_model.ProteinSpecie(name=prot.uniprot_id, uniprot_id =prot.uniprot_id,
                            entrez_id=prot.entrez_id, gene_name=prot.gene_name, length=prot.length,
                            mass=prot.mass, sequence=prot.canonical_sequence)
                elif participant == 'peptide':
                    specie = data_model.PolymerSpecie(name=protein, sequence=protein)

                participants.append(specie)

            metadata = self.metadata_dump(item)
            resource = [data_model.Resource(namespace=source.namespace, id=source._id) for source in item._metadata.resource]
            interaction = data_model.SpecieInteraction(specie_a=participants[0], specie_b=participants[1], stoichiometry_a = item.stoich_a, stoichiometry_b = item.stoich_b,
            loc_a=item.loc_a, loc_b=item.loc_b, cross_references=resource, name=item.name, confidence=item.confidence, type_a = item.type_a, type_b=item.type_b, interaction_type=item.interaction_type)
            observed_interaction.append(data_model.ObservedInteraction(interaction=interaction, metadata=metadata))


        return observed_interaction


    def get_observable_complex(self, protein_subunit):
        """ Get known protein complex that were observed for a given subunit

        Args:
            protein_subunit (:obj:`models.ProteinSubunit`): subunit to find complex for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedSpecie`: list of Protein Complexes
        """
        observed_specie = []
        complex_ = self.get_known_complex_by_subunit(protein_subunit.uniprot_id)

        for item in complex_:
            metadata= self.metadata_dump(item)
            resource = data_model.Resource(namespace = item._metadata.resource[0].namespace,
                id = item._metadata.resource[0]._id)
            complex = data_model.ProteinComplexSpecie(name = item.complex_name,
            go_id = item.go_id, go_dsc = item. go_dsc, funcat_id = item.funcat_id,
            funcat_dsc = item.funcat_dsc, su_cmt = item.su_cmt, complex_cmt = item.complex_cmt,
            disease_cmt = item.disease_cmt, class_name = item.class_name,
            family_name = item.family_name, molecular_weight= item.molecular_weight,
            cross_references = [resource])

            observed_specie.append(data_model.ObservedSpecie(specie = complex, metadata=metadata))

        return observed_specie

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


    def get_interaction_by_subunit(self, uniprot, select = models.ProteinInteraction):
        """ Get interactions that were observed for a given uniprot id

        Args:
            uniprot (:obj:`str`): uniprot id to search for

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for protein interactions that contain the uniprot_id
        """
        return self.data_source.session.query(select).filter(or_(select.protein_a == uniprot,
                select.protein_b == uniprot))

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


    def get_subunits_by_known_complex(self, complex_name, select = models.ProteinSubunit):
        """ Get known protein subunits that were observed for a given complex

        Args:
            complex_name (:obj:`str`): complex to find subunits for

        Returns:
            :obj:`sqlalchemy.orm.query.Query`: query for protein subunits that are within the given protein complex
        """
        q = self.data_source.session.query(select).join(models.ProteinComplex, select.proteincomplex)
        condition = models.ProteinComplex.complex_name == complex_name
        return q.filter(condition)
