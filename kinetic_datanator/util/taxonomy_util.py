"""
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2017-04-11
:Copyright: 2017, Karr Lab
:License: MIT
"""

from ete3 import NCBITaxa
import six

TAXON_RANKS = ('root', 'superkingdom', 'pylum', 'class', 'order', 'family', 'genus', 'species', 'strain')
# :obj:`tuple` of `str`: oredered list of taxonomic ranks


def setup_database(force_update=False):
    """ Setup a local sqllite copy of the NCBI Taxonomy database. If :obj:`force_update` is `False`, then
    only download the content from NCBI and build the sqllite database, if a local database doesn't already
    exist. If :obj:`force_update` is `True`, then always download the content from NCBI and rebuild the
    sqllite copy of the database.

    Args:
        force_update (:obj:`bool`, optional):

            * :obj:`False`: only download the content for the database and build a local sqllite database
                if a local sqllite copy of the database doesn't already exist
            * :obj:`True`: always download the content for the database from NCBI and rebuild a local sqllite
                database
    """
    ncbi_taxa = NCBITaxa()
    if force_update:
        # force downloading of latest content from NCBI and (re)building of local sqllite database
        ncbi_taxa.update_taxonomy_database()
    else:
        # run an operation on the local sqllite database to trigger NCBITaxa to setup a local sqllite
        # database if one doesn't already exist
        ncbi_taxa.get_descendant_taxa('Homo')


class Taxon(object):
    """ Represents a taxon such as a genus, species, or strain

    Attributes:
        name (:obj:`str`): name of the taxon
        id_of_nearest_ncbi_taxon (:obj:`int`): ID of the nearest parent taxon which is in the NCBI database
        distance_from_nearest_ncbi_taxon (:obj:`int`): distance from the taxon to its nearest parent which
            is in the NCBI database
        additional_name_beyond_nearest_ncbi_taxon (:obj:`str`): additional part of the taxon's beyond that of
            its nearest parent in the NCBI database
    """

    def __init__(self, id_or_name):
        """
        Args:
            id_or_name (:obj:`int` or :obj:`str`): id or name of the taxon
        """

        self.name = None
        self.id_of_nearest_ncbi_taxon = None
        self.distance_from_nearest_ncbi_taxon = None
        self.additional_name_beyond_nearest_ncbi_taxon = None

        ncbi_taxa = NCBITaxa()

        if isinstance(id_or_name, six.string_types):
            name = id_or_name
            rank_names = name.split(' ')
            for i_rank in range(len(rank_names)):
                partial_name = ' '.join(rank_names[0:len(rank_names) - i_rank])
                result = ncbi_taxa.get_name_translator([partial_name])
                if result:
                    self.id_of_nearest_ncbi_taxon = result[partial_name][0]
                    self.distance_from_nearest_ncbi_taxon = i_rank
                    self.additional_name_beyond_nearest_ncbi_taxon = ''.join(' ' + n for n in rank_names[len(rank_names) - i_rank:])
                    self.name = ncbi_taxa.translate_to_names([self.id_of_nearest_ncbi_taxon])[0] \
                        + self.additional_name_beyond_nearest_ncbi_taxon
                    return

            self.name = name

        else:
            id = id_or_name
            self.id_of_nearest_ncbi_taxon = id
            self.distance_from_nearest_ncbi_taxon = 0
            self.additional_name_beyond_nearest_ncbi_taxon = ''
            self.name = ncbi_taxa.translate_to_names([id])[0]
            if self.name == id:
                raise ValueError('The NCBI taxonomy database does not contain a taxon with id {}'.format(id))

    def get_ncbi_id(self):
        """ Get the ID of the taxon within the NCBI database

        Returns:
            :obj:`int` or :obj:`None`: ID of the taxon within the NCBI database or
                :obj:`None` if the taxon isn't in the NCBI database
        """
        if self.distance_from_nearest_ncbi_taxon == 0:
            return self.id_of_nearest_ncbi_taxon
        return None

    def get_parent_taxa(self):
        """ Get parent taxa

        Returns:
            :obj:`list` of :obj:`Taxon`: list of parent taxa
        """
        if self.id_of_nearest_ncbi_taxon is None:
            return None

        cls = self.__class__
        ncbi_taxa = NCBITaxa()
        lineage = [cls(id) for id in ncbi_taxa.get_lineage(self.id_of_nearest_ncbi_taxon)]

        if self.additional_name_beyond_nearest_ncbi_taxon:
            base_name = ncbi_taxa.translate_to_names([self.id_of_nearest_ncbi_taxon])[0]
            names = self.additional_name_beyond_nearest_ncbi_taxon[1:].split(' ')
            for i_rank, name, in enumerate(names):
                lineage.append(cls(base_name + ''.join(' ' + n for n in name[0:i_rank+1])))

        return lineage[0:-1]

    def get_rank(self):
        """ Get the rank of the taxon

        Returns:
            :obj:`str`: rank of the taxon
        """
        if self.distance_from_nearest_ncbi_taxon == 0:
            ncbi_taxa = NCBITaxa()
            rank = ncbi_taxa.get_rank([self.id_of_nearest_ncbi_taxon])[self.id_of_nearest_ncbi_taxon]
            if rank != 'no rank':
                return rank

        return None

    def get_common_ancestor(self, other):
        """ Get the lastest common ancestor of two taxa

        Args:
            other (:obj:`Taxon`): a second taxon

        Returns:
            :obj:`Taxon`: latest common ancestor
        """
        if self.id_of_nearest_ncbi_taxon is None:
            return id_of_nearest_ncbi_taxon

        ncbi_taxa = NCBITaxa()
        tree = ncbi_taxa.get_topology([self.id_of_nearest_ncbi_taxon, other.id_of_nearest_ncbi_taxon], intermediate_nodes=True)
        self_node = tree.search_nodes(name=str(self.id_of_nearest_ncbi_taxon))[0]
        other_node = tree.search_nodes(name=str(other.id_of_nearest_ncbi_taxon))[0]
        ancestor = tree.get_common_ancestor(self_node, other_node)
        cls = self.__class__
        return cls(float(ancestor.name))

    def get_distance_to_common_ancestor(self, other):
        """ Calculate the number of links in the NCBI taxonomic tree between two taxa and their latest common ancestor

        Note: This distances depends on the granularity of the lineage of the taxon. For example, there are only 7 links
        between most bacteria species and the Bacteria superkingdom. However, there are 28 links between the Homo sapiens 
        species and the Eukaryota superkingdom.

        Args:
            other (:obj:`Taxon`): a second taxon

        Returns:
             :obj:`int`: number of links between :obj:`self` and its latest common ancestor with :obj:`other` in the NCBI
                taxonomic tree
        """
        if self.id_of_nearest_ncbi_taxon is None:
            return id_of_nearest_ncbi_taxon

        ncbi_taxa = NCBITaxa()
        tree = ncbi_taxa.get_topology([self.id_of_nearest_ncbi_taxon, other.id_of_nearest_ncbi_taxon], intermediate_nodes=True)
        self_node = tree.search_nodes(name=str(self.id_of_nearest_ncbi_taxon))[0]
        other_node = tree.search_nodes(name=str(other.id_of_nearest_ncbi_taxon))[0]
        ancestor = tree.get_common_ancestor(self_node, other_node)
        return tree.get_distance(self_node, ancestor) + self.distance_from_nearest_ncbi_taxon

    def get_distance_to_root(self):
        """ Get the distance from the taxon to the root of the NCBI taxonomy tree

        Returns:
            :obj:`int`: distance from the taxon to the root
        """
        if self.id_of_nearest_ncbi_taxon is None:
            return id_of_nearest_ncbi_taxon

        ncbi_taxa = NCBITaxa()
        return len(ncbi_taxa.get_lineage(self.id_of_nearest_ncbi_taxon)) - 1 + self.distance_from_nearest_ncbi_taxon

    def get_max_distance_to_common_ancestor(self):
        """ Get the maximum distance from the taxon to a common ancestor with another taxon

        Returns:
            :obj:`int`: maximum distance from the taxon to a common ancestor with another taxon
        """
        return self.get_distance_to_root()
