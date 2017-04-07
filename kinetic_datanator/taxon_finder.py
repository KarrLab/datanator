""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from ete3 import NCBITaxa
import logging

def get_taxonomic_distance(base_species, comparedSpecies):
	#Input a base species and comparedSpecies. Outputs an int of how many nodes apart they are on the NCBI tree
	#to work on: currently, if the species is not recognized, and the genus is the same as the base species, it wont be recognized

	ncbi = NCBITaxa()

	base_species = ncbi.get_name_translator([base_species])[base_species][0]
	try:
		comparedSpecies = ncbi.get_name_translator([comparedSpecies])[comparedSpecies][0]
	except KeyError:
			try:
				genus = comparedSpecies.split(" ")[0]
				comparedSpecies = ncbi.get_name_translator([genus])[genus][0]
			except:
				print("Unrecognized Species: " + comparedSpecies)
				logging.error("Unrecognized Species: ".format(comparedSpecies))
				return ""

	tree = ncbi.get_topology([base_species, comparedSpecies],intermediate_nodes=True)
	#tree = ncbi.get_topology([base_species, comparedSpecies],intermediate_nodes=True)

	A = tree&base_species
	C = tree&comparedSpecies
	node1 = tree.get_common_ancestor(A,C)
	tree.prune(node1, "preserve_branch_length")
	#format = "preserve_branch_length")
	distance = tree.get_tree_root().get_distance(A)
	return distance

def get_taxonomic_lineage(base_species):
	ncbi = NCBITaxa()
	base_species = ncbi.get_name_translator([base_species])[base_species][0]
	lineage = ncbi.get_lineage(base_species)
	names = ncbi.get_taxid_translator(lineage)
	chain =  [names[taxid] for taxid in lineage]
	i = len(chain)
	#print(i)
	new = []
	while i > 0:
		new.append(chain[i-1])
		i = i-1
	return new

def download_ncbi_database():
	ncbi = NCBITaxa()
	ncbi.get_descendant_taxa('Homo')

def update_ncbi_database():
	ncbi = NCBITaxa()
	ncbi.update_taxonomy_database()
