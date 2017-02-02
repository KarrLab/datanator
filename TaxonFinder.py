from ete3 import NCBITaxa
import argparse


#Input a base species and comparedSpecies. Outputs an int of how many nodes apart they are on the NCBI tree
def getTaxonomicDistance(baseSpecies, comparedSpecies):
	ncbi = NCBITaxa()

	baseSpecies = ncbi.get_name_translator([baseSpecies])[baseSpecies][0]
	comparedSpecies = ncbi.get_name_translator([comparedSpecies])[comparedSpecies][0]

	tree = ncbi.get_topology([baseSpecies, comparedSpecies],intermediate_nodes=True)

	A = tree&baseSpecies
	C = tree&comparedSpecies
	node1 = tree.get_common_ancestor(A,C)
	tree.prune(node1, "preserve_branch_length")
	#format = "preserve_branch_length")
	distance = tree.get_tree_root().get_distance(A)
	return distance

def getTaxonomy(baseSpecies):
	ncbi = NCBITaxa()
	baseSpecies = ncbi.get_name_translator([baseSpecies])[baseSpecies][0]
	#print baseSpecies
	lineage = ncbi.get_lineage(baseSpecies)
	names = ncbi.get_taxid_translator(lineage)
	chain =  [names[taxid] for taxid in lineage]
	i = len(chain)
	#print i
	new = []
	while i > 0:
		new.append(chain[i-1])
		i = i-1

	return new


if __name__ == '__main__':

	getTaxonomy('mycoplasma pneumoniae')