import numpy as np
from ete3 import NCBITaxa
import os


def get_genome_size(organism_name):
	path = os.path.dirname(os.path.realpath(__file__))

	domain = get_taxonomic_lineage(organism_name)[-3:-2][0]
	if domain == "Bacteria":
		file = open(os.path.join(path, 'number_of_prokaryote_genes.txt'))
	if domain == 'Eukaryota':
		file = open(os.path.join(path, 'number_of_eukaryote_genes.txt'))
	lines = file.readlines()
	lines = [line.split("	") for line in lines]
	total = []
	for line in lines:
		if line[0] == organism_name:
			if not line[12] == '-':
				total.append(int(line[12]))

	return np.average(total)

def get_taxonomic_lineage(baseSpecies):
	ncbi = NCBITaxa()
	baseSpecies = ncbi.get_name_translator([baseSpecies])[baseSpecies][0]
	lineage = ncbi.get_lineage(baseSpecies)
	names = ncbi.get_taxid_translator(lineage)
	chain =  [names[taxid] for taxid in lineage]
	i = len(chain)
	new = []
	while i > 0:
		new.append(chain[i-1])
		i = i-1
	return new

if __name__ == '__main__':
	#print(get_genome_size('Mycoplasma pneumoniae'))
	#print(get_genome_size('Mus musculus'))

	#print(get_genome_size('Homo sapiens'))

	print(get_genome_size('Zea mays'))


