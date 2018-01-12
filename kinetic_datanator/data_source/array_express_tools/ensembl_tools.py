import os
from ete3 import NCBITaxa
import ftplib #import FTP


class EnsembleInfo():
	def __init__(self, organism_strain, download_url, full_strain_specificity):
		self.organism_strain = organism_strain
		self.download_url = download_url
		self.full_strain_specificity = full_strain_specificity



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

def format_org_name(name):

	name = name.replace("substr. ", "").replace("str. ", "").replace("subsp. ", "")
	name = name.replace("substr ", "").replace("str ", "").replace("subsp ", "")
	name = name.replace("_str","").replace('_substr', "").replace("_subsp", "")
	return name.lower()



def get_ensembl_info(sample):

	organism = ""
	strain = ""
	url = ""
	full_strain_specificity = True
	for characteristic in sample.characteristics:
		if characteristic.category.lower() == 'organism':
			organism = characteristic.value
		if (characteristic.category.lower() == 'strain') or (characteristic.category.lower() =="strain background"):
			strain = characteristic.value

	spec_name = ""
	domain = get_taxonomic_lineage(organism)[-3:-2][0]
	if domain == "Bacteria":
		try:
			end_url = ""
			if strain:
				organism = "{} {}".format(organism.lower(), strain.lower())
			org_tree = organism.split(" ")
			for num in range(len(org_tree),0,-1):
				if not(end_url) and num>=2:
					if num<len(org_tree):
						full_strain_specificity = False #this means it didnt find the specificity on the first try
					file = open("{}/find_cdna_url.txt".format(os.path.dirname(os.path.abspath(__file__))))
					try_org = ""
					for word in org_tree[:num]:
						try_org = try_org + word + " "
					try_org = try_org[:-1]
					for line in file.readlines():
						sep = line.split("	")
						if format_org_name(sep[0].lower()) == format_org_name(try_org):
							print("GOT IT!")
							spec_name = format_org_name(sep[1])
							end_url = (sep[12][:sep[12].find("collection")+10]+ "/" + sep[1])
			start_url = "ftp://ftp.ensemblgenomes.org/pub/bacteria/current//fasta/{}/cdna/".format(end_url)
			ftp = ftplib.FTP("ftp.ensemblgenomes.org")				
			ftp.login()
			ftp.cwd("/pub/bacteria/current//fasta/{}/cdna/".format(end_url)) 
			files = ftp.nlst()
			for file in files:
				if file[-14:] == "cdna.all.fa.gz":
					url = start_url + file
		except ftplib.error_perm as resp:
			if str(resp) == "550 No files found":
				print("No files in this directory")
			else:
				raise

	if domain == 'Eukaryota':
		for name in organism.split(" "):
			if name[-1:] == ".":
				name = name[:-1]
			spec_name = spec_name + name + "_"
		spec_name = spec_name[:-1].lower().replace("-","_")
		if get_taxonomic_lineage(organism)[-4:-3][0] != "Viridiplantae":
			URL = "ftp://ftp.ensembl.org/pub/current_fasta"
			cwd = "/pub/current_fasta"
			ftp = ftplib.FTP("ftp.ensembl.org")
		elif get_taxonomic_lineage(organism)[-4:-3][0] == "Viridiplantae":
			URL = "ftp://ftp.ensemblgenomes.org/pub/current/plants/fasta"
			cwd = "/pub/current/plants/fasta"
			ftp = ftplib.FTP("ftp.ensemblgenomes.org")
		try:
			ftp.login()
			ftp.cwd("{}/{}/cdna/".format(cwd, spec_name))
			files = ftp.nlst()
			for file in files:
				if file[-14:] == "cdna.all.fa.gz":
					url = ("{}/{}/cdna/{}".format(URL, spec_name, file))
		except ftplib.error_perm as resp:
			if str(resp) == "550 No files found":
				print("No files in this directory")
			else:
				raise
	return EnsembleInfo(spec_name, url, full_strain_specificity)


