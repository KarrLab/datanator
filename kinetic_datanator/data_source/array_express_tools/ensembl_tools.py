from ete3 import NCBITaxa
import ftplib
import os
import socket


class StrainInfo(object):
    """ Represents information about an ensembl reference genome

    Attributes:
        organism_strain (:obj:`str`): the ensembl strain in the reference genome
        download_url (:obj:`str`): the url for that strain's refernce genome
        full_strain_specificity (:obj:`bool`): whether or not the strain mathces the full specifity
            provided in the arra express sample
    """

    def __init__(self, organism_strain, download_url, full_strain_specificity):
        self.organism_strain = organism_strain
        self.download_url = download_url
        self.full_strain_specificity = full_strain_specificity


def get_taxonomic_lineage(base_species):
    """ Get the lineage of a species

        Args:
            base_species (:obj:`bool`): a species (e.g. escherichia coli)

        Returns:
            :`list` of :obj:`str`: a list of strings corresponding to the layer of its taxonomy
    """

    ncbi = NCBITaxa()
    base_species = ncbi.get_name_translator([base_species])[base_species][0]
    lineage = ncbi.get_lineage(base_species)
    names = ncbi.get_taxid_translator(lineage)
    chain = [names[taxid] for taxid in lineage]
    i = len(chain)
    new = []
    while i > 0:
        new.append(chain[i-1])
        i = i-1
    return new


def format_org_name(name):
    """ 
    Format the name of an organism so normalize all species names

        Args:
            name (:obj:`bool`): the name of a spcies (e.g. escherichia coli str. k12)

        Returns:
            :obj:`str`: the normalized version of the strain name (e.g. escherichia coli k12)
    """

    name = name.replace("substr. ", "").replace("str. ", "").replace("subsp. ", "")
    name = name.replace("substr ", "").replace("str ", "").replace("subsp ", "")
    name = name.replace("_str", "").replace('_substr', "").replace("_subsp", "")
    return name.lower()


def get_strain_info(sample):
    """ 
    Get information about the refernce genome that should be used for a given sample

        Args:
            sample (:obj:`array_express.Sample`): an RNA-Seq sample

        Returns:
            :obj:`EnsembleInfo`: Ensembl information about the reference genome
    """
    print(sample.experiment_id)
    organism = ""
    strain = ""
    url = ""
    spec_name = ""
    full_strain_specificity = True
    list_of_characteristics = [ch.category.lower() for ch in sample.characteristics]

    if list_of_characteristics.count('organism') == 1:
        for characteristic in sample.characteristics:
            if characteristic.category.lower() == 'organism':
                organism = characteristic.value
            if (characteristic.category.lower() == 'strain') or (characteristic.category.lower() == "strain background"):
                strain = characteristic.value
    else:
        raise LookupError("No organism single organism recorded for this sample")

    domain = get_taxonomic_lineage(organism)[-3:-2][0]
    if domain == "Bacteria":
        end_url = ""
        if strain:
            organism = "{} {}".format(organism.lower(), strain.lower())
        org_tree = organism.split(" ")
        for num in range(len(org_tree), 0, -1):
            if not(end_url) and num >= 2:
                if num < len(org_tree):
                    full_strain_specificity = False  # this means it didnt find the specificity on the first try
                file = open("{}/find_cdna_url.txt".format(os.path.dirname(os.path.abspath(__file__))))
                try_org = ""
                for word in org_tree[:num]:
                    try_org = try_org + word + " "
                try_org = try_org[:-1]
                for line in file.readlines():
                    sep = line.split("\t")
                    if format_org_name(sep[0].lower()) == format_org_name(try_org):
                        spec_name = format_org_name(sep[1])
                        end_url = (sep[12][:sep[12].find("collection")+10] + "/" + sep[1])
        url = "ftp://ftp.ensemblgenomes.org/pub/bacteria/current/fasta/{}/cdna/".format(end_url)

    if domain == 'Eukaryota':
        for name in organism.split(" "):
            if name[-1:] == ".":
                name = name[:-1]
            spec_name = spec_name + name + "_"
        spec_name = spec_name[:-1].lower().replace("-", "_")

        if get_taxonomic_lineage(organism)[-4:-3][0] != "Viridiplantae":
            url = "ftp://ftp.ensembl.org/pub/current_fasta/{}/cdna/".format(spec_name)
        elif get_taxonomic_lineage(organism)[-4:-3][0] == "Viridiplantae":
            url = "ftp://ftp.ensemblgenomes.org/pub/current/plants/fasta/{}/cdna/".format(spec_name)            
    return StrainInfo(spec_name, url, full_strain_specificity)


def get_ftp_url(url):
    ftp_url = None
    attempts = 0
    while attempts < 10:
        try:
            host = url[url.find("ftp",1):url.find("org")+3]
            directory = url[url.find("org")+3:]
            ftp = ftplib.FTP(host)
            ftp.login()
            ftp.cwd(directory)
            files = ftp.nlst()
            for file in files:
                if file[-14:] == "cdna.all.fa.gz":
                    ftp_url = "{}{}".format(url,file)
            return ftp_url
        except socket.error:
            attempts += 1
            if attempts==10:
                raise LookupError("FTP not responding")
    return ftp_url
