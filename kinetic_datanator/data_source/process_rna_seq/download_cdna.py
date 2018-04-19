from six.moves.urllib.request import urlretrieve
import os
import shutil

def run(ensembl_info, top_dir):
    """Downloads the CDNA for a given sample, and creates a kallisto index file. 
            The CDNA file is stored in a "CDNA" subdirectory within the top directory. 
            The kalliso index files are stored within "kallisto_index_files" subdirectory within the top directory

        Args:
            experiment(:obj:`array_express.Experiment`): the array express experiment
            top_dirname(:obj:`str`): the name of the directory where the overall data is being stored

        """
<<<<<<< HEAD
    download_cdna(ensembl_info, top_dir)
    process_cdna(ensembl_info, top_dir)


def download_cdna(ensembl_info, top_dir):
=======
>>>>>>> a6a94d8301672297d9b8fa4d2f015f34ea8033a3

    DIRNAME = "{}/CDNA_FILES".format(top_dir)
    if not os.path.isdir(DIRNAME):
        os.makedirs(DIRNAME)
    spec_name = ensembl_info.organism_strain
    file_name = "{}/{}.cdna.all.fa.gz".format(DIRNAME, spec_name)
    url = ensembl_info.url
    if not os.path.isfile(file_name):
        file = urlretrieve(url, '{}/{}.cdna.all.fa.gz'.format(top_dir, spec_name))
        shutil.move('{}/{}.cdna.all.fa.gz'.format(top_dir, spec_name), DIRNAME)
    os.chdir(top_dir)

def process_cdna(ensembl_info, top_dir):
    DIRNAME = "{}/CDNA_FILES".format(top_dir)
    file_name = "{}/{}.cdna.all.fa.gz".format(DIRNAME, ensembl_info.organism_strain)
    KALLISTO_DIR = "{}/kallisto_index_files".format(top_dir)
    if not os.path.isdir(KALLISTO_DIR):
        os.makedirs(KALLISTO_DIR)
<<<<<<< HEAD
    if not os.path.isfile("{}/{}.idx".format(KALLISTO_DIR, ensembl_info.organism_strain)):
        os.system("kallisto index -i {}.idx {}".format(ensembl_info.organism_strain, file_name))
        shutil.move("{}/{}.idx".format(top_dir, ensembl_info.organism_strain), KALLISTO_DIR)
=======
    if not os.path.isfile("{}/{}.idx".format(KALLISTO_DIR, spec_name)):
        os.system("kallisto index -i {}.idx {}".format(spec_name, file_name))
        shutil.move("{}/{}.idx".format(top_dir, spec_name), KALLISTO_DIR)
>>>>>>> a6a94d8301672297d9b8fa4d2f015f34ea8033a3
