"""
This codebase takes CORUM protein complexes database and formats it to an SQL database

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-08-13
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from sqlalchemy import exists, ForeignKey, Column, Integer, String, Numeric
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy.ext.declarative import declarative_base
from datanator.core import data_source
import csv
import ete3
import os
import warnings
import zipfile
from six import BytesIO


Base = declarative_base()

""" --------------------------- Table Definitions ----------------------------------"""


class Taxon(Base):
    """  Represents a species
    Attributes:
        ncbi_id      (:obj:`int`): NCBI id
        species_name (:obj:`str`): name and possibly genetic variant
    """

    __tablename__ = 'taxon'
    ncbi_id = Column(Integer, primary_key=True)
    swissprot_id = Column(String(255))


class Observation(Base):
    """  Represents an observation (entries in the original DB)
    Attributes:
        id            (:obj:`int`): internal ID for the observation entry
        cell line     (:obj:`str`): cell line (in whcih the measurement was done)
        pur_method    (:obj:`str`): purification method
        pubmed_id     (:obj:`str`): Pubmed ID of the associated publication
        taxon_ncbi_id (:obj:`str`): NCBI taxonomy id of the organism
    """

    __tablename__ = 'observation'
    id = Column(Integer, primary_key=True)
    #organism         = Column(String(255))
    cell_line = Column(String(255))
    pur_method = Column(String(255))
    pubmed_id = Column(Integer)

    taxon_ncbi_id = Column(String(255), ForeignKey('taxon.ncbi_id'))
    taxon = relationship('Taxon', backref=backref('observation'), foreign_keys=[taxon_ncbi_id])
    # FORMAT: foreign_table = relationship('foreign_class',backref=backref('self_table'),foreign_keys=[self_column])


class Complex(Base):
    """  Represents a protein complex
    Attributes:
        observation_id (:obj:`int`): ID of the  observation
        complex_id     (:obj:`int`): ID of the complex
        complex_name   (:obj:`str`): Complex name
        go_id          (:obj:`str`): GO funtinal annotation
        go_dsc         (:obj:`str`): Description of the annotation
        funcat_id      (:obj:`str`): FUNCAT functional annotation
        funcat_dsc     (:obj:`str`): Description of the annotation
        su_cmt         (:obj:`str`): Subunit comments
        complex_cmt    (:obj:`str`): Compex comments
        disease_cmt    (:obj:`str`): Disease comments
    """

    __tablename__ = 'complex'
    complex_id = Column(Integer, primary_key=True)
    complex_name = Column(String(255))
    go_id = Column(String(255))
    go_dsc = Column(String(255))
    funcat_id = Column(String(255))
    funcat_dsc = Column(String(255))
    su_cmt = Column(String(255))
    complex_cmt = Column(String(255))
    disease_cmt = Column(String(255))

    observation_id = Column(Integer, ForeignKey('observation.id'))
    observation = relationship('Observation', backref=backref('complex'), foreign_keys=[observation_id])


class Subunit(Base):
    """  Represents subunits of complexes
    Attributes:
        id           (:obj:`int`): Internal subunit ID
        complex_id   (:obj:`int`): ID of the complex to which the subunit belongs
        su_uniprot   (:obj:`int`): UNIPROT ID
        su_entrezs   (:obj:`int`): ENTREZS ID
        protein_name (:obj:`str`): Name of the protein
        gene_name    (:obj:`str`): Gene name
        gene_syn     (:obj:`str`): Synonyms of the gene name
    """

    __tablename__ = 'subunits'
    id = Column(Integer, primary_key=True)
    su_uniprot = Column(String(255))
    su_entrezs = Column(Integer)
    protein_name = Column(String(255))
    gene_name = Column(String(255))
    gene_syn = Column(String(255))

    complex_id = Column(Integer, ForeignKey('complex.complex_id'))
    complex = relationship('Complex', backref=backref('subunits'), foreign_keys=[complex_id])


""" -------------------- Extracting Content and adding to DB ------------------- """


class Corum(data_source.HttpDataSource):
    """ A local sqlite copy of the `CORUM database <https://mips.helmholtz-muenchen.de/corum/>`_ """

    base_model = Base
    ENDPOINT_DOMAINS = {
        'corum': 'http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip',
    }

    def load_content(self):
        """ Collect and parse all data from CORUM website and add to SQLite database """
        database_url = self.ENDPOINT_DOMAINS['corum']
        req = self.requests_session
        session = self.session

        # Extract All Files and save to current directory
        response = req.get(database_url)
        z = zipfile.ZipFile(BytesIO(response.content))
        z.extractall(self.cache_dirname)
        self.cwd = os.path.join(self.cache_dirname, 'allComplexes.txt')

        # create object to find NCBI taxonomy IDs
        ncbi_taxa = ete3.NCBITaxa()

        with open(self.cwd, 'r') as file:
            i_entry = 0
            for entry in csv.DictReader(file, delimiter='\t'):
                # entry/line number in file
                i_entry += 1

                # stop if the maximum desired number of entries has been reached
                if i_entry > self.max_entries:
                    break

                # replace 'None' strings with None
                for key, val in entry.items():
                    if val == 'None':
                        entry[key] = None

                # extract attributes
                complex_id = int(entry['ComplexID'])
                complex_name = entry['ComplexName']
                cell_line = entry['Cell line']
                su_uniprot = entry['subunits(UniProt IDs)']  # SETS OF STRING IDS SEPARATED BY ;\
                su_entrez = entry['subunits(Entrez IDs)']  # SETS OF INT IDS SEPARATED BY ;
                pur_method = entry['Protein complex purification method']
                go_id = entry['GO ID']  # SETS OF INT IDS SEPARATED BY ; eg. GO:0005634
                go_dsc = entry['GO description']
                funcat_id = entry['FunCat ID']
                funcat_dsc = entry['FunCat description']
                pubmed_id = int(entry['PubMed ID'])
                protein_name = entry['subunits(Protein name)']
                gene_name = entry['subunits(Gene name)']
                gene_syn = entry['Synonyms']
                disease_cmt = entry['Disease comment']
                su_cmt = entry['Subunits comment']
                complex_cmt = entry['Complex comment']
                swissprot_id = entry['SWISSPROT organism']

                """ ----------------- Apply field level corrections-----------------"""
                # Split the semicolon-separated lists of subunits into protein components,
                # ignoring semicolons inside square brackets
                su_uniprot_list = parse_list(su_uniprot)
                su_entrez_list = parse_list(su_entrez)
                protein_name_list = parse_list(correct_protein_name_list(protein_name))

                # check list lengths match
                if len(protein_name_list) != len(su_entrez_list):
                    msg = 'Unequal number of uniprot/entrez subunits at line {}\n  {}\n  {}'.format(
                        i_entry, '; '.join(protein_name_list), '; '.join(su_entrez_list))
                    raise Exception(msg)

                if len(su_uniprot_list) != len(su_entrez_list):
                    msg = 'Unequal number of uniprot/entrezs subunits at line {}\n  {}\n  {}'.format(
                        i_entry, '; '.join(su_uniprot_list), '; '.join(su_entrez_list))
                    raise Exception(msg)

                # Fix the redundancy issue with swissprot_id field
                if swissprot_id:
                    swissprot_id, _, _ = swissprot_id.partition(';')
                    ncbi_name, _, _ = swissprot_id.partition(' (')
                    result = ncbi_taxa.get_name_translator([ncbi_name])
                    ncbi_id = result[ncbi_name][0]
                else:
                    ncbi_id = None

                """ ----------------- Export the entries to the SQLite database ----------------- """

                if ncbi_id:
                    q = session.query(Taxon).filter(Taxon.ncbi_id == ncbi_id)
                    if session.query(q.exists()).scalar():
                        taxon = q.first()
                    else:
                        taxon = Taxon(ncbi_id=ncbi_id, swissprot_id=swissprot_id)
                        session.add(taxon)
                else:
                    taxon = None

                observation = Observation(cell_line=cell_line, pur_method=pur_method, pubmed_id=pubmed_id, taxon=taxon)
                session.add(observation)

                complex = Complex(complex_id=complex_id, complex_name=complex_name, go_id=go_id, go_dsc=go_dsc, funcat_id=funcat_id,
                                  funcat_dsc=funcat_dsc, su_cmt=su_cmt, complex_cmt=complex_cmt, disease_cmt=disease_cmt,
                                  observation=observation)
                session.add(complex)

                for su_uniprot, su_entrez, protein_name in zip(su_uniprot_list, su_entrez_list, protein_name_list):
                    subunit = Subunit(su_uniprot=su_uniprot, su_entrezs=su_entrez, protein_name=protein_name,
                                      gene_name=gene_name, gene_syn=gene_syn, complex=complex)
                    session.add(subunit)

        session.commit()


def parse_list(str_lst):
    """ Parse a semicolon-separated list of strings into a list, ignoring semicolons that are inside square brackets

    Args:
        str_lst (:obj:`str`): semicolon-separated encoding of a list

    Returns:
        :obj:`list` of :obj:`str`: list
    """
    if str_lst:
        lst = []
        depth = 0
        phrase = ''
        for char in str_lst:
            if char == ';' and depth == 0:
                lst.append(phrase)
                phrase = ''
            else:
                if char == '[':
                    depth += 1
                elif char == ']':
                    depth -= 1
                phrase += char
        lst.append(phrase)
        return lst
    else:
        return [None]


def correct_protein_name_list(lst):
    """ Correct a list of protein names with incorrect separators involving '[Cleaved into: ...]'

    Args:
        lst (:obj:`str`): list of protein names with incorrect separators

    Returns:
        :obj:`str`: corrected list of protein names
    """
    if lst:
        lst = lst.replace('[Cleaved into: Nuclear pore complex protein Nup98;', '[Cleaved into: Nuclear pore complex protein Nup98];')
        lst = lst.replace('[Cleaved into: Lamin-A/C;', '[Cleaved into: Lamin-A/C];')
        lst = lst.replace('[Cleaved into: Lamin-A/C ;', '[Cleaved into: Lamin-A/C ];')
        lst = lst.replace('[Includes: Maltase ;', '[Includes: Maltase ];')
    return lst
