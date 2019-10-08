""" Downloads and parses the IntAct database of protein-protein interactions
"""

import glob
import pandas
import os
import zipfile
from ftplib import FTP
from datanator.util import mongo_util
import json

class IntActNoSQL(mongo_util.MongoUtil):
    """ A local MongoDB copy of the IntAct database """

    def __init__(self, cache_dirname=None, MongoDB=None, db=None,
                 replicaSet=None, verbose=False, max_entries=float('inf'),
                 username = None, password = None, authSource = 'admin'):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        super(IntActNoSQL, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet, db=db,
                                              verbose=verbose, max_entries=max_entries, username = username,
                                              password = password, authSource = authSource)
        self.client_interaction, self.db_interaction, self.collection_interaction = self.con_db('intact_interaction')
        self.client_complex, self.db_complex, self.collection_complex = self.con_db('intact_complex')

    def load_content(self):
        """ Load the content of the local copy of the data source """

        # Download data from FTP Server
        self.download_content()

        # parse data and build MongoDB database
        self.add_complexes()
        self.add_interactions()

    def download_content(self):
        """ Download data from FTP server """
        if not os.path.exists(os.path.join(self.cache_dirname, 'intact')):
            os.makedirs(os.path.join(self.cache_dirname, 'intact'))
        if not os.path.exists(os.path.join(self.cache_dirname, 'intact', 'complextab')):
            os.makedirs(os.path.join(self.cache_dirname, 'intact', 'complextab'))
        if not os.path.exists(os.path.join(self.cache_dirname, 'intact', 'psimitab')):
            os.makedirs(os.path.join(self.cache_dirname, 'intact', 'psimitab'))

        ftp = FTP('ftp.ebi.ac.uk')
        ftp.login()

        ftp.cwd('/pub/databases/intact/complex/current/complextab/')
        rel_filenames = ftp.nlst()
        for rel_filename in rel_filenames:
            local_filename = os.path.join(self.cache_dirname, 'intact', 'complextab', rel_filename)
            if not os.path.exists(local_filename):
                with open(local_filename, 'wb') as file:
                    ftp.retrbinary('RETR ' + rel_filename, file.write)

        ftp.cwd('/pub/databases/intact/current/psimitab/')
        local_filename = os.path.join(self.cache_dirname, 'intact', 'psimitab', 'intact_negative.txt')
        with open(local_filename, 'wb') as file:
            ftp.retrbinary('RETR ' + 'intact_negative.txt', file.write)

        ftp.quit()

    def add_complexes(self):
        """ Parse complexes from data and add complexes to MongoDB """
        raw_columns = [
            '#Complex ac', 'Recommended name', 'Taxonomy identifier',
            'Identifiers (and stoichiometry) of molecules in complex', 'Experimental evidence',
            'Go Annotations', 'Description', 'Source',
        ]
        relabeled_columns = ['identifier', 'name', 'ncbi_id', 'subunits', 'evidence', 'go_annotation', 'go_description', 'source']

        filenames = glob.glob(os.path.join(self.cache_dirname, 'intact', 'complextab', '*.tsv'))
        total_operations = 0
        for filename in filenames:
            raw_data = pandas.read_csv(filename, delimiter='\t', encoding='utf-8')
            relabeled_data = raw_data.loc[:, raw_columns]
            relabeled_data.columns = relabeled_columns
            # relabeled_data = relabeled_data.set_index('identifier')

            relabeled_data_json = json.loads(relabeled_data.to_json(orient = 'records'))
            # separate string of subunits
            for j in range(len(relabeled_data_json)):
                if total_operations == self.max_entries:
                    break
                if self.verbose and total_operations%50 ==0:
                    print ('Inserting {} of {} complex document'.format(total_operations+1, min(self.max_entries,len(relabeled_data_json))))

                # separate string of subunits
                subunit_info = []
                if relabeled_data_json[j]['subunits'] is not None:
                    subunit_list_with_count = relabeled_data_json[j]['subunits'].split('|')
                    subunit_list = [item.split('(')[0] for item in subunit_list_with_count] 
                    count_list = [item.split('(')[1].split(')')[0] for item in subunit_list_with_count]
                    for unit, count in zip(subunit_list, count_list):
                        subunit_info.append( {'uniprot_id': unit, 'count': count} )
                    relabeled_data_json[j]['subunits'] = subunit_info

                # separate string of go_annotation
                annotation_list = []
                if relabeled_data_json[j]['go_annotation'] is not None:
                    go_anno_after_split = relabeled_data_json[j]['go_annotation'].split('|')
                    go_id = [item[3:10] for item in go_anno_after_split]
                    go_term = [item[11:-1] for item in go_anno_after_split]
                    for _id, term in zip(go_id, go_term):
                        annotation_list.append( {'go_id': _id, 'go_term': term} )
                    relabeled_data_json[j]['go_annotation'] = annotation_list

                self.collection_complex.replace_one({'identifier': relabeled_data_json[j]['identifier']},
                        relabeled_data_json[j],
                        upsert=True
                        )
                total_operations += 1 

    def add_interactions(self):
        """ Parse interactions from data and add interactions to mongodb database """
        data = pandas.read_csv(os.path.join(self.cache_dirname, 'intact', 'psimitab', 'intact_negative.txt'),
                               delimiter='\t', encoding='utf-8')
        for index, row in data.iterrows():
            if index == self.max_entries:
                break
            if self.verbose and index%20==0:
                print ('Inserting {} of {} intercation document'.format(index+1, min(self.max_entries,len(data.index) ) ))

            interaction = {} # one document
            interaction['_id'] = index
            interaction['protein_a'], interaction['gene_a'] = self.find_protein_gene(row['#ID(s) interactor A'], row['Alias(es) interactor A'])
            interaction['protein_b'], interaction['gene_b'] = self.find_protein_gene(row['ID(s) interactor B'], row['Alias(es) interactor B'])
            interaction['interaction_type'] = self.find_between_psi_mi_parentheses(row['Interaction type(s)'])
            interaction['method'] = self.find_between_psi_mi_parentheses(row['Interaction detection method(s)'])
            interaction['type_a'] = self.find_between_psi_mi_parentheses(row['Type(s) interactor A'])
            interaction['type_b'] = self.find_between_psi_mi_parentheses(row['Type(s) interactor B'])
            interaction['role_a'] = self.find_between_psi_mi_parentheses(row['Biological role(s) interactor A'])
            interaction['role_b'] = self.find_between_psi_mi_parentheses(row['Biological role(s) interactor B'])
            interaction['feature_a'] = row['Feature(s) interactor A']
            interaction['feature_b'] = row['Feature(s) interactor B']
            interaction['stoich_a'] = row['Stoichiometry(s) interactor A']
            interaction['stoich_b'] = row['Stoichiometry(s) interactor B']
            interaction['interaction_id'] = row['Interaction identifier(s)']
            interaction['publication'] = self.find_pubmed_id(row['Publication Identifier(s)'])
            interaction['publication_author'] = row['Publication 1st author(s)']
            interaction['confidence'] = row['Confidence value(s)']

            self.collection_interaction.replace_one({'_id': interaction['_id']},
                    interaction,
                    upsert=True
                    )

    def find_protein_gene(self, interactor, alias):
        """ Parse the protein and gene identifiers from key-value pairs of interactors and their aliases

        Args:
            interactor (:obj:`str`): key-value pairs of interactor
            alias (:obj:`str`): key-value pairs of the alias of the interactor

        Returns:
            :obj:`str`: protein identifier
            :obj:`str`: gene identifier
        """
        protein = None
        gene = None
        if 'uniprotkb' in interactor:
            protein = self.split_colon(interactor)[1]
        else:
            if 'display_short' in alias:
                protein = self.find_between(alias, 'psi-mi:', '(display_short)')
            else:
                protein = None

        for item in self.split_line(alias):
            if '(gene name)' in item:
                gene = self.find_between(item, 'uniprotkb:', '(gene name)')

        return protein, gene

    def find_pubmed_id(self, string):
        """ Parse PubMed identifier from annotated key-value pair of publication type-identifier

        Args:
            string (:obj:`str`): key-value pair of publication type-identifier

        Returns:
            :obj:`str`: PubMed identifier
        """
        for item in self.split_line(string):
            if 'pubmed:' in item:
                return self.split_colon(item)[1]
        return None

    def find_between_psi_mi_parentheses(self, string):
        """ Find the text between parentheses in values of psi-mi key-value pairs

        Args:
            string (:obj:`str`): string

        Returns:
            :obj:`str`: substring between the first occurrence of the substring :obj:`first` and the 
                last occurrence of the substring :obj:`last
        """
        if 'psi-mi:' in string:
            return self.find_between(string, '(', ')')
        else:
            return None

    def find_between(self, string, first, last):
        """ Get the substring between the first occurrence of the substring :obj:`first` and the 
        last occurrence of the substring :obj:`last`

        Args:
            string (:obj:`str`): string
            first (:obj:`str`): starting substring
            last (:obj:`str`): ending substring

        Returns:
            :obj:`str`: substring between the first occurrence of the substring :obj:`first` and the 
                last occurrence of the substring :obj:`last
        """
        return string[string.index(first) + len(first):string.index(last, string.index(first) + len(first))]

    def split_colon(self, string):
        """ Split a string into substrings separated by ':'

        Args:
            string (:obj:`str`): string

        Returns:
            :obj:`list`: substring separated by ':'
        """
        return string.split(':')

    def split_line(self, string):
        """ Split a string into substrings separated by '|'

        Args:
            string (:obj:`str`): string

        Returns:
            :obj:`list`: substring separated by '|'
        """
        return string.split('|')

