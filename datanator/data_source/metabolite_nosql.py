'''
	Generates metabolite (ECMDB, YMDB) documents for MongoDB
    Stores documents in MongoDB and stores documents as JSON

:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2019-04-02
:Copyright: 2019, Karr Lab
:License: MIT
'''


import io
import os
import json
import requests
import requests.exceptions
import warnings
import zipfile
import xmltodict
from pymongo import MongoClient


class MetaboliteNoSQL():
    '''Loads metabolite information into mongodb and output documents as JSON files for each metabolite
        Attribuites:
            source: source database e.g. 'ecmdb' 'ymdb'
            MongoDB: mongodb server address e.g. 'mongodb://localhost:27017/'
            max_entries: maximum number of documents to be processed
            output_direcotory: directory in which JSON files will be stored.
    '''
    def __init__(self,output_directory, source, MongoDB, db, verbose=True, max_entries = float('inf')):
        self.verbose = verbose
        self.source = source
        self.db = db
        self.max_entries = max_entries
        self.MongoDB = MongoDB
        self.output_directory = output_directory      

        if self.source == 'ecmdb':
            self.domain = 'http://ecmdb.ca'
            self.compound_index = self.domain + \
                '/download/ecmdb.json.zip'  # list of metabolites
            self.compound_url = self.domain + '/compounds/{}.xml'
            self.collection_dir = output_directory
        else:
            self.domain = 'http://ymdb.ca'
            self.compound_index = self.domain + \
                '/system/downloads/current/ymdb.json.zip'  # list of metabolites
            self.compound_url = self.domain + '/compounds/{}.xml'
            self.collection_dir = output_directory

    # make connections wth mongoDB
    def con_db(self):
        try:
            client = MongoClient(self.MongoDB, 400)  # 400ms max timeout
            client.server_info()
            db = client[self.db]
            collection = db[self.source]
            return collection
        except pymongo.errors.ConnectionFailure:
            return ('Server not available')

    '''Each compound's collection of imformation is written into a json file for mongodb

    '''
    def write_to_json(self):

        if self.verbose:
            print('Download list of all compounds: ...')

        response = requests.get(self.compound_index)
        response.raise_for_status()
        collection = self.con_db()

        if self.verbose:
            print('... Done!')
        if self.verbose:
            print('Unzipping and parsing compound list ...')

        with zipfile.ZipFile(io.BytesIO(response.content), 'r') as zip_file:
            json_name = self.source + '.json'
            with zip_file.open(json_name, 'r') as json_file:
                entries = json.load(json_file)

        if self.verbose:
            print('  found {} compounds'.format(len(entries)))

        if self.source == 'ecmdb':
            entries.sort(key=lambda e: e['m2m_id'])
        else:
            entries.sort(key=lambda e: e['ymdb_id'])

        if len(entries) > self.max_entries:
            entries = entries[0:self.max_entries]
        
        if self.verbose:
            print('Downloading {} compounds ...'.format(len(entries)))

        for i_entry, entry in enumerate(entries):

            if self.verbose and (i_entry % 10 == 0):
                print('  Downloading compound {} of {}'.format(
                    i_entry + 1, len(entries)))

            # actual xml file for each metabolite
            if self.source == 'ecmdb':
                response = requests.get(
                    self.compound_url.format(entry['m2m_id']))
                try:
                    response.raise_for_status()
                except requests.exceptions.HTTPError:
                    warnings.warn(
                        'Unable to download data for compound {}'.format(entry['m2m_id']))
                    continue
            else:
                response = requests.get(
                    self.compound_url.format(entry['ymdb_id']))
                try:
                    response.raise_for_status()
                except requests.exceptions.HTTPError:
                    warnings.warn(
                        'Unable to download data for compound {}'.format(entry['ymdb_id']))
                    continue

            doc = xmltodict.parse(response.text)

            # delete key "compound" but keep key's value
            new_doc = doc['compound']
            # original source spelled wikipedia wrong
            new_doc['wikipedia'] = new_doc.pop('wikipidia')

            os.makedirs(os.path.dirname(str(self.collection_dir) + '/'), exist_ok=True)
            if self.source == 'ecmdb':
                file_name = str(self.collection_dir) +'/'+ entry['m2m_id'] + '.json'
            else:
                file_name = str(self.collection_dir) + '/'+entry['ymdb_id'] + '.json'
            with open(file_name, "w") as f:
                f.write(json.dumps(new_doc, indent=4))

            
            if self.source == 'ecmdb':
                collection.replace_one(
                    {'m2m_id': new_doc['m2m_id']},
                    new_doc,
                    upsert=True
                )
            else:
                collection.replace_one(
                    {'ymdb_id': new_doc['ymdb_id']},
                    new_doc,
                    upsert=True
                )              
                

