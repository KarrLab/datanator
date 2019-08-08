import datanator.collections.core
import json 
import requests
import requests.exceptions
import zipfile
import io


class metabolite_entry()
def connect():
    client, database, collection = datanator.collections.core.get_mongo_connection("Metabolites")

   
def add_entry(): 
    pass

def download_mdb(source):

    if source == 'ecmdb':
        url = "http://ecmdb.ca/download/ecmdb.json.zip"
    elif source == 'ymdb':
        url = 'http://ymdb.ca/system/downloads/current/ymdb.json.zip'
    response = requests.get(url)
    try: 
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        print(e)
    
    with zipfile.ZipFile(io.BytesIO(response.content)) as f:
        with f.open(source+".json", 'r') as json_file: 
            entries=json.load(json_file)

    return entries
def parse_ecmdb():
    entries=download_mdb(source='ecmdb')
    print(entries[1])



if __name__ == "__main__":
    parse_ecmdb()