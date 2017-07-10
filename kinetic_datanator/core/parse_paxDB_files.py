"""
This codebase takes the txt files of the PaxDB protein abundance database
and inserts them into an SQL database

parse_paxDB.py - parses PaxDB txt files

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2017 June 3
:Copyright: 2017, Karr Lab
:License: MIT
"""

from define_tables  import Taxon, Dataset, Protein, Observation
from sqlalchemy     import exists

def parse_paxDB_files(session,report,file_id,data_files,data_folder):
    """ This function parses pax DB files and adds them to the SQL database
    Attributes:
        session     (:obj:)     : SQLalchemy object
        file_id     (:obj:`str`): internal ID of the file
        data_files  (:obj:`str`): list of the files to be processed
        data_folder (:obj:`str`): root folder of the database
    """
    file_path = data_files[file_id]

    # Get NCBI taxonomy ID from file name
    start  = file_path.find('/',len(data_folder)-1)+1
    finish = file_path.find('/',len(data_folder))
    ncbi_id = int(float(file_path[start:finish]))

    # Get file_name
    start     = file_path.find('/',len(data_folder))+1
    file_name = file_path[start:]
    #print(file_name)

    with open(file_path,'r') as f:
        lines=f.readlines()

        # Get species name
        start  = lines[0].find(':')+2
        finish = lines[0].find('-')-1
        species_name = lines[0][start:finish]

        field_name,_,_ = lines[0].partition(':')
        if field_name=='#name':
            pass
        else:
            print('Error found, see reports.txt')
            report.write('Warning: invalid #name field, excluding file form DB (file_id='+str(file_id)+'; '+file_name+')\n')
            return

        # Get score
        finish = len(lines[1])-1
        score = float(lines[1][8:finish])

        field_name,_,_ = lines[1].partition(':')
        if field_name=='#score':
            pass
        else:
            print('Error found, see reports.txt')
            report.write('Warning: invalid #score field, excluding file form DB (file_id='+str(file_id)+'; '+file_name+')\n')
            return

        # Get weight
        finish = lines[2].find('%')
        if finish==-1:
            weight = None
        else:
            weight = float(lines[2][9:finish])

        field_name,_,_ = lines[2].partition(':')
        if field_name=='#weight':
            pass
        else:
            print('Error found, see reports.txt')
            report.write('Warning: invalid #weight field, excluding file form DB (file_id='+str(file_id)+'; '+file_name+')\n')
            return

        # Get publication link
        start  = lines[3].find('http:')
        finish = lines[3].find('"',start)
        publication = lines[3][start:finish]

        field_name,_,_ = lines[3].partition(':')
        if field_name=='#description':
            pass
        else:
            print('Error found, see reports.txt')
            report.write('Warning: invalid #description field, excluding file form DB (file_id='+str(file_id)+'; '+file_name+')\n')
            return

        # Get organ
        start  = lines[4].find(':')+2
        finish = len(lines[4])-1
        organ  = lines[4][start:finish]

        field_name,_,_ = lines[4].partition(':')
        if field_name=='#organ':
            pass
        else:
            print('Error found, see reports.txt')
            report.write('Warning: invalid #organ field, excluding file form DB (file_id='+str(file_id)+'; '+file_name+')\n')
            return

        # Get coverage
        start    = lines[6].find(':')+2
        finish   = len(lines[6])-1
        coverage = float(lines[6][start:finish])

        field_name,_,_ = lines[6].partition(':')
        if field_name=='#coverage':
            pass
        else:
            print('Error found, see reports.txt')
            report.write('Warning: invalid #coverage field, excluding file form DB (file_id='+str(file_id)+'; '+file_name+')\n')
            return

        # Check column header
        column_headers = lines[10].split()
        if column_headers[0]=='#internal_id' and column_headers[1]=='string_external_id' and column_headers[2]=='abundance' and len(column_headers)<5:
            pass
        else:
            print('Error found, see reports.txt')
            report.write('Warning: invalid column headers, excluding file form DB (file_id='+str(file_id)+'; '+file_name+')\n')
            return

        """ --- Add taxon and database (metadata info) to session ---------- """

        q = session.query(Taxon).filter(Taxon.ncbi_id==ncbi_id)
        if session.query(q.exists()).scalar():
            taxon = q.first()
        else:
            taxon = Taxon(ncbi_id=ncbi_id, species_name=species_name)
            session.add(taxon)

        dataset = Dataset(publication=publication, file_name=file_name, score=score, weight=weight, coverage=coverage, taxon=taxon)
        session.add(dataset)
        #print(taxon.species_name)

        """ --- Parse individual measurements and add them to DB ----------- """

        for i in range(11,len(lines)):
            split_line = lines[i].split()
            protein_id = split_line[0]
            string_id  = split_line[1]
            abundance  = split_line[2]

            # Insert relevant table entries
            q = session.query(Protein).filter(Protein.protein_id==protein_id)
            if session.query(q.exists()).scalar():
                protein = q.first()
            else:
                protein = Protein(protein_id=protein_id, string_id=string_id)
                session.add(protein)

            observation = Observation(dataset=dataset, abundance=abundance, protein=protein)
            session.add(observation)

    return
