"""
This codebase takes the corum protein comlexes database and inserts the
entries into an SQL database

parse_corum.py - parses PaxDB txt files

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2017 June 3
:Copyright: 2017, Karr Lab
:License: MIT

-------------------------------------------------------------------------------
Column headers:
 [0]: ComplexID      [1]: ComplexName         [2]: Organism          [3]: Synonyms
 [4]: Cell line      [5]: SUs(UniProt)        [6]: SUs(Entrezs)      [7]: Complex purification method
 [8]: GO ID          [9]: GO dsc         [10]: FunCat ID        [11]: FunCat dsc
[12]: PubMed ID     [13]: SUs(Protein name)  [14]: SUs(Gene name)   [15]: SUs(Gene name syn)
[16]: Disease cmt   [17]: SUs cmt            [18]: Complex cmt      [19]: SWISSPROT organism

SU - subunit; CMT - comment; dsc - description
"""

from define_corum_tables import Taxon, Observation, Complex, Subunit
from sqlalchemy          import exists
from sqlalchemy          import create_engine
from sqlalchemy.orm      import sessionmaker

def find_ncbi_id(swissprot_id):
    # an embarassing solution, but works for now
    if swissprot_id == 'Neovison vison (American mink) (Mustela vison)':
        ncbi_id = 452646
    elif swissprot_id == 'Cricetus cricetus (Black-bellied hamster)':
        ncbi_id = 10034
    elif swissprot_id == 'Bubalus bubalis (Domestic water buffalo)':
        ncbi_id = 89462
    elif swissprot_id == 'Macaca mulatta (Rhesus macaque)':
        ncbi_id = 9544
    elif swissprot_id == 'Oryctolagus cuniculus (Rabbit)':
        ncbi_id = 9986
    elif swissprot_id == 'Cavia porcellus (Guinea pig)':
        ncbi_id = 10141
    elif swissprot_id == 'Rattus norvegicus (Rat)':
        ncbi_id = 10116
    elif swissprot_id == 'Homo sapiens (Human)':
        ncbi_id = 9606
    elif swissprot_id == 'Mus musculus (Mouse)':
        ncbi_id = 10090
    elif swissprot_id == 'Bos taurus (Bovine)':
        ncbi_id = 9913
    elif swissprot_id == 'Sus scrofa (Pig)':
        ncbi_id = 9825 #Assume here it is Sus scrofa domestica
    elif swissprot_id == 'Rattus sp.':
        ncbi_id = 10118
    else:
        ncbi_id = 0
    return int(ncbi_id)

engine = create_engine('sqlite:///corum.db')
Session = sessionmaker(bind=engine)
session = Session()

with open('/home/balazs/Desktop/aggregator/data/corum_complexes.txt','r') as f:
    lines            = f.readlines()
    column_headers   = lines[0].split('\t')
    unique_organisms = []

    for i in range(1,len(lines)):
        observation = lines[i].split('\t')
        #for j in range(0,len(observation)):
        #    print('['+str(j)+']: '+observation[j])

        """ ----------------- Collect fields ----------------- """
        #organism     = observation[2] #not used, since it encodes the species, already done by swissprot_id
        #synonyms     = observation[3] not used since really not sure what it encodes, seems a not well controlled functional annotation

        complex_id   = int(observation[0])
        assert isinstance(complex_id,int)

        complex_name = observation[1]
        assert isinstance(complex_name,str)

        cell_line    = observation[4]
        assert isinstance(cell_line,str)

        su_uniprot   = observation[5] # SETS OF STRING IDS SEPARATED BY ;
        assert isinstance(su_uniprot,str)

        su_entrezs   = observation[6] # SETS OF INT IDS SEPARATED BY ;
        assert isinstance(su_entrezs,str)

        pur_method   = observation[7]
        assert isinstance(pur_method,str)

        go_id        = observation[8] # SETS OF INT IDS SEPARATED BY ; eg. GO:0005634
        assert isinstance(go_id,str)

        go_dsc       = observation[9]
        assert isinstance(go_dsc,str)

        funcat_id    = observation[10]
        assert isinstance(funcat_id,str)

        funcat_dsc   = observation[11]
        assert isinstance(funcat_dsc,str)

        pubmed_id    = int(observation[12])
        assert isinstance(pubmed_id,int)

        protein_name = observation[13]
        assert isinstance(protein_name,str)

        gene_name    = observation[14]
        assert isinstance(gene_name,str)

        gene_syn     = observation[15]
        assert isinstance(gene_syn,str)

        disease_cmt  = observation[16]
        assert isinstance(disease_cmt,str)

        su_cmt       = observation[17]
        assert isinstance(su_cmt,str)

        complex_cmt  = observation[18]
        assert isinstance(complex_cmt,str)

        swissprot_id = observation[19]
        assert isinstance(swissprot_id,str)

        #print(len(go_id.split(';')),len(funcat_id.split(';')),len(funcat_dsc.split(';')))
        #print(len(su_uniprot.split(';')),len(su_entrezs.split(';')),len(protein_name.split(';')),'\n')

        """ ----------------- Apply field level corrections-----------------"""
        # Pull apart the subunits to protein componetns
        su_uniprot   = su_uniprot.split(';')
        su_entrezs   = su_entrezs.split(';')
        protein_name = protein_name.split(';')

        if len(protein_name)!=len(su_entrezs):
            print('Unequal number of uniprot/entrezs subunits at line: ',i)
            continue
        if len(su_uniprot)!=len(su_entrezs):
            print('Unequal number of uniprot/entrezs subunits at line: ',i)
            continue

        n_subunits = len(su_entrezs)

        # Fix the redundancy issue with swissprot_id field
        swissprot_id = swissprot_id.split(';')
        swissprot_id = [swissprot_id[element].rstrip() for element in range(len(swissprot_id))]
        swissprot_id = list(set(swissprot_id))
        swissprot_id = swissprot_id[0]

        ncbi_id = find_ncbi_id(swissprot_id)
        if ncbi_id == 0:
            print('Incorrect NCBI taxonomy ID at line: ',i)
            continue


        """ ----------------- Export the entries to the DB -----------------"""

        q = session.query(Taxon).filter(Taxon.ncbi_id==ncbi_id)
        if session.query(q.exists()).scalar():
            taxon = q.first()
            #print('taxon exists')
        else:
            taxon = Taxon(ncbi_id = ncbi_id, swissprot_id = swissprot_id)
            session.add(taxon)
            #print('new taxon')

        observation = Observation(cell_line  = cell_line, pur_method = pur_method, pubmed_id = pubmed_id, taxon=taxon)
        session.add(observation)

        complex = Complex(complex_id = complex_id, complex_name = complex_name, go_id = go_id, go_dsc = go_dsc, funcat_id = funcat_id, \
                         funcat_dsc = funcat_dsc, su_cmt = su_cmt, complex_cmt = complex_cmt, disease_cmt = disease_cmt, observation=observation)
        session.add(complex)

        for j in range(0,n_subunits):
            subunit = Subunit(su_uniprot = su_uniprot[j], su_entrezs = su_entrezs[j], protein_name = protein_name[j], gene_name = gene_name, gene_syn = gene_syn, complex=complex)
            session.add(complex)


#flat_list = [item for element in unique_organisms for item in element]
#print(list(set(flat_list))) # This is the list of all unique Swiss prot fields

session.commit()
session.close()



"""
#new_complex = Complex(complex_id = complex_id, complex_name = complex_name, go_id = go_id, go_dsc = go_dsc, funcat_id = funcat_id, funcat_dsc = funcat_dsc, \
#              disease_cmt = disease_cmt, su_cmt = su_cmt, complex_cmt = complex_cmt)

#session.add(new_complex)


Other organisms in the Swissprot field
    'Chlorocebus aethiops (Green monkey) (Cercopithecus',
    'Cricetulus griseus (Chinese hamster) (Cricetulus b'
    'Chlorocebus sabaeus (Green monkey) (Cercopithecus',
    'Canis lupus familiaris (Dog) (Canis familiaris)',
    'None'
"""
"""
q = session.query(Taxon).filter(Taxon.ncbi_id==ncbi_id)
if session.query(q.exists()).scalar():
    taxon = q.first()
else:
    taxon = Taxon(ncbi_id=ncbi_id, species_name=species_name)
    session.add(taxon)

Print column headers
for i in range(0,1): #len(lines)):
    observation = lines[i].split('\t')
    for j in range(0,len(observation)):
        print('['+str(j)+']: '+observation[j])

"""
