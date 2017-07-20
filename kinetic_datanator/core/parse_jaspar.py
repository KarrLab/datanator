# -*- coding: utf-8 -*-
"""
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017 July 12th
:Copyright: 2017, Karr Lab
:License: MIT
"""
from define_jaspar import MatrixObservation, BindingMatrix, Type, Family, Class, Resources, TranscriptionFactor, Species
from sqlalchemy.orm import sessionmaker
from csv import reader
from itertools import groupby
from operator import itemgetter
from urllib2 import urlopen

""" - - - - - - - - - - - - - Iterative Sorting and Data Type Functions Used for Parsing - - - - - - - - - - - - - - - """

def make_jaspar_int(data):
    data = list(data)
    for i in range(0,len(data)):
        data[i][0] = int(data[i][0])
    return data

def make_data_int(data):
    data = list(data)
    for i in range(0,len(data)):
        data[i][0] = int(data[i][0])
        data[i][2] = int(data[i][2])
        data[i][3] = float(data[i][3])
    return data

def sort(data):
    sortedlist = sorted(data,  key=itemgetter(0), reverse=False)
    return sortedlist

def group_by_jaspar_id(data):
    group = []
    key = []
    for k, g in groupby(data, lambda x: x[0]):
        group.append(list(g))
        key.append(k)
    return group, key

def group_by_position(data):
    group = []
    for k, g in groupby(data, lambda x: x[2]):
        group.append(list(g))
    return group

""" - - - - - - - - - - - - - -  Retrieval Functions - - - - - - - - - - - - - - - - - - """

def call_Attribute(key_list, data_list, jid):
    answer = []
    ind = 0
    while True:
        if key_list[ind] == jid:
            for i in range(0,len(data_list[ind])):
                answer.append(data_list[ind][i][1])
            key_list.pop(ind)
            data_list.pop(ind)
            return answer, key_list, data_list
            break
        elif key_list[ind] < jid:
            ind += 1
            continue
        else:
            return [], key_list, data_list
            break

def call_Data(key_list, data_list, jid):
    A = []
    C = []
    G = []
    T = []
    ind = 0
    while True:
        if key_list[ind] == jid and len(data_list[ind][0]) == 4:
            for position in range(0,len(data_list[ind])):
                A.append(data_list[ind][position][0][3])
                C.append(data_list[ind][position][1][3])
                G.append(data_list[ind][position][2][3])
                T.append(data_list[ind][position][3][3])
            key_list.pop(ind)
            data_list.pop(ind)
            return A, C, G, T, key_list, data_list
            break
        elif len(data_list[ind][0]) != 4:
            key_list.pop(ind)
            data_list.pop(ind)
            continue
        elif key_list[ind] != jid:
            ind += 1
            continue
        else:
            return [],[],[],[], key_list, data_list
            break


""" - - - -  Collecting and Parsing Data From Website then Adding to Session DB - -  - - """

def parse_Jaspar_db(session, database_url):

    f = urlopen(database_url+'MATRIX_PROTEIN.txt')
    uniprots = reader(f, delimiter = '\t')
    uniprots = make_jaspar_int(uniprots)
    sortedlist = sort(uniprots)
    uniprots, key_uniprot = group_by_jaspar_id(sortedlist)
    print 'Collected Uniprot Data'

    f = urlopen(database_url+'MATRIX_SPECIES.txt')
    NCBI = reader(f, delimiter = '\t')
    NCBI = make_jaspar_int(NCBI)
    sortedlist = sort(NCBI)
    NCBI, key_NCBI = group_by_jaspar_id(sortedlist)
    print 'Collected Species Data'

    f = urlopen(database_url+'MATRIX_ANNOTATION.txt')
    annots = reader(f, delimiter = '\t')
    annots = make_jaspar_int(annots)
    sortedlist = sort(annots)
    class_info = []
    fam_info = []
    med_info = []
    type_info = []
    for i in range(0,len(sortedlist)):
        if sortedlist[i][1] == 'class':
            class_info.append([sortedlist[i][0], sortedlist[i][2]])
        elif sortedlist[i][1] == 'family':
            fam_info.append([sortedlist[i][0], sortedlist[i][2]])
        elif sortedlist[i][1] == 'medline':
            med_info.append([sortedlist[i][0], sortedlist[i][2]])
        elif sortedlist[i][1] == 'type':
            type_info.append([sortedlist[i][0], sortedlist[i][2]])
    class_info, key_class = group_by_jaspar_id(class_info)
    fam_info, key_fam = group_by_jaspar_id(fam_info)
    med_info, key_med = group_by_jaspar_id(med_info)
    type_info, key_type = group_by_jaspar_id(type_info)
    print 'Collected Annotation Data'

    f = urlopen(database_url+'MATRIX_DATA.txt')
    matrix_data = reader(f, delimiter = '\t')
    matrix_data = make_data_int(matrix_data)
    sortedlist = sort(matrix_data)
    matrix_data, key_data = group_by_jaspar_id(sortedlist)
    for i in range(0,len(matrix_data)):
        matrix_data[i] = sorted(matrix_data[i], key=itemgetter(2), reverse=False)
        matrix_data[i] = group_by_position(matrix_data[i])
    print 'Collected Matrix Binding Data'


""" - - - - - - - - - - - - - - - - Exporting to DB - - - - - - - - - - - - - - - """

    f = urlopen(database_url+'MATRIX.txt')
    lines = f.readlines()
    print 'Table Creation Started'
    for i in range(0, len(lines)):
        unit = lines[i].split('\t')

        jaspar_matrix_ID = int(unit[0])
        assert isinstance(jaspar_matrix_ID, int)

        jaspar_tf_ID = unit[2]
        assert isinstance(jaspar_tf_ID, str)

        version = int(unit[3])
        assert isinstance(version, int)

        tf_name = unit[4]
        assert isinstance(tf_name, str)

        jaspar_collection = unit[1]
        assert isinstance(jaspar_collection, str)

        uniprot_id, key_uniprot, uniprots = call_Attribute(key_uniprot, uniprots, jaspar_matrix_ID)
        NCBI_id, key_NCBI, NCBI = call_Attribute(key_NCBI, NCBI, jaspar_matrix_ID)

        for n in range(0,len(uniprot_id)):
            assert isinstance(uniprot_id[n], str)
            transcriptionfactor = TranscriptionFactor(uniprot_id = uniprot_id[n], jaspar_matrix_ID = jaspar_matrix_ID)
            session.add(transcriptionfactor)

        for n in range(0,len(NCBI_id)):
            assert isinstance(NCBI_id[n], str)
            species = Species(NCBI_id = NCBI_id[n], jaspar_matrix_ID = jaspar_matrix_ID)
            session.add(species)

        class_name, key_class, class_info = call_Attribute(key_class, class_info, jaspar_matrix_ID)
        for n in range(0,len(class_name)):
            assert isinstance(class_name[n], str)
            proteinclass = Class(class_name = class_name[n], jaspar_matrix_ID = jaspar_matrix_ID)
            session.add(proteinclass)

        fam_name, key_fam, fam_info = call_Attribute(key_fam, fam_info, jaspar_matrix_ID)
        for n in range(0,len(fam_name)):
            u = unicode(fam_name[n], 'utf-8')
            assert isinstance(u, unicode)
            family = Family(family_name = u, jaspar_matrix_ID = jaspar_matrix_ID)
            session.add(family)

        medline_id, key_med, med_info = call_Attribute(key_med, med_info, jaspar_matrix_ID)
        for n in range(0,len(medline_id)):
            assert isinstance(medline_id[n], str)
            resources = Resources(medline_id = medline_id[n], jaspar_matrix_ID = jaspar_matrix_ID)
            session.add(resources)

        type_name, key_type, type_info = call_Attribute(key_type, type_info, jaspar_matrix_ID)
        for n in range(0,len(type_name)):
            assert isinstance(type_name[n], str)
            type_ = Type(type_name = type_name[n], jaspar_matrix_ID = jaspar_matrix_ID)
            session.add(type_)

        frequency_A, frequency_C, frequency_G, frequency_T, key_data, matrix_data = call_Data(key_data, matrix_data, jaspar_matrix_ID)
        for position in range(0,len(frequency_A)):
            A = frequency_A[position]
            assert isinstance(A, float)
            C = frequency_C[position]
            assert isinstance(C, float)
            G = frequency_G[position]
            assert isinstance(G, float)
            T = frequency_T[position]
            assert isinstance(T, float)
            position = position+1
            assert isinstance(position, int)
            bindingmatrix = BindingMatrix(position = position, frequency_A = A,
                                        frequency_C = C, frequency_G = G,
                                        frequency_T = T, jaspar_matrix_ID = jaspar_matrix_ID)
            session.add(bindingmatrix)


        matrixobservation = MatrixObservation(jaspar_matrix_ID = jaspar_matrix_ID, jaspar_tf_ID = jaspar_tf_ID,
                                                version = version, tf_name = tf_name, jaspar_collection = jaspar_collection,
                                                data = [bindingmatrix], references = [resources],
                                                type_info = [type_], tf = [transcriptionfactor],
                                                family_info = [family], class_info = [proteinclass], species = [species])

        session.add(matrixobservation)

    return
