from datanator.util import mongo_util
from datanator.util import chem_util
from datanator.util import file_util
import time
import hashlib
import numpy as np

class DataQuery(mongo_util.MongoUtil):

    '''Collection agnostic queries
    '''

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet= None, db=None,
                verbose=False, max_entries=float('inf'), username = None, 
                 password = None, authSource = 'admin'):
        super(DataQuery, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet, 
                                    db=db, verbose=verbose, max_entries=max_entries, username = username, 
                                    password = password, authSource = authSource)


    def find_text(self, v, collection=None):
        ''' Find documents containing string v
            v cannot be part of a word (e.g. 'wor' of 'word')
            v needs to be in a previously indexed field
            Args:
                v: value to be matched
        '''
        if collection == None:
            return None
        else:
            _, _, col_obj = self.con_db(collection)
            return col_obj.find({'$text': {'$search': v}},
                                { 'score': { '$meta': "textScore" } }).sort( { 'score': { '$meta': "textScore" } } )

    def doc_feeder(self,collection_str=None, sym_link = False, step=1000, 
        s=None, e=None, inbatch=False, query=None, 
        batch_callback=None, projection=None, verbose=False):
        '''An iterator for returning docs in a collection, with batch query.
           additional filter query can be passed via "query", e.g.,
           doc_feeder(collection_str, query={'taxid': {'$in': [9606, 10090, 10116]}})
           batch_callback is a callback function as fn(cnt, t), called after every batch
           fields is optional parameter passed to find to restrict fields to return.
        '''
        _, _, collection = self.con_db(collection_str)
        cur = collection.find(query, no_cursor_timeout=True, projection=projection)
        n = cur.count()
        s = s or 0
        e = e or n
        if verbose:
            print('Retrieving %d documents from collection "%s".' %
                  (n, collection_str))
        t0 = time.time()
        if inbatch:
            doc_li = []
        cnt = 0
        t1 = time.time()
        try:
            if s:
                cur.skip(s)
                cnt = s
                if verbose:
                    print("Skipping %d documents." % s)
            if e:
                cur.limit(e - (s or 0))
            cur.batch_size(step)
            if verbose:
                print("Processing %d-%d documents..." %
                  (cnt + 1, min(cnt + step, e)), end='')
            for doc in cur:
                if inbatch:
                    doc_li.append(doc)
                else:
                    yield doc
                cnt += 1
                if cnt % step == 0:
                    if inbatch:
                        yield doc_li
                        doc_li = []
                    if verbose:
                        print('Done.[%.1f%%,%s]' % (cnt * 100. / n, self.timesofar(t1)))
                    if batch_callback:
                        batch_callback(cnt, time.time()-t1)
                    if cnt < e:
                        t1 = time.time()
                        if verbose:
                            print("Processing %d-%d documents..." %
                                  (cnt + 1, min(cnt + step, e)), end='')
            if inbatch and doc_li:
                # Important: need to yield the last batch here
                yield doc_li

            # print 'Done.[%s]' % timesofar(t1)
            if verbose:
                print('Done.[%.1f%%,%s]' % (cnt * 100. / (n+1), self.timesofar(t1)))
                print("=" * 20)
                print('Finished.[total time: %s]' % self.timesofar(t0))
        finally:
            cur.close()

    def ask(self, prompt, options='YN'):
        '''Prompt Yes or No,return the upper case 'Y' or 'N'.
        '''
        options = options.upper()
        while 1:
            s = input(prompt+'[%s]' % '|'.join(list(options))).strip().upper()
            if s in options:
                break
        return s

    def timesofar(self, t0, clock=0, t1=None):
        '''return the string(eg.'4m5.32s') for the passed real time/CPU time so far
           from given t0 (clock 0 for cpu time 1 for realtime).
        '''
        t1 = t1 or time.clock() if clock else time.time()
        t = t1 - t0
        h = int(t / 3600)
        m = int((t % 3600) / 60)
        s = round((t % 3600) % 60, 2)
        t_str = ''
        if h != 0:
            t_str += '%sh' % h
        if m != 0:
            t_str += '%sm' % m
        t_str += '%ss' % s
        return t_str


class QueryMetabolitesMeta(DataQuery):
    '''Queries specific to metabolites_meta collection
    '''
    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet= None, db=None,
                collection_str='metabolites_meta', verbose=False, max_entries=float('inf'), username = None, 
                 password = None, authSource = 'admin'):
        self.collection_str = collection_str
        self.verbose = verbose
        super(DataQuery, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, 
                replicaSet= replicaSet, db=db,
                verbose=verbose, max_entries=max_entries, username = username, 
                 password = password, authSource = authSource)
        self.client, self.db_obj, self.collection = self.con_db(self.collection_str)
        self.file_manager = file_util.FileUtil()
        self.chem_manager = chem_util.ChemUtil()

    def get_metabolite_synonyms(self, compounds):
        ''' Find synonyms of a compound
            Args:
                compound: name(s) of the compound e.g. "ATP", ["ATP", "Oxygen", ...]
            Returns:
                synonyms: dictionary of synonyms of the compounds
                        {'ATP': [], 'Oxygen': [], ...}
                rxns: dictionary of rxns in which each compound is found
                    {'ATP': [12345,45678,...], 'Oxygen': [...], ...}
        '''
        synonyms = {}
        rxns = {}

        def find_synonyms_of_str(c):
            if len(c) != 0:
                query = {'synonyms.synonym': c}
                projection = {'synonyms.synonym': 1, '_id': -1, 'kinlaw_id': 1}
                collation = {'locale': 'en', 'strength': 2}
                doc = self.collection.find_one(filter = query, projection = projection, collation = collation)
                synonym = {}
                rxn = {}
                try:
                    synonym[c] = doc['synonyms']['synonym']
                    rxn[c] = doc['kinlaw_id']
                except TypeError as e:
                    synonym[c] = (c + ' does not exist in '+ self.collection_str)
                    rxn[c] = (c + ' does not exist in '+ self.collection_str)
                return rxn, synonym
            else:
                return ({'reactions': None}, {'synonyms': None})

        if len(compounds) == 0:
            return ({'reactions': None}, {'synonyms': None})
        elif isinstance(compounds, str):
            rxn, syn = find_synonyms_of_str(compounds)
            synonyms.update(syn)
            rxns.update(rxn)            
        else:
            for c in compounds:
                rxn, syn = find_synonyms_of_str(c)
                synonyms.update(syn)
                rxns.update(rxn)
        return rxns, synonyms

    def get_metabolite_inchi(self, compounds):
        '''Given a list of compound name(s)
            Return the corrensponding inchi string
            Args:
                compounds: list of compounds
                ['ATP', '2-Ketobutanoate']
            Return:
                ['....', 'InChI=1S/C4H6O3/c1-2-3(5)4(6)7']
        '''
        inchi = []
        projection = {'_id': 0, 'inchi': 1}
        collation = {'locale': 'en', 'strength': 2}
        for compound in compounds:
            cursor = self.collection.find_one({'synonyms.synonym': compound},
                                projection = projection, collation = collation)
            inchi.append(cursor['inchi'])
        return inchi

    def get_metabolite_hashed_inchi(self, compounds):
        ''' Given a list of compound name(s)
            Return the corresponding hashed inchi string
            Args:
                compounds: ['ATP', '2-Ketobutanoate']
            Return:
                hashed_inchi: ['3e23df....', '7666ffa....']
        '''
        hashed_inchi = []
        projection = {'_id': 0, 'inchi_hashed': 1}
        collation = {'locale': 'en', 'strength': 2}
        for compound in compounds:
            cursor = self.collection.find_one({'synonyms.synonym': compound},
                                projection = projection, collation = collation)
            hashed_inchi.append(cursor['inchi_hashed'])
        return hashed_inchi        

   
    def get_metabolite_name_by_hash(self, compounds):
        ''' Given a list of hashed inchi, 
            return a list of name (one of the synonyms)
            for each compound
            Args:
                compounds: list of compounds in hashed_inchi format
            Return:
                result: list of names
                    [name, name, name]
        '''
        result = []
        projection = {'_id': 0, 'synonyms.synonym': 1}
        for compound in compounds:
            cursor = self.collection.find_one({'inchi_hashed': compound},
                                            projection = projection)
            try:
                result.append(cursor['synonyms'])
            except KeyError:
                result.append('No synonyms')
        return [x['synonym'][-1] for x in result]


    def get_metabolite_similar_compounds(self, compounds, num = 0, threshold = 0):
        ''' Given a list of compound names
            Return the top num number of similar compounds
            with tanimoto score above threshold values
            Args:
                compounds: list of compound names
                num: number of similar compounds to return
                threshold: threshold tanimoto coefficient value
                return_format: return dictionary key format, either
                                hashed inchi or name
            Return:
                result: list of similar compounds and their tanimoto score
                [ {'compound1': score, 'compound2': score, ... 'compound_num': score},
                  {'compound1': score, 'compound2': score, ... 'compound_num': score}, ...]
                    compound(1-n) will be in name format
                raw: list of similar compounds and their tanimoto score
                [ {'compound1': score, 'compound2': score, ... 'compound_num': score},
                  {'compound1': score, 'compound2': score, ... 'compound_num': score}, ...]
                    compound(1-n) will be in hashed_inchi format
        '''
        result = []
        raw = []
        hashed_inchi = self.get_metabolite_hashed_inchi(compounds)
        projection = {'_id': 0, 'similar_compounds': 1}

        for item in hashed_inchi:
            cursor = self.collection.find_one({'inchi_hashed': item},
                                                projection = projection)
            compounds = cursor['similar_compounds']

            scores = list(compounds.values())
            hashes = list(compounds.keys())
            names = self.get_metabolite_name_by_hash(hashes[:num])
            # convert to numpy object for faster calculations
            scores = np.asarray(scores)
            indices = np.nonzero(scores >= threshold)
            size = indices[0].size 
            if size == 0:
                raw.append(['No similar compound above threshold'])
                result.append(['No similar compound above threshold'])
            elif size < num:
                raw.append(compounds)
                replaced = self.file_manager.replace_dict_key(compounds, names)
                result.append(replaced)
            else:
                first_num = self.file_manager.access_dict_by_index(compounds, num)
                raw.append(first_num)
                replaced = self.file_manager.replace_dict_key(first_num, names[:num])
                result.append(replaced)

        return raw, result


class QuerySabio(DataQuery):
    '''Queries specific to sabio_rk collection
    '''
    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet= None, db=None,
                collection_str='sabio_rk', verbose=False, max_entries=float('inf'), username = None, 
                 password = None, authSource = 'admin'):
        self.collection_str = collection_str
        super(DataQuery, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, 
                replicaSet= replicaSet, db=db,
                verbose=verbose, max_entries=max_entries, username = username, 
                 password = password, authSource = authSource)
        self.chem_manager = chem_util.ChemUtil()
        self.file_manager = file_util.FileUtil()
        self.client, self.db_obj, self.collection = self.con_db(self.collection_str)

    def find_reaction_participants(self, kinlaw_id):
        ''' Find the reaction participants defined in sabio_rk using kinetic law id
            Args:
                kinlaw_id: list of kinlaw_id to search for
            Return:
                rxns: list of dictionaries containing names of reaction participants
                [{'substrates': [], 'products': [] }, ... {} ]
        '''
        if isinstance(kinlaw_id, list):
            query = {'kinlaw_id': {'$in': kinlaw_id} }
        else:
            query = {'kinlaw_id': kinlaw_id}
        docs = self.doc_feeder(collection_str = self.collection_str, query=query)
        rxns = []
        i = 0
        for doc in docs:
            if i == self.max_entries:
                break 
            if i % 10 == 0:
                print ('Finding reaction participants for kinlaw_id {} ...'.format(doc['kinlaw_id']))
            doc.pop('_id', None)
            doc_flat = self.file_manager.flatten_json(doc)

            substrates = []
            products = []
            substrate_identifier = ['reaction_participant', 'substrate', 'name']
            product_identifier = ['reaction_participant', 'product', 'name']

            for k, v in doc_flat.items():
                if all(x in k for x in substrate_identifier) and 'compartment' not in k:
                    substrates.append(v)
                elif all(x in k for x in product_identifier) and 'compartment' not in k:
                    products.append(v)
                else:
                    continue
            rxn = {'substrates': substrates, 'products': products}

            rxns.append(rxn)
            i += 1 


        return rxns

    def get_kinlawid_by_inchi_slow(self, inchi):
        ''' Find the kinlaw_id defined in sabio_rk using 
            rxn participants' inchi string
            Args:
                inchi: list of inchi, all in one rxn
            Return:
                rxns: list of kinlaw_ids that satisfy the condition
                [id0, id1, id2,...,  ]
        '''
        short_inchi = [self.chem_manager.simplify_inchi(s) for s in inchi]
        inchi_exp = ['\"' + s + '\"' for s in short_inchi]
        inchi_str = ''
        for s in inchi_exp:
            inchi_str = inchi_str + s + ' '
        condition = { '$text': {'$search': inchi_str} }
        projection = {'kinlaw_id': 1, '_id': 0}
        col = self.db_obj[self.collection_str]
        if self.verbose:
            print("\nQuerying text {} in collection {} ...".format(inchi_str, self.collection_str))
        cursor = col.find(filter = condition, projection = projection)
        _id = []
        for doc in cursor:
            _id.append(doc['kinlaw_id'])

        return _id

    def get_kinlawid_by_inchi(self, inchi):
        ''' Find the kinlaw_id defined in sabio_rk using 
            rxn participants' inchi string
            Args:
                inchi: list of inchi, all in one rxn
            Return:
                rxns: list of kinlaw_ids that satisfy the condition
                [id0, id1, id2,...,  ]
        '''
        short_inchi = [self.chem_manager.simplify_inchi(s) for s in inchi]
        hashed_inchi = [hashlib.sha224(s.encode()).hexdigest() for s in short_inchi]
        substrate = 'reaction_participant.substrate.hashed_inchi'
        product = 'reaction_participant.product.hashed_inchi'
        projection = {'kinlaw_id': 1}
        
        id_tally = []
        for inchi in hashed_inchi:
            ids = []
            query = {'$or': [ {substrate: inchi}, {product: inchi} ] }
            cursor = self.collection.find(filter = query, projection = projection)
            for doc in cursor:
                ids.append(doc['kinlaw_id'])
            id_tally.append(ids)

        return list(set(id_tally[0]).intersection(*id_tally))

    def get_kinlawid_by_rxn(self, substrates, products):
        ''' Find the kinlaw_id defined in sabio_rk using 
            rxn participants' inchi string
            Args:
                substrates: list of substrates' inchi
                products: list of products' inchi
            Return:
                rxns: list of kinlaw_ids that satisfy the condition
                [id0, id1, id2,...,  ]
        '''

        def get_kinlawid(inchi, side = 'substrate'):
            ''' Find the kinlaw_id defined in sabio_rk using 
                rxn participants' inchi string
                Args:
                    inchi: list of inchi, all in one rxn, on one side
                Return:
                    rxns: list of kinlaw_ids that satisfy the condition
                    [id0, id1, id2,...,  ]
            '''
            short_inchi = [self.chem_manager.simplify_inchi(s) for s in inchi]
            hashed_inchi = [hashlib.sha224(s.encode()).hexdigest() for s in short_inchi]

            substrate = 'reaction_participant.substrate.hashed_inchi'
            product = 'reaction_participant.product.hashed_inchi'
            projection = {'kinlaw_id': 1, '_id': 0}
            
            id_tally = []
            if side == 'substrate':
                for inchi in hashed_inchi:
                    ids = []
                    query = {substrate: inchi}
                    cursor = self.collection.find(filter = query, projection = projection)
                    for doc in cursor:
                        ids.append(doc['kinlaw_id'])
                    id_tally.append(ids)

                return list(set(id_tally[0]).intersection(*id_tally))
            else:

                for inchi in hashed_inchi:
                    ids = []
                    query = {product: inchi}
                    cursor = self.collection.find(filter = query, projection = projection)
                    for doc in cursor:
                        ids.append(doc['kinlaw_id'])
                    id_tally.append(ids)

                return list(set(id_tally[0]).intersection(*id_tally))

        sub_id = get_kinlawid(substrates, side = 'substrate')
        pro_id = get_kinlawid(products, side = 'product')
        result = list(set(sub_id) & set(pro_id))

        return result


class QueryTaxonTree(DataQuery):
    '''Queries specific to taxon_tree collection
    '''
    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet= None, db=None,
                collection_str='taxon_tree', verbose=False, max_entries=float('inf'), username = None, 
                 password = None, authSource = 'admin'):
        self.collection_str = collection_str
        super(DataQuery, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, 
                replicaSet= replicaSet, db=db,
                verbose=verbose, max_entries=max_entries, username = username, 
                 password = password, authSource = authSource)
        self.chem_manager = chem_util.ChemUtil()
        self.file_manager = file_util.FileUtil()
        self.client, self.db_obj, self.collection = self.con_db(self.collection_str)

    def get_anc_id_by_name(self, names):
        ''' Get organism's ancestor ids by
            using organism's names
            Args:
                names: list of organism's names e.g. Candidatus Diapherotrites
            Return:
                result: list of ancestors in order of the farthest to the closest
        '''
        result = []
        
        collation = {'locale': 'en', 'strength': 2}
        projection = {'_id': 0, 'anc_id':1}
        for name in names:
            query = {'tax_name': name}
            cursor = self.collection.find_one(query, collation = collation,
                                             projection = projection)
            result.append(cursor['anc_id'])
        return result

    def get_anc_id_by_id(self, ids):
        ''' Get organism's ancestor ids by
            using organism's ids
            Args:
                ids: list of organism's ids e.g. Candidatus Diapherotrites
            Return:
                result: list of ancestors in order of the farthest to the closest
        '''
        result = []
        
        collation = {'locale': 'en', 'strength': 2}
        projection = {'_id': 0, 'anc_id':1}
        for _id in ids:
            query = {'tax_id': _id}
            cursor = self.collection.find_one(query, collation = collation,
                                             projection = projection)
            result.append(cursor['anc_id'])
        return result

    def get_common_ancestor(self, org1, org2, org_format = 'name'):
        ''' Get the closest common ancestor between
            two organisms and their distances to the 
            said ancestor
            Args:
                org1: organism 1
                org2: organism 2
                org_format: the format of organism eg tax_id or tax_name
            Return:
                ancestor: closest common ancestor's id
                distance: each organism's distance to the ancestor
        '''
        if org_format == 'name':
            anc_ids = self.get_anc_id_by_name([org1, org2])
        else:
            anc_ids = self.get_anc_id_by_id([org1, org2])

        org1_anc = anc_ids[0]
        org2_anc = anc_ids[1]
        
        ancestor = self.file_manager.get_common(org1_anc, org2_anc)
        idx_org1 = org1_anc.index(ancestor)
        idx_org2 = org2_anc.index(ancestor)

        distance1 = len(org1_anc) - (idx_org1)
        distance2 = len(org2_anc) - (idx_org2)

        return (ancestor, [distance1, distance2])
