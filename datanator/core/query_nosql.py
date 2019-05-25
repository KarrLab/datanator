from datanator.util import mongo_util
import time

class DataQuery(mongo_util.MongoUtil):

    '''Collection agnostic queries
    '''

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet= None, db=None,
                verbose=False, max_entries=float('inf')):
        super(DataQuery, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet, 
                                    db=db, verbose=verbose, max_entries=max_entries)

    '''TODO:1. Make query language more user friendly
            2. Full text search for all key:value pairs
    '''

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
        s=None, e=None, inbatch=False, query=None, batch_callback=None, projection=None):
        '''An iterator for returning docs in a collection, with batch query.
           additional filter query can be passed via "query", e.g.,
           doc_feeder(collection_str, query={'taxid': {'$in': [9606, 10090, 10116]}})
           batch_callback is a callback function as fn(cnt, t), called after every batch
           fields is optional parameter passed to find to restrict fields to return.
        '''
        collection = self.fill_db(collection_str, sym_link=sym_link)
        cur = collection.find(query, no_cursor_timeout=False, projection=projection)
        n = cur.count()
        s = s or 0
        e = e or n
        if self.verbose:
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
                if self.verbose:
                    print("Skipping %d documents." % s)
            if e:
                cur.limit(e - (s or 0))
            cur.batch_size(step)
            if self.verbose:
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
                    if self.verbose:
                        print('Done.[%.1f%%,%s]' % (cnt * 100. / n, self.timesofar(t1)))
                    if batch_callback:
                        batch_callback(cnt, time.time()-t1)
                    if cnt < e:
                        t1 = time.time()
                        if self.verbose:
                            print("Processing %d-%d documents..." %
                                  (cnt + 1, min(cnt + step, e)), end='')
            if inbatch and doc_li:
                # Important: need to yield the last batch here
                yield doc_li

            # print 'Done.[%s]' % timesofar(t1)
            if self.verbose:
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


class QuerySabio(DataQuery):
    '''Queries specific to sabio_rk collection
    '''
    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet= None, db=None,
                collection_str='sabio_rk', verbose=False, max_entries=float('inf')):
        self.collection_str = collection_str
        super(DataQuery, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, 
                replicaSet= replicaSet, db=db,
                verbose=verbose, max_entries=max_entries)

    def find_reaction_participants(self, kinlaw_id):
        ''' Find the reaction participants defined in sabio_rk using kinetic law id
            Args:
                kinlaw_id: list of kinlaw_id to search for
            Return:
                rxns: list of dictionaries containing names of reaction participants
                [{'substrates': [], 'products': [] }, ... {} ]
        '''
        query = {'kinlaw_id': {'$in': kinlaw_id} }
        docs = self.doc_feeder(collection_str = self.collection_str, query=query)
        rxns = []
        i = 0
        for doc in docs:
            if i == self.max_entries:
                break 
            doc.pop('_id', None)
            doc_flat = self.flatten_json(doc)

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