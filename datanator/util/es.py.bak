'''Because elasticsearch is inadequate

:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2019-04-02
:Copyright: 2019, Karr Lab
:License: MIT
'''

from __future__ import print_function
import sys
import time
import re

from elasticsearch import Elasticsearch
from elasticsearch import helpers

from config import ES_HOST, ES_INDEX_NAME
from util.mongo import doc_feeder, ask, timesofar

import logging
formatter = logging.Formatter("%(levelname)s:%(name)s:%(message)s")
es_logger = logging.getLogger('elasticsearch')
es_logger.setLevel(logging.WARNING)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
es_logger.addHandler(ch)

es_tracer = logging.getLogger('elasticsearch.trace')
es_tracer.setLevel(logging.WARNING)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
es_tracer.addHandler(ch)


def get_es(es_host=None):
    es_host = es_host or ES_HOST
    es = Elasticsearch(es_host, timeout=600, max_retries=10)
    return es


def lastexception():
    exc_type, exc_value, tb = sys.exc_info()
    if exc_type is None:
        print("No exception occurs.")
        return
    print(exc_type.__name__ + ':', end='')
    try:
        excArgs = exc_value.__dict__["args"]
    except KeyError:
        excArgs = ()
    return str(exc_type)+':'+''.join([str(x) for x in excArgs])


class ESIndexer(object):
    def __init__(self, es_index_name=None, mapping=None, es_host=None, step=5000):
        self.conn = get_es(es_host)
        self.ES_INDEX_NAME = es_index_name or ES_INDEX_NAME

        self.step = step
        #self.conn.bulk_size = self.step
        self.number_of_shards = 5      # set number_of_shards when create_index
        self.s = None   # optionally, can specify number of records to skip,
                        # useful to continue indexing after an error.
        self.use_parallel = False
        self._mapping = mapping

    def _get_es_version(self):
        info = self.conn.info()
        return info['version']['number']

    def check(self):
        '''print out ES server info for verification.'''
        # print "Servers:", self.conn.servers
        print("Servers:", '|'.join(["{host}:{port}".format(**h) for h in self.conn.transport.hosts]))
        print("Default indices:", self.ES_INDEX_NAME)

    def create_index(self, index_name=None, mapping=None):
        index_name = index_name or self.ES_INDEX_NAME
        if not self.conn.indices.exists(index_name):
            body = {
                'settings': {
                    'number_of_shards': self.number_of_shards,
                    "number_of_replicas": 0,    # set this to 0 to boost indexing
                                                # after indexing, set "auto_expand_replicas": "0-all",
                                                #   to make additional replicas.
                }
            }
            if mapping:
                mapping = {"mappings": mapping}
                body.update(mapping)
            print(self.conn.indices.create(index=index_name, body=body))

    def exists_index(self, index):
        return self.conn.indices.exists(index)

    def verify_mapping(self, update_mapping=False):
        '''Verify if index and mapping exist, update mapping if mapping does not exist,
           or "update_mapping=True" explicitly
        '''
        conn = self.conn
        index_name = self.ES_INDEX_NAME

        # Test if index exists
        if not self.exists_index(index_name):
            print('Error: index "%s" does not exist. Create it first.' % index_name)
            return -1

        mapping = conn.indices.get_mapping(index=index_name)
        empty_mapping = mapping == {}

        if empty_mapping:
            #if no existing mapping available for index_type
            #force update_mapping to True
            update_mapping = True

        if update_mapping:
            print("Updating mapping...", end='')
            if not empty_mapping:
                print("\n\tRemoving existing mapping...", end='')
                print(conn.indices.delete_mapping(index=index_name))
            self.get_field_mapping()
            print(conn.indices.put_mapping(index=index_name,
                                           body=self._mapping))

    def update_mapping_meta(self, meta):
        index_name = self.ES_INDEX_NAME

        if isinstance(meta, dict) and len(set(meta) - set(['_meta', '_timestamp'])) == 0:
            body = {index_name: meta}
            print(self.conn.indices.put_mapping(
                index=index_name,
                body=body)
            )
        else:
            raise ValueError('Input "meta" should have and only have "_meta" field.')

    def count(self, query=None):
        conn = self.conn
        index_name = self.ES_INDEX_NAME
        if query:
            if isinstance(query, dict) and 'query' in query:
                return conn.count(index_name, query)
            else:
                raise ValueError("Not a valid input query")
        else:
            return conn.count(index_name)

    def get(self, id, **kwargs):
        '''get a specific doc by its id.'''
        conn = self.conn
        index_name = self.ES_INDEX_NAME
        return conn.get(index_name, id, **kwargs)

    def get_docs(self, ids, step=None, **kwargs):
        '''return matching docs for given ids, if not found return None.
           A generator is returned and the order is perserved.
        '''
        conn = self.conn
        index_name = self.ES_INDEX_NAME
        step = step or self.step
        for i in range(0, len(ids), step):
            _ids = ids[i:i + step]
            body = {'ids': _ids}
            res = conn.mget(body=body, index=index_name, **kwargs)
            for doc in res['docs']:
                if doc['found']:
                    yield doc['_source']
                else:
                    yield None

    def index(self, doc, id=None):
        '''add a doc to the index. If id is not None, the existing doc will be
           updated.
        '''
        return self.conn.index(self.ES_INDEX_NAME, doc, id=id)

    def index_bulk(self, docs, step=None):
        index_name = self.ES_INDEX_NAME
        step = step or self.step

        def _get_bulk(doc):
            doc.update({
                "_index": index_name,
            })
            return doc
        actions = (_get_bulk(doc) for doc in docs)
        return helpers.bulk(self.conn, actions, chunk_size=step)

    def add_docs(self, docs, step=None):
        self.index_bulk(docs, step=step)
        self.conn.indices.flush()
        self.conn.indices.refresh()

    def delete_doc(self, id):
        '''delete a doc from the index based on passed id.'''
        return self.conn.delete(self.ES_INDEX_NAME, id)

    def delete_docs(self, ids, step=None):
        index_name = self.ES_INDEX_NAME
        step = step or self.step

        def _get_bulk(_id):
            doc = {
                '_op_type': 'delete',
                "_index": index_name,
                "_id": _id
            }
            return doc
        actions = (_get_bulk(_id) for _id in ids)
        return helpers.bulk(self.conn, actions, chunk_size=step, stats_only=True, raise_on_error=False)

    def update(self, id, extra_doc, bulk=False):
        '''update an existing doc with extra_doc.'''
        conn = self.conn
        index_name = self.ES_INDEX_NAME

        if not bulk:
            body = {'doc': extra_doc}
            return conn.update(index_name, id, body)
        else:
            raise NotImplementedError

    def update_docs(self, partial_docs, **kwargs):
        index_name = self.ES_INDEX_NAME
        def _get_bulk(doc):
            doc = {
                '_op_type': 'update',
                "_index": index_name,
                "_id": doc['_id'],
                "doc": doc
            }
            return doc
        actions = (_get_bulk(doc) for doc in partial_docs)
        return helpers.bulk(self.conn, actions, chunk_size=self.step, **kwargs)

    def wait_till_all_shards_ready(self, timeout=None, interval=5):
        raise NotImplementedError

    def optimize(self):
        '''optimize the default index.'''
        return self.conn.indices.optimize(self.ES_INDEX_NAME,
                                          wait_for_merge=False,   # True,
                                          max_num_segments=5)

    def optimize_all(self):
        """optimize all indices"""
        return self.conn.indices.optimize('', wait_for_merge=False,  # True,
                                          max_num_segments=5)

    def get_field_mapping(self):
        # raise NotImplementedError
        return self._mapping

    def build_index(self, collection, update_mapping=False, verbose=False, query=None, bulk=True):
        conn = self.conn
        index_name = self.ES_INDEX_NAME

        self.verify_mapping(update_mapping=update_mapping)
        #update some settings for bulk indexing
        body = {
            "index": {
                #"refresh_interval": "-1",              # disable refresh temporarily
                "auto_expand_replicas": "0-all",
                #"number_of_replicas": 0,
                "refresh_interval": "30s",
            }
        }
        conn.indices.put_settings(body, index_name)
        try:
            print("Building index...")
            if self.use_parallel:
                cnt = self._build_index_parallel(collection, verbose)
            else:
                cnt = self._build_index_sequential(collection, verbose, query=query, bulk=bulk)
        finally:
            #restore some settings after bulk indexing is done.
            body = {
                "index": {
                    "refresh_interval": "1s"              # default settings
                }
            }
            conn.indices.put_settings(body, index_name)

            try:
                print("Flushing...", conn.indices.flush())
                print("Refreshing...", conn.indices.refresh())
            except:
                pass

            time.sleep(1)
            print("Validating...", end='')
            target_cnt = collection.find(query).count()
            es_cnt = self.count()['count']
            if target_cnt == es_cnt:
                print("OK [total count={}]".format(target_cnt))
            else:
                print("\nWarning: total count of gene documents does not match [{}, should be {}]".format(es_cnt, target_cnt))

        if cnt:
            print('Done! - {} docs indexed.'.format(cnt))
            print("Optimizing...", self.optimize())

    def _build_index_sequential(self, collection, verbose=False, query=None, bulk=True):

        def rate_control(cnt, t):
            delay = 0
            if t > 90:
                delay = 30
            elif t > 60:
                delay = 10
            if delay:
                print("\tPausing for {}s...".format(delay), end='')
                time.sleep(delay)
                print("done.")

        src_docs = doc_feeder(collection, step=self.step, s=self.s, batch_callback=rate_control, query=query)
        if bulk:
            res = self.index_bulk(src_docs)
            if len(res[1]) > 0:
                print("Error: {} docs failed indexing.".format(len(res[1])))
            return res[0]
        else:
            cnt = 0
            for doc in src_docs:
                self.index(doc)
                cnt += 1
                if verbose:
                    print(cnt, ':', doc['_id'])
            return cnt

    def _build_index_parallel(self, collection, verbose=False):
        raise NotImplementedError
        from utils.parallel import (run_jobs_on_ipythoncluster,
                                    collection_partition,
                                    require)
        kwargs_common = {'ES_HOST': ES_HOST,
                         'ES_INDEX_NAME': self.ES_INDEX_NAME
                         }
        task_list = []
        for kwargs in collection_partition(collection, step=self.step):
            kwargs.update(kwargs_common)
            task_list.append(kwargs)

        @require('mongokit', 'pyes')
        def worker(kwargs):
            import mongokit
            import pyes
            server = kwargs['server']
            port = kwargs['port']
            src_db = kwargs['src_db']
            src_collection = kwargs['src_collection']
            skip = kwargs['skip']
            limit = kwargs['limit']

            mongo_conn = mongokit.Connection(server, port)
            src = mongo_conn[src_db]

            ES_HOST = kwargs['ES_HOST']
            ES_INDEX_NAME = kwargs['ES_INDEX_NAME']

            es_conn = pyes.ES(ES_HOST, default_indices=[ES_INDEX_NAME],
                              timeout=120.0, max_retries=10)

            cur = src[src_collection].find(skip=skip, limit=limit, timeout=False)
            cur.batch_size(1000)
            cnt = 0
            try:
                for doc in cur:
                    es_conn.index(doc, ES_INDEX_NAME, doc['_id'], bulk=True)
                    cnt += 1
            finally:
                cur.close()
            es_conn.indices.flush()   # this is important to avoid missing docs
            es_conn.indices.refresh()
            return cnt

        job_results = run_jobs_on_ipythoncluster(worker, task_list)
        if job_results:
            cnt = sum(job_results)
            return cnt

    def doc_feeder(self, index_name=None, step=10000, verbose=True, query=None, scroll='10m', **kwargs):
        conn = self.conn
        index_name = index_name or self.ES_INDEX_NAME

        n = self.count(query=query)['count']
        cnt = 0
        t0 = time.time()
        if verbose:
            print('\ttotal docs: {}'.format(n))

        _kwargs = kwargs.copy()
        _kwargs.update(dict(size=step, index=index_name))
        res = helpers.scan(conn, query=query, scroll=scroll, **_kwargs)
        t1 = time.time()
        for doc in res:
            if verbose and cnt % step == 0:
                if cnt != 0:
                    print('done.[%.1f%%,%s]' % (cnt*100./n, timesofar(t1)))
                print('\t{}-{}...'.format(cnt+1, min(cnt+step, n)), end='')
                t1 = time.time()
            yield doc
            cnt += 1
        if verbose:
            print('done.[%.1f%%,%s]' % (cnt*100./n, timesofar(t1)))
            print("Finished! [{}]".format(timesofar(t0)))

    def get_id_list(self, index_name=None, step=100000, verbose=True):
        cur = self.doc_feeder(index_name=index_name, step=step, fields=[], verbose=verbose)
        id_li = [doc['_id'] for doc in cur]
        return id_li

    def get_id_list_parallel(self, taxid_li, index_name=None, step=1000, verbose=True):
        '''return a list of all doc ids in an index_type.'''
        raise NotImplementedError
        from utils.parallel import run_jobs_on_ipythoncluster

        def _get_ids_worker(args):
            from utils.es import ESIndexer
            from pyes import MatchAllQuery
            es_kwargs, start, step = args
            q = MatchAllQuery().search()
            q.sort = [{'entrezgene': 'asc'}, {'ensembl.gene': 'asc'}]
            q.fields = []
            q.start = start
            q.size = step
            esi = ESIndexer(**es_kwargs)
            cnt = esi.count()['count']
            res = esi.conn.search_raw(q)
            assert res['hits']['total'] == cnt
            return [doc['_id'] for doc in res['hits']['hits']]

        def _get_ids_worker_by_taxid(args):
            from utils.es import ESIndexer
            from pyes import TermQuery
            es_kwargs, taxid, step = args
            q = TermQuery()
            q.add('taxid', taxid)
            q.fields = []
            q.size = step
            esi = ESIndexer(**es_kwargs)
            res = esi.conn.search(q)
            xli = [doc['_id'] for doc in res]
            assert len(xli) == res.total
            return xli

        es_kwargs = {'es_index_name': self.ES_INDEX_NAME, 'es_host': 'su02:9200'}
        task_li = [(es_kwargs, taxid, step) for taxid in taxid_li]
        #print task_li
        job_results = run_jobs_on_ipythoncluster(_get_ids_worker_by_taxid, task_li)
        return job_results

    def clone_index(self, src_index, target_index, target_es_host=None, step=10000, scroll='10m',
                    target_index_settings=None, number_of_shards=None):
        '''clone src_index to target_index on the same es_host, or another one given
           by target_es_host.

           This method can now be replaced by helpers.reindex
        '''


def es_clean_indices(keep_last=2, es_host=None, verbose=True, noconfirm=False, dryrun=False):
    '''clean up es indices, only keep last <keep_last> number of indices.'''
    conn = get_es(es_host)
    index_li = list(conn.indices.get_aliases().keys())
    if verbose:
        print("Found {} indices".format(len(index_li)))

    for prefix in ('genedoc_mygene', 'genedoc_mygene_allspecies'):
        pat = prefix + '_(\d{8})_\w{8}'
        _li = []
        for index in index_li:
            mat = re.match(pat, index)
            if mat:
                _li.append((mat.group(1), index))
        _li.sort()   # older collection appears first
        index_to_remove = [x[1] for x in _li[:-keep_last]]   # keep last # of newer indices
        if len(index_to_remove) > 0:
            print ("{} \"{}*\" indices will be removed.".format(len(index_to_remove), prefix))
            if verbose:
                for index in index_to_remove:
                    print ('\t', index)
            if noconfirm or ask("Continue?") == 'Y':
                for index in index_to_remove:
                    if dryrun:
                        print("dryrun=True, nothing is actually deleted")
                    else:
                        conn.indices.delete(index)
                print("Done.[%s indices removed]" % len(index_to_remove))
            else:
                print("Aborted.")
        else:
            print("Nothing needs to be removed.")


def get_latest_indices(es_host=None):
    conn = get_es(es_host)
    index_li = list(conn.indices.get_aliases().keys())
    print("Found {} indices".format(len(index_li)))

    latest_indices = []
    for prefix in ('genedoc_mygene', 'genedoc_mygene_allspecies'):
        pat = prefix + '_(\d{8})_\w{8}'
        _li = []
        for index in index_li:
            mat = re.match(pat, index)
            if mat:
                _li.append((mat.group(1), index))
        if not _li:
            print("Warning: no matching indices found!")
            continue
        latest_indices.append(sorted(_li)[-1])

    if len(latest_indices) == 2:
        if latest_indices[0][0] != latest_indices[1][0]:
            print ("Warning: unmatched timestamp:")
            print ('\n'.join([x[1] for x in latest_indices]))
        latest_indices = [x[1] for x in latest_indices]
        return latest_indices
