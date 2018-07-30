import time

def timeit(method):

    def timed(*args, **kw):
        if args[0].verbose:
            print('\n------------------------ Initializing %r ------------------------' % (method.__name__))
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        if args[0].verbose:
            print('%r took %2.2f sec' % \
                      (method.__name__, te-ts))
            print('%r completed' % (method.__name__))
        return result

    return timed
