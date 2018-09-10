import time

def timemethod(method):

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

def timeloadcontent(method):

    def timed(*args, **kw):
        if args[0].verbose:
            print(''' \n
                ===================================
                |                                 |
                |                                 |
                |    Starting Datanator Build     |
                |                                 |
                |                                 |
                ===================================

                ''')

        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        if args[0].verbose:
            print(''' \n
                =============================================
                |                                           |
                |             Finished Build                |
                    Total time taken for build: %2.2f secs
                |                                           |
                =============================================
                ''' % (te - ts))

        return result

    return timed


def continuousload(method):

    def continuous(*args, **kw):
        try:
            result = method(*args, **kw)
            return result
        except Exception as e:
            print(e)

    return continuous
