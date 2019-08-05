import hashlib
from .base26 import base26_triplet_1, base26_triplet_2, base26_triplet_3, base26_triplet_4, \
                         base26_dublet_for_bits_56_to_64, base26_dublet_for_bits_28_to_36
                         
class ChemUtil:

    def __init__(self):
        self.INCHI_STRING_PREFIX = "InChI="
        self.LEN_INCHI_STRING_PREFIX = len(self.INCHI_STRING_PREFIX)

    def simplify_inchi(self, inchi= 'InChI = None'):
        '''Remove molecules's protonation state
        "InChI=1S/H2O/h1H2" = > "InChI=1S/H2O"
        '''
        # if self.verbose:
        #     print('Parsing inchi by taking out protonation state')
        try:
            inchi_neutral = inchi.split('/h')[0]
            return inchi_neutral
        except AttributeError:
            return 'InChI = None'

    def hash_inchi(self, inchi = 'InChI = None'):
        ''' Hash inchi string using sha224
        '''
        try:
            hashed_inchi = hashlib.sha224(inchi.encode()).hexdigest()
            return hashed_inchi
        except AttributeError:
            return 'InChI = None'

    def get_sha256(self, text):
        hash = hashlib.sha256()
        hash.update(text.encode('ascii'))
        digest = hash.digest()
        digest_bytes_list = [ord(digest_byte) for digest_byte in digest] if isinstance(digest, str) else list(digest)
        return digest_bytes_list


    def inchi_to_inchikey(self, szINCHISource):
        '''
            fork from git@github.com:mnowotka/chembl_ikey.git
        '''
        flagstd = 'S'
        flagnonstd = 'N'
        flagver = 'A'
        flagproto = 'N'
        pplus = "OPQRSTUVWXYZ"
        pminus = "MLKJIHGFEDCB"

        if not szINCHISource:
            return None
        slen = len(szINCHISource)

        if slen < self.LEN_INCHI_STRING_PREFIX + 3:
            return None

        if not szINCHISource.startswith(self.INCHI_STRING_PREFIX):
            return None

        if szINCHISource[self.LEN_INCHI_STRING_PREFIX] != '1':
            return None

        pos_slash1 = self.LEN_INCHI_STRING_PREFIX + 1

        if szINCHISource[pos_slash1] == 'S':
            bStdFormat = 1
            pos_slash1 += 1

        if szINCHISource[pos_slash1] != '/':
            return None

        if not szINCHISource[pos_slash1+1].isalnum() and szINCHISource[pos_slash1+1] != '/':
            return None

        string = szINCHISource[self.LEN_INCHI_STRING_PREFIX:].strip()

        if not string:
            return None

        aux = string[(pos_slash1 - self.LEN_INCHI_STRING_PREFIX) + 1:]
        slen = len(aux)
        proto = False
        end = 0

        for idx, ch in enumerate(aux):
            if ch == '/':
                cn = aux[idx+1]
                if cn == 'c' or cn == 'h' or cn == 'q':
                    continue
                if cn == 'p':
                    proto = idx
                    continue
                if cn == 'f' or cn == 'r':
                    return None
                end = idx
                break
        else:
            end = slen

        #if end == (slen - 1):
        #    end += 1

        if not proto:
            smajor = aux[:end]
        else:
            smajor = aux[:proto]

        if proto:
            nprotons = int(aux[proto + 2:end])
            if nprotons > 0:
                if nprotons > 12:
                    flagproto = 'A'
                else:
                    flagproto = pplus[nprotons-1]

            elif nprotons < 0:
                if nprotons < -12:
                    flagproto = 'A'
                else:
                    flagproto = pminus[-nprotons-1]
            else:
                return None

        sminor = ''
        if end != slen:
            sminor = aux[end:]
        if len(sminor) < 255:
            sminor += sminor

        flag = flagstd if bStdFormat else flagnonstd

        digest_major = self.get_sha256(smajor)
        digest_minor = self.get_sha256(sminor)
        major = base26_triplet_1(digest_major) + base26_triplet_2(digest_major) + base26_triplet_3(digest_major) + \
                                    base26_triplet_4(digest_major) + base26_dublet_for_bits_56_to_64(digest_major)
        minor = base26_triplet_1(digest_minor) + base26_triplet_2(digest_minor) + \
                                    base26_dublet_for_bits_28_to_36(digest_minor)
        return "%s-%s%s%s-%s" % (major, minor, flag, flagver, flagproto)