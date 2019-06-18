import hashlib

class ChemUtil:

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