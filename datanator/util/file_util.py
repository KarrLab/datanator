from itertools import chain

class FileUtil:

    def extract_values(self, obj, key):
        """Pull all values of specified key from nested JSON.
        """
        arr = []

        def extract(obj, arr, key):
            """Recursively search for values of key in JSON tree."""
            if isinstance(obj, dict):
                for k, v in obj.items():
                    if isinstance(v, (dict, list)):
                        extract(v, arr, key)
                    elif k == key:
                        arr.append(v)
            elif isinstance(obj, list):
                for item in obj:
                    extract(item, arr, key)
            return arr

        results = extract(obj, arr, key)

        return results

    def flatten_json(self, nested_json):
        '''
            Flatten json object with nested keys into a single level.
            e.g. 
            {a: b,                      {a: b,  
             c: [                        d: e,
                {d: e},    =>            f: g }
                {f: g}]}
            Args:
                nested_json: A nested json object.
            Returns:
                The flattened json object if successful, None otherwise.
        '''
        out = {}

        def flatten(x, name=''):
            if type(x) is dict:
                for a in x:
                    flatten(x[a], name + a + '_')
            elif type(x) is list:
                i = 0
                for a in x:
                    flatten(a, name + str(i) + '_')
                    i += 1
            else:
                out[name[:-1]] = x

        flatten(nested_json)
        return out

    def unpack_list(self, _list):
        ''' Unpack sublists in a list
            Args:
                _list: a list containing sublists  e.g. [ [...], [...], ...  ]
            Return:
                result: unpacked list e.g. [ ....  ]
        '''
        return list(chain.from_iterable(_list))

    def access_dict_by_index(self, _dict, count):
        ''' Assuming dict has an order, return 
            the first num of elements in dictionary
            Args:
                _dict: { 'a':1, 'b':2, 'c':3, ... }
                count: number of items to return
            Return:
                result: a dictionary with the first count 
                        from _dict
                        {'a':1}
        '''
        result = {}
        tuples = _dict.items()
        i = 0
        for item in tuples:
            if i == count:
                continue
            result[item[0]] = item[1]
            i += 1
        return result

    def replace_dict_key(self, _dict, replacements):
        ''' Replace keys in a dictionary with the order
            in replacements e.g.,
            {'a': 0, 'b': 1, 'c': 2}, ['d', 'e', 'f'] =>
            {'d': 0, 'e': 1, 'f': 2}            
            Args:
                _dict: dictionary whose keys are to be replaced
                replacement: list of replacement keys

            Return:
                result: dictionary with replaced keys
        '''
        result = {}
        i = 0

        for k, v in _dict.items():
            result[replacements[i]] = v
            i += 1
        return result

    def get_common(self, list1, list2):
        ''' Given two lists, find the closest
            common ancestor
            Args:
                list1: [a, b, c, f, g] 
                list2: [a, b, d, e]
            Return:
                result: the closest common ancestor, in
                        the above example would be b
        '''
        ancestor = ''
        for a, b in zip(list1, list2):
            if a == b:
                ancestor = a
            else:
                return ancestor
        return ancestor

    def make_dict(self, keys, values):
        ''' Give two lists, make a list of 
            dictionaries
            Args:
                keys: [a, b, c, d, ...]
                values: [1, 2, 3, 4]
            Return:
                dic: {'a': 1, 'b': 2, 'c': 3, ...} 
        '''
        result = {}
        for k, v in zip(keys, values):
            result[k] = v
        return result

    def search_dict_list(self, dict_list, key, value):
        ''' Find the dictionary with 
            key/value pair in a list of dictionaries

            Args:
                dict_list (:obj: `list`): list of dictionaries
                key (:obj: `string`): key in the dictionary
                value (:obj: `string`): value to be matched
            Returns:
                result (:obj: `dictionary`): list of dictionaries with the key/value pair
        '''
        return list(filter(lambda search: search[key] == value, dict_list))