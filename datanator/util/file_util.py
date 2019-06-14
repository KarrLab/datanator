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