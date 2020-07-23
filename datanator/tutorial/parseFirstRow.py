import json

class ParsePsort:

    def __init__(self, max_entries=float('inf')):
        self.max_entries = max_entries

    def parse_psortdb(self):
        """To parse datbase psortdb

        Args:
            max_entries(:obj:`int`): description of what it is

        Return:
            ()
        """
        x = {"a": 1, "b": 2}
        return json.dumps(x)
