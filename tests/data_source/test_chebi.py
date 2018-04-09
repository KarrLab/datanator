""" Tests of Chebi
:Author:  Saahith Pochiraju  <saahith116@gmail.com>
:Date: 2018-04-09
:Copyright: 2018, Karr Lab
:License: MIT
"""

from kinetic_datanator.data_source import chebi
import datetime
import dateutil
import os
import shutil
import tempfile
import unittest


class TestChebiDownload(unittest.TestCase):

    def test_download(self):
        pass
        # chub = chebi.Chebi(download_backups=False, load_content=True, clear_content=False)
