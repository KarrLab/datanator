"""
Downloads and parse the latest release for 
ortholog -> uniprot mappings from omabrowser.org
"""

import urllib2
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm

mapping = urllib2.urlopen("https://omabrowser.org/All/oma-uniprot.txt.gz")
