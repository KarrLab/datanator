Installation
============

*Dependencies*

Install Openbabel from <https://packages.debian.org/wheezy/amd64/python-openbabel/download>

*Install KineticDatanator*::

    $ pip install KineticDatanator

*Install the NCBI taxonomic database*::

    $ python -c 'from KineticDatanator import TaxonFinder; TaxonFinder.downloadNCBIDatabase()'
