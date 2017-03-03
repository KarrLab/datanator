Installation
============

*Dependencies*

Install Openbabel from <http://openbabel.org/wiki/Category:Installation>

*Install KineticDatanator*::

    $ pip install KineticDatanator

*Install the NCBI taxonomic database*::

    $ python -c 'from KineticDatanator import TaxonFinder; TaxonFinder.downloadNCBIDatabase()'
