``kinetic_datanator``
=====================
``kinetic_datanator`` is a software tool for finding experimental data for building and calibrating dynamical models of cellular biochemistry such as metabolite, RNA, and protein abundances; protein complex compositions; transcription factor binding motifs; and kinetic parameters. ``kinetic_datanator`` is particularly useful for building large models, such as whole-cell models, that require large amounts of data to constrain large numbers of parameters. ``kinetic_datanator`` was motivated by the need for large amounts of data to constrain whole-cell models and the fact that this data is hard to utilize because it is scattered across numerous siloed repositories.

``kinetic_datanator`` currently supports the following data types and data sources:

* Metabolite concentrations: `ECMDB <http://www.ecmdb.ca>`_ and `YMBD <http://www.ymdb.ca>`_
* RNA abundance: `ArrayExpress <https://www.ebi.ac.uk/arrayexpress>`_
* Protein abundance: `PaxDb <http://pax-db.org>`_
* Protein complex composition: `CORUM <http://mips.helmholtz-muenchen.de/corum>`_
* Transcription factor binding motifs: `JASPAR <http://jaspar.genereg.net>`_
* Reaction kinetics: `SABIO-RK <http://sabio.h-its.org>`_
* Taxonomy: `NCBI Taxonomy <https://www.ncbi.nlm.nih.gov/taxonomy>`_

``kinetic_datanator`` (1) downloads these repositories; (2) normalizes their data to a common ontology and units; (3) stores their data to a local SQLite database; and (4) provides a Python API for (a) finding relevant data to model a specific organism and environmental condition from similar species, reactions, genotypes (taxon, variant), and environments (temperature, pH, media), and (b) reducing multiple relevant observations to a single consensus recommended parameter value, and (c) exporting these consensus recommendations and their provenance to an Excel workbook. To make ``kinetic_datanator`` easier to use, we plan to develop user-friendly command line and web-based interfaces for finding data for SBML-encoded models.

``kinetic_datanator`` is under active development and is not yet ready for end users. Please check back soon for updates.

This website contains detailed documentation of the ``kinetic_datanator`` source code. Going forward, this website will also contain detailed instructions and tutorials on how to use ``kinetic_datanator``.

Contents
--------

.. toctree::
   :maxdepth: 3

   intro
   installation
   tutorial
   API documentation <source/modules.rst>
   about
