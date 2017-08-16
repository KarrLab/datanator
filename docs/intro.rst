Overview
========
``kinetic_datanator`` is a software tool for finding experimental data for building and calibrating dynamical models of cellular biochemistry such as metabolite, RNA, and protein abundances; protein complex compositions; transcription factor binding motifs; and kinetic parameters. ``kinetic_datanator`` is particularly useful for building large models such as whole-cell models that require large amounts of data to constrain large numbers of parameters.


Motivation
----------
Large models such as whole-cell models are needed to help researchers, bioengineers, and physicians predict how genotype and the environment determine phenotype. In particular, large models that represent the function of each individual gene are needed to help bioengineers rationally design microorganisms to perform specific functions such as producing drugs and neutralizing pathogens and to help physicians interpret personal genomes and design personalized therapies tailored to each patient's unique genome.

Despite their potential, large models are hard to build, in part, because it is hard to aggregate the large amount of data needed to build large models because this data is scattered across numerous siloed repositories. Without data aggregation tools such as ``kinetic_datanator``, this data must be manually aggregated:

1. Modelers must identify appropriate data sources for their model.
2. Modelers must identify the relevant subset of this data for their model.
3. Modelers must merge inconsistently annotated data from multiple sources.
4. Modelers must reduce multiple experimental observations to individual consensus parameter value recommendations.
5. Modelers must record the provenance of each recommendation so their model is comprehensible and reproducible.

This requires extensive time and effort, is unscalable, introduces substantial selection bias into the data underlying models, and is irreproducible.

To address these problems, we have developed ``kinetic_datanator`` to systemize and accelerate the aggregation of data for biomodels.


Methodology
----------- 
1. ``kinetic_datanator`` downloads data from the data repositories listed below.
2. ``kinetic_datanator`` parses this data, normalizes this data to a common schema and common identifiers, and converts this data to common units.
3. ``kinetic_datanator`` stores the normalized data and its provenance to a local SQLite database.
4. ``kinetic_datanator`` helps researchers find relevant data to model a specific organism and environmental condition from similar species, reactions, genotypes, and environments according to the following filters:

    * Chemical similarity: ``kinetic_datanator`` helps researchers identify data observed for (a) molecularly-similar species as
      determined by the Tanimoto similarity of their molecular fingerprints and (b) chemically-similar reactions that involve similar
      reaction centers and reaction mechanisms.
    * Taxonomic similarity: ``kinetic_datanator`` helps researchers identify data observed for taxonomically-close taxa and similar
      genetic variants.
    * Environmental similarity: ``kinetic_datanator`` helps researchers identify data observed for similar temperatures, pHs, and growth
      media.

5. ``kinetic_datanator`` reduces multiple relevant observations to a single consensus recommended parameter value.
6. ``kinetic_datanator`` records the provenance of the data underlying each consensus recommendation.
7. ``kinetic_datanator`` exports the consensus recommendations and their provenance to Excel workbooks.

Currently, ``kinetic_datanator`` provides researchers a Python API to programmatically aggregate data for models. To make ``kinetic_datanator`` easier to use, we plan to develop user-friendly command line and web-based interfaces which will help users find data for SBML-encoded models, review this data, and generate consensus recommendations for parameter values.


Supported data types and data sources
--------------------------------------
``kinetic_datanator`` currently supports the following data type and data sources. We are also actively integrating additional data types and data sources into ``kinetic_datanator``.

* Databases

    * Metabolite concentrations: `ECMDB <http://www.ecmdb.ca>`_ and `YMBD <http://www.ymdb.ca>`_
    * RNA abundance: `ArrayExpress <https://www.ebi.ac.uk/arrayexpress>`_
    * Protein abundance: `PaxDb <http://pax-db.org>`_
    * Protein complex composition: `CORUM <http://mips.helmholtz-muenchen.de/corum>`_
    * Transcription factor binding motifs: `JASPAR <http://jaspar.genereg.net>`_
    * Reaction kinetics: `SABIO-RK <http://sabio.h-its.org>`_
    * Taxonomy: `NCBI Taxonomy <https://www.ncbi.nlm.nih.gov/taxonomy>`_

* Prediction tools

    * Metabolite properties (charge, pK\ :subscript:`a`, protonation): `Open Babel <http://openbabel.org>`_
    * EC number: `E-zyme <http://www.genome.jp/tools/e-zyme>`_
    

Advantages of ``kinetic_datanator`` versus manual data aggregation
------------------------------------------------------------------
* ``kinetic_datanator`` accelerates data aggregation by merging data from multiple repositories into a single location.
* ``kinetic_datanator`` reduces selection bias in the data used to build models by systemizing how researchers find data for models.
* ``kinetic_datanator`` helps researchers find more data for models by helping researchers find data from chemically-similar species and reactions, taxonomically-similar organisms, and similar environmental conditions based on their temperature, pH, and media.
* ``kinetic_datanator`` increases the comprehensibility and reproducibility of models by automatically tracking the provenance of each recommended parameter value.
