Overview
================

What is ``kinetic_datanator``
-------------------------

``kinetic_datanator`` is a software tool that helps biomodelers find relevant experimental data for building and calibating dynamical models of cellular biochemistry. ``kinetic_datanator`` is particularly useful for large models, such as whole-cell models, that require large amounts of data to constrain large numbers of parameters.

A user can query kinetic_datanator with a list of biological parameters (e.g. the intracellular concentration of a ATP, ADP, and AMP) and a list of filters tailored to the model’s conditions (e.g. Organism: Bacillus subtilis, pH: 6.7, Temp: 23 °C). For each parameter, kinetic_datanator outputs a consensus value and the underlying data or prediction tools used to calculate the consensus.  

Motivation
----------

Many biological applications, such as personalized therapy and rational microbial engineering, require large-scale models. Construction of these large-scale models requires extensive heterogeneous data. The modeler is tasked with:

1. aggregating and merging heterogeneous data
2. identifying the most relevant subset of data for his model conditions
3. selecting a single value from the experimental distribution to describe each parameter
4. recording the provenance of the aggregated data. 

Performing the data-collection process manually is slow. The need to cross reference distinct databases further slows data collection. For example, a modeler must consult a taxonomic database to identify the relatedness of experimental organisms in a reaction kinetics database. This time cost restricts the possible size of biomodels. Furthermore, provenance is tedious to record for manually aggregated data. 

To solve these problems, kinetic_datanator automates the data-collection process and records provenance.

kinetic_datanator's data-selection process
-------------------------------------------

**1. Aggregating and merging heterogenous data**

Kinetic_datanator integrates heterogeneous data types by downloading multiple data sources to a central database, unifying the biological definitions, and normalizing the units. 

**2. Identifying the most relevant data for a given model**

Kinetic_datanator filters to the most relevant data data based on the user’s specifications. 

Possible filters include: 

a. biological similarity - identify experiments with organisms most closely related to the model’s organism 
b. chemical similarity - identify experiments with the most similar reaction mechanism 
c. environmental similarity - identify the experiments with the closest pH, temperature, and growth media to the model’s environmental conditions. 

The user can control how the filters are weighted.

**3. Reducing the distribution of data to a single consensus value**

Kinetic_datanator reduces the distribution of data to a single consensus value. The user has the option to specify how the distribution should be reduced. 

**4. Recording the provenance of the data**

For each parameter value, kinetic_datanator records the underlying data-sources and recommended consensus value. The user has the option to review the recommendations and modify them. Kinetic_datanator also records the entirety of the query to allow reproducible and updateable model construction.  


Advantages to using kinetic_datanator:
--------------------------------------
1. Kinetic_datanator enables the creation of large-scale models by accelerating the data-collection process 
2. kinetic_datanator automatically tracks provenance
3. Models queried in kinetic_datanator are able to be easily reproduced with the same data or updated to incorporate new experiments. 
4. Models queried in kinetic_datanator are able to be easily studied under different conditions (e.g. analyzing the same biological parameters in a different organism).
5. By making the parameter selection process explicit,  we hope a methodology on parameter selection can emerge.

Data sources currently in kinetic_datanator
-----------------------------------------------------

**Databases**

* Metabolite concentrations: `ECMDB <http://www.ecmdb.ca>`_ and `YMBD <http://www.ymdb.ca>`_
* RNA abundance: `ArrayExpress <https://www.ebi.ac.uk/arrayexpress>`_
* Protein abundance: `PaxDb <http://pax-db.org>`_
* Protein complex composition: `CORUM <http://mips.helmholtz-muenchen.de/corum>`_
* Transcription factor binding motifs: `JASPAR <http://jaspar.genereg.net>`_
* Reaction kinetics: `SABIO-RK <http://sabio.h-its.org>`_
* Taxonomy: `NCBI Taxonomy <https://www.ncbi.nlm.nih.gov/taxonomy>`_

**Prediction tools**

* Protonation: `Open Babel <http://openbabel.org>`_
* EC number: `E-zyme <http://www.genome.jp/tools/e-zyme>`_
