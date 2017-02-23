Introduction
============

What Is KineticDatanator
------------------------


KineticDatanator provides a software and methodology for large scale kinetic data collection in a user friendly way.
The input is an excel sheet with reaction query information, and the output is an excel sheet with information about the Km and Vmax of the reation.

The methodology of the how KineticDatanator searches, as well as additional search setting options will be discussed later, but at the base, using KineticDatanator is very simple: Excel sheet in, excel sheet out.


Motivation
-----------

With recent advances in biological modelling, researchers increasingly rely upon kinetic data to set the kinetic 
parameters for their models. Sabio-RK, and other such databases, have much of this data. However, the data is often not found 
due to inconsistent cateloging or time constraints. 

A dynamical modeler is faced with the following tasks when searching for kinetic data to build a model for a particular species. 

1. The scientist must format a query string using the reaction information
2. From the expiremental entries that match the query, the scientist must narrow the entries to expirements done in the species being modelled, or to expirements done in organisms most closely related to the species being modelled. 
3. From the narrowed set of expirments, the scientiest must use a statistical method to decide which value is the most meaningful.
4. If, after step 1, no resulsts showed up, then the scientist must look for similar reactions, and then repeat steps 2 and 3

And all this for a only a single reaction!

KineticDatanator solves, and automates all these steps. 



Step 1 - Formulating a Reaction Query
-------------------------------------

*Why This is Difficult*

One approach is to search a reaction using external identifiers (KEGG Reaction ID, KEGG Compound ID, PubChemID, etc).
However, a modeller may not know this information, and aggregating it for each reaction or compound may be too time-consuming.

Another second approach is to search for a reaction using canonical compound names for the subtrates and products. However, canonical compound names are more of a myth. Sure, everyone knows what ATP is, but is it “N-Methylanthranilic acid”  or “2-(Methylamino)benzoic acid”. For many compounds, no canonical names exist. 

A third approach is to search using molecular/structural information (ie. molecular formula, inchi string, SMILES).One could formulate a search using molecular formula, InchiString, or SMILES. However, each of these is problematic. The problem with using molecular formula is that it contains no structural information, and therefore identifies more compounds than it should. For example, ethanol and dimethyl ether are both C2H60, but are different molecules. The problem with using Inchi or SMILES is that both are sensitive to the particular pH of a molecule because the structural (number of attached hydrogens) changes based on ph. For example,……..


*How KineticDatanator Solves This*



Step 2 - Narrowing Expiremental Results By Taxonomic Proximity  
---------------------------------------------------------------

*Why This Is Hard*

This problem is simply way too time consuming. No one has time to look through 102 entries and catalogue them by relationship. Do you know offhand if E. Coli is more closely related to Asticcacaulis sp. AC402 or to Bryobacter aggregatus MPL3???? We certainly don’t. 


*How KineticDatanator Solves This*

KineticDatanator uses the NCBI taxonomic tree to judge proximity of species. When a reaction query is made, KineticDatanator assigns a “proximity” to each entry. The proximity corresponds to the number of nodes one needs to go up on the taxonomic tree to get to the  lowest common node between the model’s species, and the experimental entries species. For example…


Step 3 - Finding Most Relevant Values Among Entries
-----------------------------------------------------

*Why This Is Hard*

Again, this is simply time consuming

*How Kinetic Datanator Solves This*

For both Vmax and Km, Kinetic Datanator finds the median value among the entries, and records it in the excel output.





Databases Used
------------

Sabio. Kegg. NCBI
