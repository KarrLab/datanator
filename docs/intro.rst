Introduction
================

What Is kinetic_datanator
-------------------------

**Target Audience**

Scientists who would like to aggregate kinetic data (Vmax and Km) for a particular species. kinetic_datanator is especially useful 
for dynamical modelers as it allows for large-scale data collection tailored to the species being modelled. 

**What Does kinetic_datanator Do?**

kinetic_datanator takes inputted reaction information, searches through Sabio-RK's kinetic database, and returns the most relevant kinetic results. 

**What Are The Input and Outputs?**

The input is an excel sheet containing rection information. The output is an excel sheet containing kinetic information associated with each reaction. 

**What will I Need To Generate The Input File?**

In order to generate the input file, you will need two sets of information:

1. A stoichiometric string for each reaction
2. Structural information (either an Inchi or SMILES string) for each metabolite that participates in your reactions

(I should site some databases where people can look these up)




Motivation
----------

Dynamical modelers rely on expiremental data to set kinetic parameters. Becuase these biological models are species specific, the expiremental entries must be tailored to the particular species being modelled. Current kinetic databases do not provide modellers with search tools to tailor the kinetic data. Furthermore, modellers must somehow find the most relevant kinetic value from a varied list of expiremental values.

Kinetic Datanator seeks to overcome these problems by providing scientists three tools:

1. Software to automate the search process
2. A default methodolgy for parsing through large quantities of kinetic data to create customized "species-specific" kinetic values
3. Tools for the scientist to define his own search methodology

We hope that by automating and systematizing the search process, scientists can design biological models more efficiently and accurately.



Searching Methodology
---------------------

kinetic_datanator creates customised kinetic data in four steps

1. A query string is generated from the reaction information and used to search Sabio-RK's reaction database
2. The expiremental results found in Sabio-RK are narrowed to include only the expiremental entries most closely related to the species being modelled.
3. From the narrowed set of entries, the most relevant kinetic values are highlighted to the user.
4. If after step 1 no results were found, then kinetic_datanator searches for similar reactions, and repeates steps 2 and 3



Step 1 - Formulating a Reaction Query
-------------------------------------

kinetic_datanator searches for kinetic data using the structural information (inchi or SMILES) of the metabolites involved in a reaction. In general, searching with structural information is difficult becuase varying protonation states can differentiate otherwise identical molecules. kinetic_datanator solves this by exploiting a unique characteristic of Inchi strings: all of their protonation information is stored at the end and is clearly delinated with a "/h" tag. Removal of the protonation information allows for the creation "generic inchi's" while keeping the necessary structural information.

For example, here are two potential Inchi strings for ATP::

    InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1
    
    InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-4/t4-,6-,7-,10-/m1/s1

Notice how the end of the strings vary. 

However, if all the protonation information is removed, the two Inchi strings are identical::

    InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20
    
    InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20

kinetic_datanator can be used with Inchi or SMILES. If metabolites are defined with SMILES, the program will simply convert it to Inchi string. 


Step 2 - Narrowing Expiremental Results By Taxonomic Proximity  
--------------------------------------------------------------


kinetic_datanator uses the NCBI taxonomic tree to judge proximity of species. When a reaction query is made, kinetic_datanator assigns a “proximity” to each entry. The proximity corresponds to the number of nodes one needs to go up on the taxonomic tree to get to the  lowest common node between the model’s species, and the experimental entries species. For example…


Step 3 - Finding Most Relevant Values Among Entries
----------------------------------------------------





Databases Used
--------------

Sabio. Kegg. NCBI
