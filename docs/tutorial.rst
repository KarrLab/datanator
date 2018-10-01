Tutorial - Getting Familiar With datanator's Command Line Interface
===========================================================================

Generate Template Doc
---------------------

datanator can be used run directly from the command line.


Let's do an example

Create a new directory

In that directory run the following command::

    $ datanator generate-template

This will output a template document to your directory. This template is an example of an excel sheet that can be entered into datanator

Collecting Kinetic Data
-----------------------

In order to collect kinetic information on the reaction in the template, we will use "find-kinetics"

There are three required arguments necesssary to run the data collection method "find-kinetics"

1. A path to an excel sheet with the reactions you would like to search (we will discuss the necessary formatting later)
2. A name for the results file (the results are saved as an excel file, so the name should end in '.xlsx')
3. The name of the species you are searching for 

The first argument is a path to an excel sheet, we will use the template document we just generated labelled InputTemplate.xlsx

The second argument is a name for the output file. We will use `Results.xlsx`, but you can name it anything. 

The third argument is the name of the species we are searching for. We will use 'homo sapiens'

Run "get-kinetics"::

    $ datanator get-kinetics InputTemplate.xlsx Results.xlsx 'homo sapiens'


$ datanator get-taxonomic-lineage 'Escherichia coli'


Entering Reactions Into datanator
-----------------------------------------

Open up the template document InputTemplate.xlsx


This document has two worksheets: Reactions, and Metabolites.

*Reactions Worksheet*

+------------------------------------------------+-------------------------+
|Reaction ID                                     |Stoichimetry             |
+================================================+=========================+
|A Reacion Name #1 (ATP + UMP <==> UDP + ADP)    |ATP + UMP <==> UDP + ADP |
+------------------------------------------------+-------------------------+
|A Reaction Name #2 (GMP + ATP <==> ADP + GDP)   |GMP + ATP ==> ADP + GDP  |
+------------------------------------------------+-------------------------+


Look at the Reactions worksheet. There are two columns. The first is an identifier for the reaction. You can name this whatever you
like as long as it is unique. A good suggestion is just to use the reaction string (in the template document, we did not follow this practice in order to illustrate that the name can be anything). The second column is a the stoichiometric string. 

Each metabolite used in the stoichiometric string needs to be structurally defined in the "Metabolites" worksheet. Open up the "Metabolites" worksheet


+------------+----------------------------------------------------------------------------+
|Compound ID |Structure                                                                   |
+============+============================================================================+
|ATP         |NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1  |
+------------+----------------------------------------------------------------------------+
|UMP         |OC1C(O)C(OC1COP([O-])([O-])=O)N1C=CC(=O)NC1=O                               |              
+------------+----------------------------------------------------------------------------+
|UDP         |OC1C(O)C(OC1COP([O-])(=O)OP([O-])([O-])=O)N1C=CC(=O)NC1=O                   |
+------------+----------------------------------------------------------------------------+
|AMP         |NC1=C2N=CN(C3OC(COP([O-])([O-])=O)C(O)C3O)C2=NC=N1                          | 
+------------+----------------------------------------------------------------------------+
|GMP         |NC1=NC2=C(N=CN2C2OC(COP([O-])([O-])=O)C(O)C2O)C(=O)N1                       |  
+------------+----------------------------------------------------------------------------+
|GDP         |NC1=NC2=C(N=CN2C2OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C2O)C(=O)N1           | 
+------------+----------------------------------------------------------------------------+


The compound ID corresponds to the ID of the compound in the stoichiometric string. The structure is 
either an Inchi or a SMILES string (SMILES in this example). 

Let's try adding a reaction:

Lets say we wanted to add the reaction::

   Adenosine 3',5'-bisphosphate + H2O ==> phosphate + AMP 

The first step is to give each metabolite in the stoichiometric string and ID without spaces. So let's 
change it to::

    A-3-5-bisphosphate + H2O ==> phosphate + AMP

The second step is to add this reaction string to the second column in "Reactions" worksheet. 
In the first column, give the reaction some distinct name.


+------------------------------------------------+-----------------------------------------------------+
|Reaction ID                                     |Stoichimetry                                         |
+================================================+=====================================================+
|A Reacion Name #1 (ATP + UMP <==> UDP + ADP)    |ATP + UMP <==> UDP + ADP                             |
+------------------------------------------------+-----------------------------------------------------+
|A Reaction Name #2 (GMP + ATP <==> ADP + GDP)   |GMP + ATP ==> ADP + GDP                              |
+------------------------------------------------+-----------------------------------------------------+
|A distinct name of your choosing                |A-3-5-bisphosphate + H2O ==> phosphate + AMP         |
+------------------------------------------------+-----------------------------------------------------+


The third step is to structurally define each metabolite in the reaction string. We already defined structurally defined AMP (we used it in the previous reactions), so we will have to structurally define A-3-5-bisphosphate, H2O, and phosphate. 

The structural information is here::
    
    A-3-5-bisphosphate - NC1=C2N=CN(C3OC(COP([O-])([O-])=O)C(OP([O-])([O-])=O)C3O)C2=NC=N1
    
    H2O - O

    phosphate - OP([O-])([O-])=O


Now we need to add this information to the "Metabolites" worksheet.

Open up the "Metabolites" worksheet. Add the name of the compound used in the stoichiometric string (ex: phosphate) to the first column, and add the structure in the second. 

+-------------------+----------------------------------------------------------------------------+
|Compound ID        |Structure                                                                   |
+===================+============================================================================+
|ATP                |NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1  |
+-------------------+----------------------------------------------------------------------------+
|UMP                |OC1C(O)C(OC1COP([O-])([O-])=O)N1C=CC(=O)NC1=O                               |              
+-------------------+----------------------------------------------------------------------------+
|UDP                |OC1C(O)C(OC1COP([O-])(=O)OP([O-])([O-])=O)N1C=CC(=O)NC1=O                   |
+-------------------+----------------------------------------------------------------------------+
|AMP                |NC1=C2N=CN(C3OC(COP([O-])([O-])=O)C(O)C3O)C2=NC=N1                          | 
+-------------------+----------------------------------------------------------------------------+
|GMP                |NC1=NC2=C(N=CN2C2OC(COP([O-])([O-])=O)C(O)C2O)C(=O)N1                       |  
+-------------------+----------------------------------------------------------------------------+
|GDP                |NC1=NC2=C(N=CN2C2OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C2O)C(=O)N1           | 
+-------------------+----------------------------------------------------------------------------+
|A-3-5-bisphosphate |NC1=C2N=CN(C3OC(COP([O-])([O-])=O)C(OP([O-])([O-])=O)C3O)C2=NC=N1           | 
+-------------------+----------------------------------------------------------------------------+
|H2O                |O                                                                           |  
+-------------------+----------------------------------------------------------------------------+
|phosphate          |OP([O-])([O-])=O                                                            | 
+-------------------+----------------------------------------------------------------------------+


Once again, run "get-kinetics"::

    $ datanator get-kinetics InputTemplate.xlsx Results.xlsx 'Escherichia coli'


You should see results for the new reaction you inputted. 


Set Maximum Proximity Limit
------------------------------

The ideal kinetic data is information about the reaction being studied, collected from an expirement done in the species you are studying. However, often you will have to rely on kinetic data from different species. At a certain taxonomic distance, you might decide it's better to collect data from similar reactions taken from expirements in more closely related organisms. 

There are two dimensions of granularity - reaction variation and species variation - and the user can decide which data is preferred. 

The user can do this by setting the "proximLimit"

Let's try an example:

First, you want to get taxonomic infomation about the organism you are stuyding. Run::

    $ datanator get-taxonomic-lineage 'Escherichia coli'

You should see::

    1: Escherichia coli
    2: Escherichia
    3: Enterobacteriaceae
    4: Enterobacterales
    5: Gammaproteobacteria
    6: Proteobacteria
    7: Bacteria
    8: cellular organisms
    9: root

This is the number of nodes as you start from your organism, and climb up to the top of the taxonomic tree. Each number corresponds to a node. 

So let's say you are studying Escherichia coli. Maybe you think that anything outside the phylum protobacteria 
is too distantly related to be useful. In that case, you will want to run the "get-kinetics" argument, with an 
optional argument --proximit-limit. The number given after --proxim-limit is the highest node that will be considered useful. Since we have chosen Protobacteria, that number is 6

So, run::

    datanator get-kinetics InputTemplate.xlsx Results.xlsx 'Escherichia coli' --proxim-limit 6
