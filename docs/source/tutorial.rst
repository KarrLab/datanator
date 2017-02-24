Tutorial - Getting Familiar With KineticDatanator's Command Line Interface
========================================================================

Generate Template Doc
--------------------

KineticDatanator can be used run directly from the command line.


Let's do an example

Create a new directory

In that directory run the following command::

    $ python -mKineticDatanator generate-template

This will output a template document to your directory. This template is an example of an excel sheet that can be entered into KineticDatanator

Collecting Kinetic Data
-------------------

In order to collect kinetic information on the reaction in the template, we will use "find-kinetics"

There are three required arguments necesssary to run the data collection method "find-kinetics"

1. A path to an excel sheet with the reactions you would like to search (we will discuss the necessary formatting later)
2. A name for the results file (the results are saved as an excel file, so the name should end in '.xlsx')
3. The name of the species you are searching for 

The first argument is a path to an excel sheet, we will use the template document we just generated labelled TheTemplateDocument.xlsx

The second argument is a name for the output file. We will use ExampleResults.xlsx, but you can name it anything. 

The third argument is the name of the species we are searching for. We will use 'homo sapiens'

Run "find-kinetics"::

    $ python -mKineticDatanator get-kinetics TheTemplateDocument.xlsx ExampleResults.xlsx 'homo sapiens'


$ python -mKineticDatanator get-taxonomic-lineage 'homo sapiens'


Entering Reactions Into KineticDatanator
---------------------------------------

Open up the template document TheTemplateDocument.xlsx


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

Each metabolite used in the stoichiometric string needs to be structurally defined in the "Metabolites" worksheet. Open up 
the "Metabolites" worksheet


+-----------------------------------------------------------------------------------------+
|Compound ID |Structure                                                    |
+=========================================================================================+
|ATP         |NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1
+-----------------------------------------------------------------------------------------+
|UMP         |OC1C(O)C(OC1COP([O-])([O-])=O)N1C=CC(=O)NC1=O
+-----------------------------------------------------------------------------------------+
|UDP         |OC1C(O)C(OC1COP([O-])(=O)OP([O-])([O-])=O)N1C=CC(=O)NC1=O
+-----------------------------------------------------------------------------------------+
|AMP         |NC1=C2N=CN(C3OC(COP([O-])([O-])=O)C(O)C3O)C2=NC=N1
+-----------------------------------------------------------------------------------------+
|GMP         |NC1=NC2=C(N=CN2C2OC(COP([O-])([O-])=O)C(O)C2O)C(=O)N1
+-----------------------------------------------------------------------------------------+
|GDP         |NC1=NC2=C(N=CN2C2OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C2O)C(=O)N1
+-----------------------------------------------------------------------------------------+




















Set Maximum Proximity Limit
--------------------------
