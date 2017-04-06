.. KineticDatanator documentation master file, created by
   sphinx-quickstart on Wed Feb 22 12:28:35 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to KineticDatanator's documentation!
============================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:




Introduction
============

Motivation
----------

How it Works
--------------

Databases Used
---------------

Installation
============

*Dependencies*

Install Openbabel from <https://packages.debian.org/wheezy/amd64/python-openbabel/download>

*Install KineticDatanator*::

    $ pip install KineticDatanator

*Install the NCBI taxonomic database*::

    $ python -c 'from KineticDatanator import TaxonFinder; TaxonFinder.downloadNCBIDatabase()'




Tutorial - Getting Familiar With KineticDatanator's Command Line Interface
============================================================================

Generate Template Doc
------------------------

KineticDatanator can be used run directly from the command line.

There are three required arguments necesssary to run the data collection method "find-kinetics"

1. A path to an excel sheet with the reactions you would like to search (we will discuss the necessary formatting later)
2. A name for the results file (the results are saved as an excel file, so the name should end in '.xlsx')
3. The name of the species you are searching for 

Let's do an example

Create a new directory

In that directory run the following command::

    $ python -mKineticDatanator generate-template

This will output a template document to your directory. This template is an example of an excel sheet that can be entered into KineticDatanator

In order to collect kinetic information on the reaction in the template, we will use "find-kinetics"
The first argument is a path to an excel sheet, we will use the template document we just generated labelled TheTemplateDocument.xlsx

The second argument is a name for the output file. We will use ExampleResults.xlsx, but you can name it anything. 

The third argument is the name of the species we are searching for. We will use 'homo sapiens'

Run "find-kinetics"::

    $ python -mKineticDatanator get-kinetics TheTemplateDocument.xlsx ExampleResults.xlsx 'homo sapiens'




$ python -mKineticDatanator get-taxonomic-lineage 'homo sapiens'


Set Maximum Proximity Limit
-----------------------------



Formatting The Reaction Entries
=================================



Understanding The Output
===========================










