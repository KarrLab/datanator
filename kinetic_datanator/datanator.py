""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import io
from . import sabio_rk
from .core import data_structs
from .data_source import ezyme
from .util import reaction_util


class Datanator(object):
    """
    Attributes:
        max_taxon_dist (:obj:`float`, optional): maximum taxonomic distance from the target taxon to its latest common ancestor with the observed taxon
        include_mutants (:obj:`bool`, optional): if :obj:`True`, include observations from mutants
        min_temp (:obj:`float`, optional): minimum observed temperature
        max_temp (:obj:`float`, optional): maximum observed temperature
        min_ph (:obj:`float`, optional): minimum observed pH
        max_ph (:obj:`float`, optional): maximum observed pH
    """

    def __init__(self, max_taxon_dist=None, include_mutants=False,
                 min_temp=None, max_temp=None, min_ph=None, max_ph=None):
        """
        Args:
            max_taxon_dist (:obj:`float`, optional): maximum taxonomic distance from the target taxon to its latest common ancestor with the observed taxon
            include_mutants (:obj:`bool`, optional): if :obj:`True`, include observations from mutants
            min_temp (:obj:`float`, optional): minimum observed temperature
            max_temp (:obj:`float`, optional): maximum observed temperature
            min_ph (:obj:`float`, optional): minimum observed pH
            max_ph (:obj:`float`, optional): maximum observed pH        
        """
        self.max_taxon_dist = max_taxon_dist
        self.include_mutants = include_mutants
        self.min_temp = min_temp
        self.max_temp = max_temp
        self.min_ph = min_ph
        self.max_ph = max_ph

    def run(self, input_filename, output_filename=''):
        """
        Args:
            input_filename (:obj:`str`): path to an Excel workbook which contains a list of reactions to find kinetic data for
            output_filename (:obj:`str`, optional): path to save the retrieved kinetic data as an Excel workbook

        Returns:
            :obj:``
        """
        # read in preliminary model definition (compartments, molecules, reactions)
        model = self.read_model(input_filename)

        # annotate model
        self.annotate_model(model)

        # get relevant data for model
        data = self.get_data(model)

        # save data to a file
        if output_filename:
            self.save_data(model, data, output_filename)

        # return results
        return data

    def read_model(self, filename):
        """
        Args:
            filename (:obj:`str`): path to an Excel workbook which contains a list of reactions to find kinetic data for
        """
        return io.InputReader.run(filename)

    def annotate_model(self, model):
        """ Annotate a model

        Args:
            model
        """
        self.annotate_molecules(model)
        self.annotate_reactions(model)

    def annotate_molecules(self, model):
        """ Annotate the molecules of a model

        * Cross-reference molecule to SABIO-RK compound IDs and names

        Args:
            model
        """
        taxon, compartments, molecules, reactions = model

        sabio_db = sabio_rk.SabioRkUtil()
        for mol in molecules:
            for sabio_compound in sabio_db.get_compounds_by_structure(mol.to_inchi()):
                mol.cross_references.append(data_structs.CrossReference(source='sabio-id', id=sabio_compound.id))
                mol.cross_references.append(data_structs.CrossReference(source='sabio-name', id=sabio_compound.name))

    def annotate_reactions(self, model):
        """ Annotate the reactions of a model

        * Predict missing EC numbers

        Args:
            model
        """
        taxon, compartments, molecules, reactions = model

        for rxn in reactions:
            if not rxn.get_ec_number():
                for result in ezyme.Ezyme.run(rxn):
                    rxn.cross_references.append(
                        data_structs.CrossReference(source='ec-code',
                                                    id=result.ec_number,
                                                    relevance=result.score,
                                                    assignment_method=data_structs.CrossReferenceAssignmentMethod.predicted))

    def get_data(self, model):
        """
        Args:
            model

        Returns:
            :obj:``
        """
        return self.get_sabio_data()

    def get_sabio_data(self, model):
        """
        Args:
            model

        Returns:
            :obj:``
        """
        taxon, compartments, molecules, reactions = model

        results = []
        for rxn in reactions:
            results.append(sabio_rk.Query(rxn, taxon, self.max_taxon_dist, self.include_mutants,
                                                 self.min_temp, self.max_temp, self.min_ph, self.max_ph).run())
        return results

    def save_data(self, model, data, filename):
        """ Save 

        Args:
            model
            data
            filename (:obj:`str`): path to save the retrieved kinetic data as an Excel workbook
        """
        taxon, compartments, molecules, reactions = model
        io.ResultsWriter.run(taxon, data, filename)
