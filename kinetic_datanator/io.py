""" Utilities for reading and writing data to/from files

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-14
:Copyright: 2017, Karr Lab
:License: MIT
"""

from . import data_structs
from . import inchi_generator
from .util import taxonomy_util
from wc_utils.workbook.core import Row, Workbook, Worksheet
from wc_utils.workbook.io import WorkbookStyle, WorksheetStyle, write
import openpyxl


class InputReader(object):

    @classmethod
    def read_compounds(cls, filename):
        """
        Args:
            filename

        Returns:
            :obj:`list` of `Compound`: list of compounds
        """

        sabio_name_to_inchi_dict = inchi_generator.getSabioNameToInchiDict()
        wb = openpyxl.load_workbook(filename=filename)
        ws = wb.get_sheet_by_name('Metabolites')

        # instantiate a compound object for each metabolite in the excel sheet
        compound_list = []
        for i in range(2, ws.max_row + 1):
            id = ws.cell(row=i, column=1).value

            smiles_or_inchi = ""
            structure = ws.cell(row=i, column=2).value
            sabio_names = []
            if structure != None:
                generic_inchi = inchi_generator.generateGenericInchi(structure)
                for name in sabio_name_to_inchi_dict:
                    if sabio_name_to_inchi_dict[name] == generic_inchi:
                        sabio_names.append(name)
                smiles_or_inchi = structure

            comp = data_structs.Compound(id, smiles_or_inchi, sabio_names)
            compound_list.append(comp)

        return compound_list
    
    @classmethod
    def read_reactions(cls, filename):
        wb = openpyxl.load_workbook(filename=filename)
        ws = wb.get_sheet_by_name('Reactions')
        reactions = []
        for i in range(2, ws.max_row + 1):
            reactions.append({
                'id': ws.cell(row=i, column=1).value, 
                'stoichiometry': cls.parse_reaction_stoichiometry(ws.cell(row=i, column=2).value),
                })
        return reactions

    @classmethod
    def parse_reaction_stoichiometry(cls, reaction_string):
        # takes in a search string and parses it into two arrays.
        # one array for substrate IDs, and one array for product IDs

        balancedMetab = []
        substrates = []
        products = []

        if "<==>" in reaction_string:
            bothSides = reaction_string.split("<==>")
        elif "==>" in reaction_string:
            bothSides = reaction_string.split("==>")
        elif "-->" in reaction_string:
            bothSides = reaction_string.split("-->")
        elif "<->" in reaction_string:
            bothSides = reaction_string.split("<->")
        elif "<=>" in reaction_string:
            bothSides = reaction_string.split("<=>")
        elif "=" in reaction_string:
            bothSides = reaction_string.split("=")

        substrates = bothSides[0]
        products = bothSides[1]

        parsedSubstrates = []
        parsedProducts = []

        # IMPORTANT!!!!! we are getting rid of hydrogens here.
        for entry in substrates.split(" "):
            if not (cls.is_number(entry)) and entry != "[c]:" and entry != "H" and entry != "H[c]" and entry != "h_m" and entry != "h_x" and entry != "h_c" and entry != "H[e]" and entry != "[m]:" and entry != "[e]:"and entry != "(2)" and entry != "+" and entry != "==>" and entry != "":
                # get rid of the [c] tag on some molecules
                if entry.find("[") != -1:
                    entry = entry[:entry.find("[")]
                parsedSubstrates.append(entry)
        for entry in products.split(" "):
            if not (cls.is_number(entry)) and entry != "[c]:" and entry != "H" and entry != "H[c]" and entry != "h_m" and entry != "h_x" and entry != "h_c" and entry != "H[e]" and entry != "[m]:" and entry != "[e]:"and entry != "(2)" and entry != "+" and entry != "==>" and entry != "":
                # get rid of the [c] tag on some molecules
                if entry.find("[") != -1:
                    entry = entry[:entry.find("[")]
                parsedProducts.append(entry)

        balancedMetab.append(parsedSubstrates)
        balancedMetab.append(parsedProducts)

        return balancedMetab

    @classmethod
    def is_number(cls, s):
        try:
            float(s)
            return True
        except ValueError:
            return False


class ResultsWriter(object):
    """ Save reaction kinetic data to an Excel workbook """

    @classmethod
    def run(cls, taxon, reactions, filename):
        """ Save a list of reaction kinetic data to an Excel workbook or set of csv/tsv files

        Args:
            taxon (:obj:`str`): taxon
            reactions (:obj:`list` of :obj:`kinetic_datanator.datanator.FormattedData`): list of reactions and their kinetic data
            filename (:obj:`str`): filename to store the list reaction kinetic data
        """
        wb = Workbook()
        style = WorkbookStyle()

        # taxon
        style['Taxon'] = WorksheetStyle(head_rows=1, head_columns=0,
                                        head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)

        ws = wb['Taxon'] = Worksheet()
        ws.append(Row(['Name', 'NCBI ID']))
        ws.append(Row([taxon, taxonomy_util.Taxon(taxon).get_ncbi_id()]))

        # kinetics
        style['Kinetics'] = WorksheetStyle(head_rows=1, head_columns=0,
                                           head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)

        ws = wb['Kinetics'] = Worksheet()
        ws.append(Row([
            'Name', 'Sabio Reaction IDs',
            'Vmax Median Value', 'Vmax Median Entry', 'Vmax Min Entry', 'Vmax Max Entry', 'Vmax Proximity', 'Vmax Lift Info' 'Closest Vmax', 'Sabio Entry IDs', 'Closest Vmax Values',
            'Km Median Value', 'Km Median Entry', 'Km Min Enry', 'Km Max Entry', 'Km Proximity', 'Km Lift Info', 'Closest Km Sabio Entry IDs', 'Closest Km Values'
        ]))
        for i_row, rxn in enumerate(reactions):
            ws.append(Row([
                rxn.id,
                cls._format_list(rxn.reaction_ids),

                # vmax information
                cls._format_object(rxn.vmax_data.median_entry, "vmax"),
                cls._format_object(rxn.vmax_data.median_entry),
                cls._format_object(rxn.vmax_data.min_entry),
                cls._format_object(rxn.vmax_data.max_entry),
                cls._format_proximity(rxn.vmax_data.closest_entries),
                cls._format_object(rxn.vmax_data, "lift_info"),
                cls._format_list(rxn.vmax_data.closest_entry_ids),
                cls._format_list(rxn.vmax_data.closest_values),

                # km information
                cls._format_object(rxn.km_data.median_entry, "km"),
                cls._format_object(rxn.km_data.median_entry),
                cls._format_object(rxn.km_data.min_entry),
                cls._format_object(rxn.km_data.max_entry),
                cls._format_proximity(rxn.km_data.closest_entries),
                cls._format_object(rxn.km_data, "lift_info"),
                cls._format_list(rxn.km_data.closest_entry_ids),
                cls._format_list(rxn.km_data.closest_values),
            ]))

        # save to file
        write(filename, wb, style=style)

    @classmethod
    def run2(cls, taxon, reactions, filename):
        """ Save a list of reaction kinetic data to an Excel workbook or set of csv/tsv files

        Args:
            taxon (:obj:`str`): taxon name
            reactions (:obj:`list` of :obj:`kinetic_datanator.datanator.FormattedData`): list of reactions and their kinetic data
            filename (:obj:`str`): filename to store the list reaction kinetic data
        """
        # kinetics
        style['Kinetics'] = WorksheetStyle(head_rows=2, head_columns=0,
                                           head_row_font_bold=True, head_row_fill_fgcolor='CCCCCC', row_height=15)

        ws = wb['Kinetics'] = Worksheet()
        ws.append(Row([
            '',
            'Vmax', '', '', '', '', '' '', '', '',
            'Km', '', '', '', '', '' '', '', '',
            'Cross references', '', '',
        ]))
        ws.append(Row([
            'ID',
            'Median value', 'Median entry', 'Proximity', 'Lift Info' 'Closest Vmax', 'Min', 'Max', 'Sabio Entry IDs', 'Closest Values',
            'Median', 'Median Entry', 'Proximity', 'Lift Info' 'Closest Vmax', 'Min', 'Max', 'Sabio Entry IDs', 'Closest Values',
            'EC number', 'Predicted EC number', 'SABIO IDs',
        ]))

        for i_row, rxn in enumerate(reactions):
            row = [
                rxn.id,

                # vmax
                cls._format_object(rxn.vmax_data.median_entry, "vmax"),
                cls._format_object(rxn.vmax_data.median_entry),
                cls._format_object(rxn.vmax_data.min_entry),
                cls._format_object(rxn.vmax_data.max_entry),
                cls._format_proximity(rxn.vmax_data.closest_entries),
                cls._format_object(rxn.vmax_data, "lift_info"),
                cls._format_list(rxn.vmax_data.closest_entry_ids),
                cls._format_list(rxn.vmax_data.closest_values),

                # km
                cls._format_object(rxn.km_data.median_entry, "km"),
                cls._format_object(rxn.km_data.median_entry),
                cls._format_object(rxn.km_data.min_entry),
                cls._format_object(rxn.km_data.max_entry),
                cls._format_proximity(rxn.km_data.closest_entries),
                cls._format_object(rxn.km_data, "lift_info"),
                cls._format_list(rxn.km_data.closest_entry_ids),
                cls._format_list(rxn.km_data.closest_values),

                # cross references
                rxn.ec_number,
                cls._format_list(rxn.predicted_ec_numbers),
                cls._format_list(rxn.reaction_ids),
            ]

            ws.append(Row(row))

    @classmethod
    def _format_object(cls, obj, attr=None):
        """ Generate a string representation of an object or of an attribute of an object

        Args:
            obj (:obj:`object`): object
            attr (:obj:`str`, optional): attribute to format

        Returns:
            :obj:`str`: string representation of the object or of the attribute of the object
        """
        if obj is not None:
            if attr:
                return str(getattr(obj, attr))
            else:
                return str(obj.__dict__)

        return ""

    @classmethod
    def _format_proximity(cls, entries):
        """ Generate a string representation of the proximity of the first of a list of entry

        Args:
            entries (:obj:`list of :obj:`kinetic_datanator.sabio_interface.Entry`): list of entries

        Returns:
            :obj:`str`: string representation of the proximity of the first of a list of entry
        """
        if entries:
            return str(entries[0].proximity)
        else:
            return ""

    @classmethod
    def _format_list(cls, values):
        """ Generate a comma-separated string representation of a list of values

        Args:
            values (:obj:`list`): list of values

        Returns:
            :obj:`str`: comma-separated string representation of a list of values
        """
        return ', '.join(str(val) for val in values)
