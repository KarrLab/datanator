""" Utilities for reading and writing data to/from files

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-14
:Copyright: 2017, Karr Lab
:License: MIT
"""

from .core import data_structs
from .util import compartment_util
from .util import molecule_util
from .util import reaction_util
from .util import taxonomy_util
from wc_utils.workbook.core import Row, Workbook, Worksheet
from wc_utils.workbook.io import WorkbookStyle, WorksheetStyle, write
import re
import openpyxl


class InputReader(object):

    @classmethod
    def run(cls, filename):
        """ Read input data from an Excel workbook

        Args:
            filename (:obj:`str`): filename of Excel workbook

        Returns:
            :obj:`tuple`: 

                * :obj:`taxonomy_util.Taxon`: taxon
                * :obj:`list` of :obj:`compartment_util.Compartment`: list of compartments
                * :obj:`list` of :obj:`molecule_util.Molecule`: list of molecules
                * :obj:`list` of :obj:`reaction_util.Reaction`: list of reactions
        """
        wb = openpyxl.load_workbook(filename=filename)
        taxon = cls._read_taxon(wb.get_sheet_by_name('Taxon'))
        compartments = []
        molecules = cls._read_molecules(wb.get_sheet_by_name('Molecules'))
        reactions = cls._read_reactions(wb.get_sheet_by_name('Reactions'))
        
        cls.link_objects(compartments, molecules, reactions)

        return (taxon, compartments, molecules, reactions)

    @classmethod
    def _read_taxon(cls, ws):
        """ Read taxon from an Excel worksheet

        Args:
            ws (:obj:`openpyxl.Worksheet`): worksheet

        Returns:
            :obj:`taxonomy_util.Taxon`: taxon
        """
        return taxonomy_util.Taxon(name=ws.cell(row=2, column=1).value)

    @classmethod
    def _read_molecules(cls, ws):
        """ Read molecules from an Excel worksheet

        Args:
            ws (:obj:`openpyxl.Worksheet`): worksheet

        Returns:
            :obj:`list` of `molecule_util.Molecule`: list of molecules
        """
        molecules = []
        for i in range(2, ws.max_row + 1):
            molecules.append(molecule_util.Molecule(
                id=ws.cell(row=i, column=1).value,
                structure=ws.cell(row=i, column=2).value,
            ))

        return molecules

    @classmethod
    def _read_reactions(cls, ws):
        """ Read reactions from an Excel worksheet

        Args:
            ws (:obj:`openpyxl.Worksheet`): worksheet

        Returns:
            :obj:`list` of `reaction_util.Reaction`: list of reactions
        """
        reactions = []
        for i in range(2, ws.max_row + 1):
            rxn = cls.parse_reaction_equation(ws.cell(row=i, column=2).value)
            rxn.id = ws.cell(row=i, column=1).value
            reactions.append(rxn)
        return reactions

    @classmethod
    def parse_reaction_equation(cls, equation):
        """ Parse a reaction equation, e.g.

        * [c]: ATP + H2O ==> ADP + PI + H
        * GLC[e] + ATP[c] + H2O[c] ==> GLC[c] + ADP[c] + PI[c] + H[c]

        Args:
            equation (:obj:`str`): reaction equation

        Returns:
            :obj:`reaction_util.Reaction': reaction
        """
        global_comp = '\[(?P<comp>[a-z0-9_]+)\]: *'
        global_part = ' *(([0-9\.]+) )?([a-z0-9_]+)'
        global_lhs = '(?P<lhs>{}( *\+ *{})*)'.format(global_part, global_part)
        global_rhs = '(?P<rhs>{}( *\+ *{})*)'.format(global_part, global_part)
        local_part = ' *(([0-9\.]+) )?([a-z0-9_]+)\[([a-z0-9_]+)\]'
        local_lhs = '(?P<lhs>{}( *\+ *{})*)'.format(local_part, local_part)
        local_rhs = '(?P<rhs>{}( *\+ *{})*)'.format(local_part, local_part)
        sep = ' *(?P<sep><{0,1})[-=]{1,2}> *'

        global_pattern = global_comp + global_lhs + sep + global_rhs
        local_pattern = local_lhs + sep + local_rhs

        global_match = re.match(global_pattern, equation.rstrip(), re.IGNORECASE)
        local_match = re.match(local_pattern, equation.rstrip(), re.IGNORECASE)

        if global_match:
            global_comp = global_match.groupdict()['comp']
            sep = global_match.groupdict()['sep']
            lhs = global_match.groupdict()['lhs']
            rhs = global_match.groupdict()['rhs']

            participants = []
            for part in re.findall(global_part, lhs, re.IGNORECASE):
                participants.append(reaction_util.ReactionParticipant(
                    molecule=part[2],
                    compartment=global_comp,
                    coefficient=-float(part[0][0:-1] or 1.),
                ))
            for part in re.findall(global_part, rhs, re.IGNORECASE):
                participants.append(reaction_util.ReactionParticipant(
                    molecule=part[2],
                    compartment=global_comp,
                    coefficient=float(part[0][0:-1] or 1.),
                ))

        elif local_match:
            sep = local_match.groupdict()['sep']
            lhs = local_match.groupdict()['lhs']
            rhs = local_match.groupdict()['rhs']

            participants = []
            for part in re.findall(local_part, lhs, re.IGNORECASE):
                participants.append(reaction_util.ReactionParticipant(
                    molecule=part[2],
                    compartment=part[3],
                    coefficient=-float(part[0][0:-1] or 1.),
                ))
            for part in re.findall(local_part, rhs, re.IGNORECASE):
                participants.append(reaction_util.ReactionParticipant(
                    molecule=part[2],
                    compartment=part[3],
                    coefficient=float(part[0][0:-1] or 1.),
                ))
        else:
            raise ValueError('Reaction is not parseable: {}'.format(equation))

        # save specified participant order
        for i_part, part in enumerate(participants):
            part.order = i_part

        return reaction_util.Reaction(participants=participants, reversible=sep == '<')

    @classmethod
    def link_objects(cls, compartments, molecules, reactions):
        """ Link objects to build object graph

        Args:
            compartments (:obj:`list` of :obj:`compartment_util.Compartment`): list of compartments
            molecules (:obj:`list` of :obj:`molecule_util.Molecule`): list of molecules
            reactions (:obj:`list` of :obj:`reaction_util.Reaction`): list of reactions
        """
        for rxn in reactions:
            for part in rxn.participants:
                m_str = part.molecule
                m_obj = next((m for m in molecules if m.id == m_str), None)
                if not m_obj:
                    raise ValueError('Participant {} of reaction {} is not defined'.format(m_str, rxn.id))
                part.molecule = m_obj

                c_str = part.compartment
                c_obj = next((c for c in compartments if c.id == c_str), None)
                if not c_obj:
                    c_obj = compartment_util.Compartment(id=c_str)
                    compartments.append(c_obj)
                part.compartment = c_obj


class ResultsWriter(object):
    """ Save reaction kinetic data to an Excel workbook """

    @classmethod
    def run(cls, taxon, reactions, filename):
        """ Save a list of reaction kinetic data to an Excel workbook or set of csv/tsv files

        Args:
            taxon (:obj:`str`): taxon
            reactions (:obj:`list` of :obj:`kinetic_datanator.datanator.SabioResult`): list of reactions and their kinetic data
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
            reactions (:obj:`list` of :obj:`kinetic_datanator.datanator.SabioResult`): list of reactions and their kinetic data
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
            entries (:obj:`list of :obj:`kinetic_datanator.sabio_rk.Entry`): list of entries

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
