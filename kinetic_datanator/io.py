""" Utilities for reading and writing data to/from files

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-14
:Copyright: 2017, Karr Lab
:License: MIT
"""

from .util import taxonomy_util
from wc_utils.workbook.core import Row, Workbook, Worksheet
from wc_utils.workbook.io import WorkbookStyle, WorksheetStyle, write


class InputReader(object):
    pass


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
