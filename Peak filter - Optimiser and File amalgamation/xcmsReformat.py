from __future__ import division
from past.utils import old_div
import pandas as pds
import os


def __secondsToMinutes__(seconds):
    minutes = old_div(seconds, 60)
    return round(minutes, 2)


def reformat(fileName):
    """Imports and reformats the XCMS csv to be compatible with the LipidFinder process

    Args:
        fileName (str): The XCMS file name

    Returns:
        pandas DataFrame: The reformatted XCMS file in a DataFrame
    """
    input_frame = pds.read_csv(fileName)

    # dropping and renaming the columns based on their index rather than their
    # names, as the columns from the input file can have unexpected names.
    columns_to_drop = [1, 2, 3, 4, 6, 7, 9, 10, 11, 12, 13]
    input_frame = input_frame.drop(
        input_frame.columns[columns_to_drop], axis=1)

    input_frame_columns = input_frame.columns.values
    input_frame_columns[0] = 'id'
    input_frame_columns[1] = 'MZ'
    input_frame_columns[2] = 'Time'
    input_frame.columns = input_frame_columns

    input_frame['Time'] = input_frame['Time'].apply(__secondsToMinutes__)

    return input_frame
