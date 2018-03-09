import pandas
import os


def __secondsToMinutes__(seconds):
    minutes = seconds / 60.0
    return float('{:.2f}'.format(minutes))


def reformat(fileName, firstRepOffset):
    """Imports and reformats the XCMS csv to be compatible with the LipidFinder process

    Args:
        fileName (str): The XCMS file name
        firstRepOffset (int): Index of the first intensity column

    Returns:
        pandas DataFrame: The reformatted XCMS file in a DataFrame
    """
    input_frame = pandas.read_csv(fileName)
    # Compare columns in lowercase
    colNames = [x.lower() for x in list(input_frame)]
    colsToKeep = [0]
    if ("mz" in colNames) :
        colsToKeep.append(colNames.index("mz"))
    else:
        colsToKeep.append(colNames.index("mzmed"))
    if ("rt" in colNames) :
        colsToKeep.append(colNames.index("rt"))
    else:
        colsToKeep.append(colNames.index("rtmed"))

    # dropping and renaming the columns based on their index rather than their
    # names, as the columns from the input file can have unexpected names.
    colsToDrop = [x for x in range(0, firstRepOffset) if x not in colsToKeep]
    input_frame = input_frame.drop(input_frame.columns[colsToDrop], axis=1)

    input_frame_columns = input_frame.columns.values
    input_frame_columns[0] = 'id'
    input_frame_columns[1] = 'MZ'
    input_frame_columns[2] = 'Time'
    input_frame.columns = input_frame_columns

    input_frame['Time'] = input_frame['Time'].apply(__secondsToMinutes__)

    return input_frame
