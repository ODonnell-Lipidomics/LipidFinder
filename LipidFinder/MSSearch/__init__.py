# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Search the lipid-alike features on the selected LIPID MAPS database
for bulk structure identification.

The output XLSX file will include every feature with its matched lipid
bulk structure and its relevant information such as lipid category or
formula. It will also include the unmatched features for completeness.
Optionally, it will generate a lipid category scatter plot from the
putative identification. Additionally, a log file is created to keep
track of the number of m/z values in the input data that got at at least
one match on the selected database for the given parameters.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> import MSSearch
    >>> parameters = LFParameters('mssearch', 'parameters.json')
    >>> data = LFDataFrame('filtered_data.csv', parameters)
    >>> MSSearch.bulk_structure_search(data, parameters)

    It can also work with a pandas.DataFrame as input:
    >>> import pandas
    >>> from Configuration import LFParameters
    >>> import MSSearch
    >>> parameters = LFParameters('mssearch', 'parameters.json')
    >>> data = pandas.read_csv('filtered_data.csv')
    >>> MSSearch.bulk_structure_search(data, parameters)
"""

from collections import defaultdict
import json
import logging
import os
import re
import requests
import time
import warnings

import numpy
import pandas
import pkg_resources
from requests_toolbelt.multipart.encoder import MultipartEncoder

from LipidFinder.MSSearch import DataPlots
from LipidFinder.MSSearch import Summary
from LipidFinder._py3k import viewitems, StringIO, quote_plus
from LipidFinder._utils import print_progress_bar


# Ignore Future Warnings from pandas library
warnings.simplefilter(action='ignore', category=FutureWarning)
# Deactivate pandas warnings
pandas.options.mode.chained_assignment = None
LIPIDMAPS_URL = 'https://www.lipidmaps.org/tools/ms/py_bulk_search.php'
# Maximum number of m/z values to send at once to LIPID MAPS
BATCH_SIZE = 150


def bulk_structure_search(data, parameters, dst=''):
    # type: (object, LFParameters, str) -> None
    """Search in LIPID MAPS for matches of the m/z values in the input
    dataframe.

    'data' must have, at least, m/z, retention time (RT) and "Polarity"
    columns. The adducts included in the search as well as the specific
    in-house lipidomics database, the mass tolerance and the lipid
    categories are provided in 'parameters'.
    The resulting dataframe will include every bulk structure match for
    each m/z, including its RT, main class, category and other relevant
    information. If 'dst' is not an absolute path, the current working
    directory will be used as starting point. If "mssearch_<db>.xslx"
    already exists, it will be overwritten without warning.
    "<db>" stands for the selected LIPID MAPS database.

    Keyword arguments:
        data       -- LFDataFrame or pandas.DataFrame instance
        parameters -- LipidFinder's MS Search parameters instance
        dst        -- destination directory where the log file, the
                      output XSLX file and the category scatter plot
                      figure (if selected) will be saved
                      [default: current working directory]
    """
    # Set the log file where the information about the steps performed
    # is saved
    logFilePath = 'mssearch.log'
    if (dst):
        logFilePath = os.path.join(dst, logFilePath)
    # Create logger and its file handler
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(logFilePath)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(asctime)s] %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    # Write initial information in log file
    logger.info('Starting MS Search on %s. Input dataframe has %d rows.',
                 parameters['database'], len(data.index))
    # Start progress bar
    progress = 0
    print_progress_bar(progress, 100, prefix='MSSearch progress:')
    mzCol = parameters['mzCol']
    rtCol = parameters['rtCol']
    # Get the list of unique m/z values from 'data'
    mzList = data[mzCol].unique().tolist()
    numMZ = len(mzList)
    # Write initial information in log file
    logger.info('%d unique m/z values found.', numMZ)
    # Get the list of target adducts from the parameters
    targetAdducts = parameters['targetAdducts']
    if (not targetAdducts):
        # If the list is empty, use the complete list of ion adducts
        targetAdducts = parameters._parameters['targetAdducts']['options']
    # Keep only the adduct information between brackets
    targetAdducts = [x[x.find('[') + 1 : x.find(']')] for x in targetAdducts]
    targetAdducts = ','.join(targetAdducts)
    if (parameters['mzToleranceUnit'] == 'Daltons'):
        tolerance = parameters['mzTolerance']
    # Get matches in batches to balance the number of requests and the
    # amount of information requested
    matches = pandas.DataFrame()
    # Calculate progress increment for each batch
    increment = 63.0 / numpy.ceil(float(numMZ) / BATCH_SIZE)
    for start in range(0, numMZ, BATCH_SIZE):
        mzBatch = mzList[start : start + BATCH_SIZE]
        # Get a string with one m/z per line (text file alike)
        mzStr = os.linesep.join(map(str, mzBatch))
        if (parameters['mzToleranceUnit'] == 'PPM'):
            # Calculate maximum tolerance in Da from tolerance in parts
            # per million (ppm)
            tolerance = mzBatch[-1] * parameters['mzTolerance'] / 1e6
        # Create the data package with the query
        if (parameters['categories']):
            mpData = MultipartEncoder(
                    fields={'CHOICE': parameters['database'], 'sort': 'DELTA',
                            'file': ('file', StringIO(mzStr), 'text/plain'),
                            'tol': str(tolerance), 'ion': targetAdducts,
                            'even': '2' if parameters['evenChains'] else '1',
                            'category': ','.join(parameters['categories'])})
        else:
            mpData = MultipartEncoder(
                    fields={'CHOICE': parameters['database'], 'sort': 'DELTA',
                            'file': ('file', StringIO(mzStr), 'text/plain'),
                            'tol': str(tolerance), 'ion': targetAdducts,
                            'even': '2' if parameters['evenChains'] else '1'})
        # Request the table containing the matches from LIPID MAPS
        try:
            response = requests.post(
                    LIPIDMAPS_URL, data=mpData,
                    headers={'Content-Type': mpData.content_type})
        except:
            raise Exception(("Connection error with the database. Please check "
                             "your network and try again after a few minutes."))
        # Go to next batch if this one returned nothing
        if (len(response.text) == 0):
            # Update progress bar
            progress += increment
            print_progress_bar(progress, 100, prefix='MSSearch progress:')
            continue
        # Process the response to create a dataframe
        batchMatches = pandas.read_csv(StringIO(response.text), sep='\t',
                                         engine='python', index_col=False)
        if (batchMatches.empty):
            # Update progress bar
            progress += increment
            print_progress_bar(progress, 100, prefix='MSSearch progress:')
            continue
        # Join all the information already gathered
        matches = matches.append(batchMatches, ignore_index=True)
        # Update progress bar
        progress += increment
        print_progress_bar(progress, 100, prefix='MSSearch progress:')
    if (matches.empty):
        matches = pandas.DataFrame(
                columns=[mzCol, 'Matched MZ', 'Delta', 'Bulk Structure',
                         'Formula', 'Adduct', 'Main Class', 'Category'])
    else:
        # Rename m/z column
        matches.rename(columns={'Input Mass': mzCol}, inplace=True)
    # Round 'Input Mass' values that might have been altered by LIPID
    # MAPS server
    matches[mzCol] = matches[mzCol].apply(round, ndigits=data._resolution)
    # Calculate the delta PPM of each row and add it to the dataframe
    dPPM = abs(matches[mzCol] - matches['Matched MZ']) * 1e6 / matches[mzCol]
    matches.insert(2, 'Delta_PPM', dPPM)
    if (parameters['mzToleranceUnit'] == 'PPM'):
        # Make sure all the matches comply with the m/z tolerance in ppm
        matches = matches[matches['Delta_PPM'] <= parameters['mzTolerance']]
    # Add RT and polarity values to each existing record and include the
    # rows in 'data' that did not have a match
    matches.insert(3, rtCol, 0.0)
    matches.insert(4, 'Polarity', '')
    # Calculate progress increment for each batch
    increment = 33.0 / numpy.ceil(len(data) / float(BATCH_SIZE))
    # Create result dataframe with all the columns in that dataframe
    colNames = [x for x in list(data) if x not in [mzCol, rtCol, 'Polarity']]
    extraCols = []
    for column in colNames:
        if (column not in list(matches)):
            extraCols.append(column)
        else:
            # Keep all columns from source dataset, adding prefix "src_"
            # if that column name is already in the dataframe
            extraCols.append('src_' + column)
            data.rename(columns={column: 'src_' + column}, inplace=True)
    result = pandas.DataFrame(columns=list(matches) + extraCols)
    # Ensure the polarity column contains only strings so the
    # conditional test in the next loop works as expected
    data['Polarity'].replace(numpy.nan, '', regex=True, inplace=True)
    # For those m/z values with more than one RT, the whole set of
    # matches is replicated for every RT
    for index, row in data.iterrows():
        mzMatches = matches.loc[matches[mzCol] == row[mzCol]]
        # Remove positive adduct matches for m/z found in negative mode,
        # and negative adduct matches for m/z found in positive mode
        if (row['Polarity'].lower().startswith('n')):
            mzMatches = mzMatches.loc[mzMatches['Adduct'].str[-1] != '+']
        elif (row['Polarity'].lower().startswith('p')):
            mzMatches = mzMatches.loc[mzMatches['Adduct'].str[-1] != '-']
        if (mzMatches.empty):
            # Unmatched m/z from 'data'
            mzMatches = mzMatches.append(row[[mzCol, rtCol, 'Polarity']],
                                         ignore_index = True)
        else:
            # Copy RT and polarity values to each matched m/z
            mzMatches[rtCol] = row[rtCol]
            mzMatches['Polarity'] = row['Polarity']
        # Copy the extra columns (if any) to each matched m/z
        for col in extraCols:
            mzMatches[col] = row[col]
        result = result.append(mzMatches, ignore_index=True)
        if ((index + 1) % BATCH_SIZE == 0):
            # Update progress bar
            progress += increment
            print_progress_bar(progress, 100, prefix='MSSearch progress:')
    # Sort the results by m/z, delta PPM and matched m/z to ease the
    # manipulation of the output XLSX file
    result.sort_values([mzCol, 'Delta_PPM', 'Matched MZ'], inplace=True,
                       kind='mergesort')
    # Create the XLSX file with the whole putative profiling dataframe
    outPath = os.path.join(
            dst, 'mssearch_{0}.xlsx'.format(parameters['database'].lower()))
    result.to_excel(outPath, index=False, engine='xlsxwriter')
    if (parameters['summary']):
        # Create summary XLSX file from the putative profiling
        # dataframe, keeping only one row per m/z and RT with the most
        # frequent lipid category
        Summary.create_summary(result, parameters, dst)
    # Update progress bar
    print_progress_bar(98, 100, prefix='MSSearch progress:')
    # Generate the category scatter plot of the most common lipid
    # category per m/z and RT
    if (parameters['plotCategories']):
        DataPlots.category_scatterplot(result, parameters, dst)
    # Update progress bar
    print_progress_bar(100, 100, prefix='MSSearch progress:')
    # Write the final information in log file and close handler
    matches = result[result['Category'].notna()]
    logger.info('MS Search completed. %d matches found for %d m/z values.\n',
                 len(matches), len(matches[mzCol].unique()))
    handler.close()
    logger.removeHandler(handler)
