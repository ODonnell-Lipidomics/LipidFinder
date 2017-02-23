"""This module takes the input dataframe and websearch parameters and searches
the LIPIDMAPS database.
"""
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
import urllib.request
import urllib.error
import urllib.parse
import urllib.request
import urllib.parse
import urllib.error
import pandas as pd
import numpy as np
import time
import io


def lipmapSearch(inputFile, lipidmapsTol, LMstatus):
    """Main function that searches the LIPIDMAPS database. 
    2 output files are generated: one has the search results, the other
    lists which input data were not found in the database.

    Args:
        inputFile (dataframe): input data file 
        lipidmapsTol (float): m/z tolerance 
        LMstatus (str): database search option ("COM", "CUR" or "ALL")

    Returns:
        int: flag is returned indicating whether the LIPIDMAPS search
        completed successfully (=1) or not (=0)
    """
    print("********************")
    failureTag = 0

    delta = 1.007281
    inputFile.insert(0, 'lipidmaps_mass', np.nan)
    inputFile.insert(0, 'ADDUCT_ION', np.nan)

    # set up empty dataframes
    lipmapOutputFile = pd.DataFrame(index=[0])
    lipmapNoOutputFile = pd.DataFrame(index=[0])
    lipmapNoResultOutputDF = pd.DataFrame(index=[0])

    # convert POS and NEG mz values to NEUTRAL for LIPID MAPS search
    for i1, r1 in inputFile.iterrows():
        if r1['Polarity'] in ['POS', 'POSC']:
            inputFile.loc[i1, 'lipidmaps_mass'] = inputFile.loc[
                i1, 'MZ'] - delta
            inputFile.loc[i1, 'ADDUCT_ION'] = "[M+H]+"
        elif r1['Polarity'] in ['NEG', 'NEGC']:
            inputFile.loc[i1, 'lipidmaps_mass'] = inputFile.loc[
                i1, 'MZ'] + delta
            inputFile.loc[i1, 'ADDUCT_ION'] = "[M-H]-"

    mass = inputFile['MZ']
    retentionTime = inputFile['Time']
    polarity = inputFile['Polarity']
    columnType = inputFile['Column-type']
    adductIon = inputFile['ADDUCT_ION']
    lipmapMass = inputFile['lipidmaps_mass']

    try:
        # LIPID MAPS url
        lipmapURL = "http://www.lipidmaps.org/data/structure/LMSDSearch.php"

        for index, row in inputFile.iterrows():
            print("Searching LIPID MAPS for mass (converted to NEUTRAL) " +
                  str(lipmapMass[index]))

            lipmapQueryArgs = {'Mode': 'ProcessTextSearch', 'status': LMstatus, 'OutputMode': 'File',
                               'OutputType': 'TSV', 'ExactMass': lipmapMass[index], 'ExactMassOffSet': lipidmapsTol}
            lipmapData = urllib.parse.urlencode(lipmapQueryArgs)
            lipmapGetUrl = lipmapURL + "?" + lipmapData
            lipmapResponse = urllib.request.urlopen(lipmapGetUrl, timeout=10)
            lipmapHtmlOutput = lipmapResponse.read().decode('utf-8')
            lipmapResponseDataframe = pd.read_csv(
                io.StringIO(lipmapHtmlOutput), sep='\t')

            if lipmapResponseDataframe.empty:    # no database matches
                lipmapNoResultOutputDF['MZ'] = mass[index]
                lipmapNoResultOutputDF['Time'] = retentionTime[index]
                lipmapNoResultOutputDF['Polarity'] = polarity[index]
                lipmapNoResultOutputDF['Column-type'] = columnType[index]
                lipmapNoOutputFile = lipmapNoOutputFile.append(
                    lipmapNoResultOutputDF, ignore_index=True)

            else:
                lipmapResponseDataframe['ORIGINAL_MASS'] = mass[index]
                lipmapResponseDataframe[
                    'RETENTION_TIME'] = retentionTime[index]
                lipmapResponseDataframe['POLARITY'] = polarity[index]
                lipmapResponseDataframe['ADDUCT_ION'] = adductIon[index]
                lipmapResponseDataframe['COLUMN_TYPE'] = columnType[index]
                lipmapOutputFile = lipmapOutputFile.append(
                    lipmapResponseDataframe, ignore_index=True)

    # HTTP error (must be BEFORE URL error), usually server error of some sort
    except urllib.error.HTTPError:
        print ("LIPID MAPS connection failure")
        failureTag = 1

    # URL error, usually network connection failure
    except urllib.error.URLError:
        print("LIPID MAPS connection failure")
        failureTag = 1

    # IO error, usually file open error or something along those lines
    except IOError:
        print("LIPID MAPS connection failure")
        failureTag = 1

# regardless of whether search finished successfully, print out results
# (and 'no results') so far to output files
    # add date and time stap to output file name
    dt = time.strftime("%Y%m%d-%H%M%S")

    # first check that at least one result was found. If no results, then
    # don't get a 'results' file
    if not lipmapOutputFile.empty:
        lipmapOutputFile = lipmapOutputFile.dropna(axis=0, how='all')

        # have to convert the resulting masses back to positive or negative:
        # this for loop makes the processing very slow, see if can find a
        # quicker way of doing it
        for i2, r2 in lipmapOutputFile.iterrows():
            if r2['POLARITY'] in ['POS', 'POSC']:
                lipmapOutputFile.loc[i2, 'MASS'] = lipmapOutputFile.loc[
                    i2, 'MASS'] + delta
            elif r2['POLARITY'] in ['NEG', 'NEGC']:
                lipmapOutputFile.loc[i2, 'MASS'] = lipmapOutputFile.loc[
                    i2, 'MASS'] - delta

        # rename column headers, calculate DELTA_PPM, drop unwanted columns,
        # add new DATABASE column
        lipmapOutputFile.rename(columns={'LM_ID': 'DATABASE_ID', 'SYSTEMATIC_NAME': 'NAME', 'MASS': 'RESULT_MASS',
                                         'COMMON_NAME': 'LIPIDMAPS_COMMON_NAME', 'MAIN_CLASS': 'LIPIDMAPS_MAIN_CLASS'}, inplace=True)
        lipmapOutputFile['DELTA_PPM'] = (lipmapOutputFile.ORIGINAL_MASS.astype(float).fillna(
            0.0) - lipmapOutputFile.RESULT_MASS.astype(float).fillna(0.0)) * (1000000) / (lipmapOutputFile.ORIGINAL_MASS.astype(float).fillna(0.0))
        lipmapOutputFile = lipmapOutputFile.drop('SUB_CLASS', 1)
        lipmapOutputFile.insert(0, 'DATABASE', 'LIPIDMAPS')

        lipmapOutputFile.to_csv(
            "lipidmaps_results_" + dt + ".csv", sep=',', float_format='%.6f', index=False)

    # only save the no results output as a .csv if there is data in the 'no
    # results'
    if not lipmapNoOutputFile.empty:
        lipmapNoOutputFile = lipmapNoOutputFile.dropna(axis=0, how='all')
        lipmapNoOutputFile.to_csv("lipidmaps_no_results_" + dt + ".csv", sep=',', index=False, columns=[
                                  'MZ', 'Time', 'Polarity', 'Column-type'], float_format='%.6f', encoding="utf8")

    print("********************")

    return failureTag
