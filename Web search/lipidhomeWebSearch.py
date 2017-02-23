"""This module takes the input dataframe and websearch parameters and searches
the LipidHome database.
"""
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
import urllib.request
import urllib.error
import urllib.parse
import urllib.request
import urllib.parse
import urllib.error
import pandas as pd
import time
import json


def liphomeSearch(inputFile, lipidhomeTol):
    """Main function that searches the LipidHome database. 
    2 output files are generated: one has the search results, the other
    lists which input data were not found in the database.

    Args:
        inputFile (dataframe): input data file 
        lipidhomeTol (float): m/z tolerance 

    Returns:
        int: flag is returned indicating whether the LipidHome search
        completed successfully (=1) or not (=0)
    """
    print("********************")
    failureTag = 0

    numberMasses = len(inputFile)
    mass = inputFile['MZ']
    retentionTime = inputFile['Time']
    polarity = inputFile['Polarity']
    columnType = inputFile['Column-type']

    # set up empty dataframes
    liphomeOutputFile = pd.DataFrame(index=[0])
    liphomeResponseDataframe = pd.DataFrame(index=liphomeOutputFile.index)
    liphomeNoOutputFile = pd.DataFrame(index=[0])
    liphomeNoResultOutputDF = pd.DataFrame(index=[0])

    try:
        # Lipidhome website
        liphomeURL = "http://www.ebi.ac.uk/metabolights/lipidhome/service/tools/ms1search"

        for count in range(numberMasses):
            # query arguments in dictionary format
            if polarity[count] in ['POS', 'POSC']:
                liphomeQueryArgs = {'masses': mass[
                    count], 'level': 'specie', 'tolerance': lipidhomeTol, 'identified': '', 'adductIons': '[2,3,4,5]'}
            elif polarity[count] in ['NEG', 'NEGC']:
                liphomeQueryArgs = {'masses': mass[
                    count], 'level': 'specie', 'tolerance': lipidhomeTol, 'identified': '', 'adductIons': '[6,7,8]'}

            print("Searching LipidHome for mass " + str(mass[count]))

        # using POST (i.e. include data as well as url), rather than GET:
            liphomeData = urllib.parse.urlencode(liphomeQueryArgs)
            binary_data = liphomeData.encode('utf8')
            liphomeReq = urllib.request.Request(liphomeURL, binary_data)
            liphomeResponse = urllib.request.urlopen(liphomeReq, timeout=10)
            str_response = liphomeResponse.read().decode('utf-8')
            liphomeResults = json.loads(str_response)

            countResults = len(liphomeResults['list'])

        # No results found, so print mass, etc. to no results output file
            if countResults == 0:
                liphomeNoResultOutputDF['Mass'] = mass[count]
                liphomeNoResultOutputDF['Time'] = retentionTime[count]
                liphomeNoResultOutputDF['Polarity'] = polarity[count]
                liphomeNoResultOutputDF['Column-type'] = columnType[count]
                liphomeNoOutputFile = liphomeNoOutputFile.append(
                    liphomeNoResultOutputDF, ignore_index=True)

            else:
                for index in range(0, countResults):
                    liphomeResponseDataframe['DATABASE_ID'] = liphomeResults[
                        'list'][index]['itemId']
                    liphomeResponseDataframe['NAME'] = liphomeResults[
                        'list'][index]['name']
                    liphomeResponseDataframe['mass'] = liphomeResults[
                        'list'][index]['mass']
                    liphomeResponseDataframe['identified'] = liphomeResults[
                        'list'][index]['identified']
                    liphomeResponseDataframe['faCarbons'] = liphomeResults[
                        'list'][index]['faCarbons']
                    liphomeResponseDataframe['faDoubleBonds'] = liphomeResults[
                        'list'][index]['faDoubleBonds']
                    liphomeResponseDataframe['RESULT_MASS'] = liphomeResults[
                        'list'][index]['resMass']
                    liphomeResponseDataframe['ADDUCT_ION'] = liphomeResults[
                        'list'][index]['adductIon']
                    liphomeResponseDataframe['CATEGORY'] = liphomeResults[
                        'list'][index]['code']
                    liphomeResponseDataframe['type'] = liphomeResults[
                        'list'][index]['type']
                    liphomeResponseDataframe['ORIGINAL_MASS'] = mass[count]
                    liphomeResponseDataframe[
                        'RETENTION_TIME'] = retentionTime[count]
                    liphomeResponseDataframe['POLARITY'] = polarity[count]
                    liphomeResponseDataframe['COLUMN_TYPE'] = columnType[count]
                    liphomeResponseDataframe['DELTA_PPM'] = (liphomeResponseDataframe.ORIGINAL_MASS.astype(float).fillna(
                        0.0) - liphomeResponseDataframe.RESULT_MASS.astype(float).fillna(0.0)) * (1000000) / (liphomeResponseDataframe.ORIGINAL_MASS.astype(float).fillna(0.0))

                # This part is using the itemId to search again (different url)
                # for the corresponding formula (not a POST request!)
                    liphomeURL2 = "http://www.ebi.ac.uk/metabolights/lipidhome/service/specie/summary"
                    liphomeQuery2 = {'id': liphomeResults[
                        'list'][index]['itemId']}
                    liphomeData2 = urllib.parse.urlencode(liphomeQuery2)
                    liphomeReq2 = liphomeURL2 + "?" + liphomeData2
                    liphomeResponse2 = urllib.request.urlopen(
                        liphomeReq2, timeout=10)
                    str_response2 = liphomeResponse2.read().decode('utf-8')
                    liphomeResults2 = json.loads(str_response2)
                    liphomeResponseDataframe[
                        'FORMULA'] = liphomeResults2['data']['formula']

                    liphomeOutputFile = liphomeOutputFile.append(
                        liphomeResponseDataframe, ignore_index=True)

    # HTTP error (must be BEFORE URL error), usually server error of some sort
    except urllib.error.HTTPError:
        print ("LipidHome connection failure")
        failureTag = 1

    # URL error, usually network connection failure
    except urllib.error.URLError:
        print("LipidHome connection failure")
        failureTag = 1

    # IO error, usually file open error or something along those lines
    except IOError:
        print("LipidHome connection failure")
        failureTag = 1

# regardless of whether search finished successfully, print out results
# (and 'no results') so far to output files
    # add date and time stap to output file name
    dt = time.strftime("%Y%m%d-%H%M%S")

    # only save the results output as a .csv if there are actual results!
    if not liphomeOutputFile.empty:
        liphomeOutputFile = liphomeOutputFile.dropna(axis=0, how='all')
        liphomeOutputFile = liphomeOutputFile.drop(
            ['faCarbons', 'faDoubleBonds', 'type', 'mass', 'identified'], 1)
        liphomeOutputFile.insert(0, 'DATABASE', 'LIPIDHOME')
        liphomeOutputFile.to_csv(
            "lipidhome_results_" + dt + ".csv", sep=',', float_format='%.6f', index=False)

    # only save the no results output as a .csv if there is data in the 'no
    # results'
    if not liphomeNoOutputFile.empty:
        liphomeNoOutputFile = liphomeNoOutputFile.dropna(axis=0, how='all')
        liphomeNoOutputFile.to_csv("lipidhome_no_results_" + dt + ".csv", sep=',', index=False, columns=[
                                   'Mass', 'Time', 'Polarity', 'Column-type'], float_format='%.6f', encoding="utf8")

    print("********************")

    return failureTag
