"""This is the main module for the online database search of LIPIDMAPS, LipidHome and HMDB.
Each database is searched separately and with different parameter settings, taken from the
websearch_parameters.csv file. 
The input file is a .csv format and is generated from SIEVE or XCMS and processed through 
Amalgamator (if combining positive and negative mode data). 
If a database search does not complete successfully, data up to that point is saved and 
the program moves to the next database search. At the end of the run, there is a listing on
screen to indicate which database searches (if any) did not complete successfully.

Attributes:
    endTime (float): Time the program ended running
    failTaghm (int): If equal to 1, the HMDB database search did not complete successfully
    failTaglh (int): If equal to 1, the LipidHome database search did not complete successfully
    failTaglm (int): If equal to 1, the LIPIDMAPS database search did not complete successfully
    inputFile (dataframe): The data from the input file (entered by the user)
    startTime (float): Time the program started running
"""
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
import time
import pandas as pd
import io

import getFile as gif
import lipidMapsWebSearch as lms
import lipidhomeWebSearch as lhs
import hmdbWebSearch as hs


# Disable chained assignment warning
pd.options.mode.chained_assignment = None

# gets input file as a DF
inputFile = gif.readFile()
lipidmapsTol, lipidhomeTol, hmdbTol, searchCriteria = gif.readParameters()

startTime = time.time()
# At the end each of these flags will tell us if a database search completed successfully (=1)
# Initialised to zero here
failTaglm = 0
failTaglh = 0
failTaghm = 0

# new option ot search by 'computational', 'curated' or both
# set this option in the websearch_parameters.csv

if searchCriteria == 'ALL':
    LMstatus = 'all'
    # comment out whichever website you DON'T want to search
    failTaglm = lms.lipmapSearch(inputFile, lipidmapsTol, LMstatus)
    failTaglh = lhs.liphomeSearch(inputFile, lipidhomeTol)
    failTaghm = hs.hmdbSearch(inputFile, hmdbTol)

if searchCriteria == 'COM':
    LMstatus = 'computational'
    # comment out whichever website you DON'T want to search
    failTaglm = lms.lipmapSearch(inputFile, lipidmapsTol, LMstatus)
    failTaglh = lhs.liphomeSearch(inputFile, lipidhomeTol)
    failTaghm = hs.hmdbSearch(inputFile, hmdbTol)

if searchCriteria == 'CUR':
    LMstatus = 'curated'
    failTaglm = lms.lipmapSearch(inputFile, lipidmapsTol, LMstatus)

endTime = time.time()

print("\n\nSUMMARY:")

print("Processing complete in", int((endTime - startTime)), "seconds\n")

if failTaglm != 0:
    print('LIPIDMAPS search did not complete successfully\n')
if failTaglh != 0:
    print('LipidHome search did not complete successfully\n')
if failTaghm != 0:
    print('HMDB search did not complete successfully\n')
