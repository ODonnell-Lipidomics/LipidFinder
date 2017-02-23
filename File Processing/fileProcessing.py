"""Main module that is used for processing and filtering the results
from the websearch darabase searches

Attributes:
    combinedFiles (dataframe): original merged (but unfiltered/processed) database results files
    deltaPPM (float): delta ppm cutoff value specified by user in websearch_parameters.csv
    endTime (float): Time the program ended running
    fileCount (int): Number of database files that need to be processed
    finalDF (dataframe): dataframe of filtered/processed database results files
    masscatsDF (dataframe): dataframe showing the most common lipid category per mass
    mergeDF (dataframe): interim dataframe in processing
    startTime (float): Time the program started running
"""
from __future__ import print_function
import pandas as pd
import time

import processingSteps as ps
import selectData as sd

# Disable chained assignment warning
pd.options.mode.chained_assignment = None


fileCount = ps.numFiles()
combinedFiles = ps.inputFiles(fileCount)
originalDF, originalFile = ps.inputOriginalfile()
deltaPPM = ps.readParameters()

startTime = time.time()

mergeDF = sd.deltaElim(combinedFiles, deltaPPM)
mergeDF = sd.adductRename(mergeDF)
mergeDF = sd.adductElim(mergeDF)
mergeDF = sd.categoryRename(mergeDF)
mergeDF, sortedFile = sd.removeDuplicateNames(mergeDF)
finalDF = mergeDF
masscatsDF = sd.getMasscategorycounts(finalDF, originalDF)

massdbDF1, massesDF2 = sd.getDatabasecategorycounts(originalDF, finalDF)

endTime = time.time()
print("\n\n    *** PROCESSING COMPLETE in",
      int((endTime - startTime)), "s ***\n")
