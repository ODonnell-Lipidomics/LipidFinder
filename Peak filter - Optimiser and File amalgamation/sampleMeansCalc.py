from builtins import range
import pandas as pd
import numpy as np
import re


def processSampleMeans(lfd):
    """Calculates and inserts sample replicate means into 'processedDataFrame'

    Args:
        lfd (LipidFinderData): A LipidFinderData object that has been previously had all clean steps performed on sample replicates.

    """
    lfd.meanColList = []
    if lfd.numberOfTechReps > 1:
        for firstRep in range(lfd.firstRepOffset, lfd.lastRepOffset + lfd.numberOfSamples, lfd.numberOfTechReps + 1):
            meanName = __colName__(lfd.processedDataFrame, firstRep)

            # Get means (not taking into account zeros) of the solvent samples
            rawMeans = lfd.processedDataFrame.ix[:, firstRep:firstRep + lfd.numberOfTechReps].apply(
                lambda x: x[np.nonzero(x)[0]].mean(), axis=1).fillna(0)

            # Round to nearest integer and convert means to integer types
            means = np.rint(rawMeans).astype(int)

            # Insert solvent mean into processedDataFrame
            lfd.processedDataFrame.insert(
                firstRep + lfd.numberOfTechReps, meanName + '_mean', means)
            lfd.meanColList.append(firstRep + lfd.numberOfTechReps)
        # increase to allow for mean field
        lfd.numberOfTechReps += 1
    else:
        lfd.meanColList = [i for i in range(
            lfd.firstRepOffset, (lfd.firstRepOffset + lfd.numberOfSamples))]
        # Cannot reassign names to columns.value
        colMeanNames = lfd.processedDataFrame.columns.values
        colMeanNames[lfd.meanColList] = colMeanNames[lfd.meanColList] + '_mean'
        lfd.processedDataFrame.columns = colMeanNames

    # To maintain correct indexing
    lfd.setLastRepOffset()
    lfd.MeanedData = lfd.processedDataFrame.copy()


def __colName__(rawData, nameRepStartInd):
    nameRaw = rawData.iloc[
        :, nameRepStartInd:nameRepStartInd + 1].columns.tolist()[0]
    return re.sub('\d+$', "", nameRaw)
