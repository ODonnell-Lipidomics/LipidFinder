import pandas as pd
import numpy as np
import re

import outlierCorrect as oc


def performSolCalcs(lfd, saveFrame=True):
    """Removes the effect of extracted blanks (solvent samples) on the the sample
    replictates

    Args:
        lfd (LipidFinderData): A LipidFinderData object
        saveFrame (bool, optional): Should the 'processedDataFrame' be archived following solvents calculations.

    """
    # Column index of first solvent replicate
    indFirst = lfd.getFirstSolRepOffset()

    # Column index of last solvent replicate
    indLast = lfd.getLastSolRepOffset()

    # Outlier correct solvent
    lfd.processedDataFrame = oc.process(
        lfd, indFirst, indLast, lfd.numberOfSolventReps)

    # Insert mean column for solvent into dataFrame
    solvMeanName = getSolvMeanName(lfd)

    solvMeans = __calcSolvMean__(indFirst, indLast, lfd.processedDataFrame)

    # Insert solvent mean into processedDataFrame
    lfd.processedDataFrame.insert(indLast, solvMeanName, solvMeans)

    lfd.numberOfSolventReps += 1

    removeMeanSolv(lfd, solvMeanName)

    if saveFrame:
        lfd.solventData = lfd.processedDataFrame.copy()


def getSolvMeanName(lfd):
    """Assigns a meaningful name for the solvent mean based on the name of the first solvent sample

    Args:
        lfd (LipidFinderData): A LipidFinderData object

    Returns:
        str: The name of the solvent mean field
    """
    return __colName__(lfd.processedDataFrame, lfd.getFirstSolRepOffset()) + '_mean'


def __calcSolvMean__(indFirst, indLast, processedDataFrame):

    # Get means (not taking into account zeros) of the solvent samples
    rawSolvMeans = processedDataFrame.ix[:, indFirst:indLast].apply(
        lambda x: x[np.nonzero(x)[0]].mean(), axis=1).fillna(0)

    # Round to nearest integer and convert means to integer types
    return np.rint(rawSolvMeans).astype(int)


def removeMeanSolv(lfd, solvMeanName):
    """Subtracts the solvent mean intensity value from each sample replicate

    Args:
        lfd (LipidFinderData): A LipidFinderData object
        solvMeanName (str): The name of the solvent mean field

    """
    # Solvent fold elimination, remove frames where all tech reps of all
    # samples are <solventFoldCutOff fold solvent mean after after outlier
    # correction
    lfd.processedDataFrame.drop(np.where(lfd.processedDataFrame.ix[:, lfd.firstRepOffset:lfd.lastRepOffset].max(
        axis=1) < lfd.solventFoldCutOff * lfd.processedDataFrame[solvMeanName])[0], inplace=True)

    # Remove solvent from all remaining samples, lowest value intensity == 0
    lfd.processedDataFrame.ix[:, lfd.firstRepOffset:lfd.lastRepOffset] = np.maximum(0, lfd.processedDataFrame.ix[
                                                                                    :, lfd.firstRepOffset:lfd.lastRepOffset].sub(lfd.processedDataFrame[solvMeanName], axis=0))

    # Reset the index of processDataFrame after the removals
    lfd.processedDataFrame.reset_index(inplace=True, drop=True)


def removeLowIntesityFrames(lfd, saveFrame=True):
    """Discard DataFrame rows (ions) where intensity is too low

    Args:
        lfd (LipidFinderData): A LipidFinderData object
        saveFrame (bool, optional): Should the 'processedDataFrame' be archived following solvents calculations.

    """
    # Set all replicate intensities that are <intensitySignificanceCutOff to 0
    temp = lfd.processedDataFrame.ix[:, lfd.firstRepOffset:lfd.lastRepOffset]
    temp[temp < lfd.intensitySignificanceCutOff] = 0
    lfd.processedDataFrame.ix[:, lfd.firstRepOffset:lfd.lastRepOffset] = temp

    # Intensity cut off - Remove frames where all tech reps of all samples are
    # <intensitySignificanceCutOff
    lfd.processedDataFrame.drop(np.where(lfd.processedDataFrame.ix[:, lfd.firstRepOffset:lfd.lastRepOffset].max(
        axis=1) < lfd.intensitySignificanceCutOff)[0], inplace=True)

    # Reset the index of processDataFrame after the removals
    lfd.processedDataFrame.reset_index(inplace=True, drop=True)

    if saveFrame:
        lfd.lowIntData = lfd.processedDataFrame.copy()


def __colName__(rawData, nameRepStartInd):
    nameRaw = rawData.iloc[
        :, nameRepStartInd:nameRepStartInd + 1].columns.tolist()[0]
    return re.sub('\d+$', "", nameRaw)
