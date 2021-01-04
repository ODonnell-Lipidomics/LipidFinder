# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Amalgamate negative and positive ion polarity datasets.

Both datasets have to match the same column layout as the output files
from LipidFinder's PeakFilter module.

Every unique feature is kept in the output file. For those features with
matching m/z and retention time (within a tolerance) in both polarities,
the one with the lowest total intensity mean is discarded. This default
behavior can be changed to the combination (sum) of the intensities of
both features. Additionally, a log file is created to keep track of the
matches found.

Examples:
    >>> import pandas
    >>> from Configuration import LFParameters
    >>> import Amalgamator
    >>> negData = pandas.read_csv('peakfilter_negative.csv')
    >>> posData = pandas.read_csv('peakfilter_positive.csv')
    >>> parameters = LFParameters('amalgamator', 'parameters.json')
    >>> Amalgamator.amalgamate_data(negData, posData, parameters)
"""

import logging
from math import sqrt
import os
import warnings

import numpy
import pandas

from LipidFinder.LFDataFrame import LFDataFrame
from LipidFinder._utils import mz_delta, mz_tol_range, rt_delta, rt_tol_range, \
                               print_progress_bar


# Ignore Future Warnings from pandas library
warnings.simplefilter(action='ignore', category=FutureWarning)
# Calculate the mass differences produced by the protonation:
# Proton (H) mass
PROTON = 1.00727646
# Hydrogen (H2) mass
HYDROGEN = 2 * PROTON
# Electron mass = 0.00054858 (source:
# http://fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator)
# Carbon 12 mass = 12
# [CH3]- mass = 1 * C12 + 3 * H + 2 * electron
METHYL_ANION = 15.02292654
# CH4 mass = H + [CH3]-
METHANE = PROTON + METHYL_ANION


def amalgamate_data(negData, posData, parameters, dst=''):
    # type: (object, object, LFParameters, str) -> None
    """Amalgamate negative and positive ion polarity dataframes.

    'negData' and 'posData' have to match the same column layout as the
    output files from LipidFinder's PeakFilter module.
    For those frames with matching m/z and retention time, the one with
    the lowest total intensity mean is discarded. Both files must have
    the same column headings. If 'dst' is not an absolute path, the
    current working directory will be used as starting point. If
    "amalgamated.csv" file already exists, it will be overwritten.

    Keyword Arguments:
        negData    -- negative polarity LFDataFrame or pandas.DataFrame
                      instance
        posData    -- positive polarity LFDataFrame or pandas.DataFrame
                      instance
        parameters -- LipidFinder's Amalgamator parameters instance
        dst        -- destination directory where the log file and the
                      amalgamated data CSV file will be saved
                      [default: current working directory]
    """
    # Set the log file where the information about the steps performed
    # is saved
    logFilePath = 'amalgamator.log'
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
    logger.info(("Starting Amalgamator. Negative dataframe has %d rows and "
                  "Positive dataframe has %d rows."), len(negData.index),
                 len(posData.index))
    mzCol = parameters['mzCol']
    rtCol = parameters['rtCol']
    # Check if columns in both dataframes are the same
    if (set(negData.columns) != set(posData.columns)):
        diffCols = set(negData.columns).symmetric_difference(posData.columns)
        raise IOError(("Input dataframes do not share the same column names: "
                       "{0}").format(', '.join(diffCols)))
    # Check for misspelling errors in m/z or retention time column names
    if ((mzCol not in negData.columns) or (rtCol not in negData.columns)):
        raise KeyError("Missing '{0}' or '{1}' column(s)".format(mzCol, rtCol))
    # Get the indices for intensity columns
    firstIndex = parameters['firstSampleIndex'] - 1
    lastIndex = firstIndex + parameters['numSamples']
    # Calculate the mean of every non-zero value of the mean columns of
    # each input dataframe and round it to the nearest integer. Replace
    # any NaN output from mean() by zero.
    totalMean = lambda x: numpy.rint(numpy.nan_to_num(
            x[numpy.where(x>0)[0]].mean())).astype(int)
    negData['TotalMean'] = negData.iloc[:, firstIndex : lastIndex].apply(
            totalMean, axis=1)
    posData['TotalMean'] = posData.iloc[:, firstIndex : lastIndex].apply(
            totalMean, axis=1)
    nind = negData.index.values
    nmz = negData[mzCol].values
    nrt = negData[rtCol].values
    nmeans = negData['TotalMean'].values
    negCol = list(negData.columns.values)
    posCol = list(posData.columns.values)
    # Empty results dataframe
    results = pandas.DataFrame(columns=negCol)
    polColIndex = results.columns.get_loc('Polarity')
    # Start progress bar
    progress = 0
    total = len(nind) + 1
    print_progress_bar(progress, total, prefix='Amalgamator progress:')
    # Loop through indices in negative file
    for i in nind:
        # Update progress bar
        progress += 1
        print_progress_bar(progress, total, prefix='Amalgamator progress:')
        negMass = nmz[i]
        negRT = nrt[i]
        pmz = posData[mzCol].values
        prt = posData[rtCol].values
        pmeans = posData['TotalMean'].values
        negMassH2 = negMass + HYDROGEN
        mzRange = mz_tol_range(negMassH2, parameters['mzFixedError'],
                               parameters['mzPPMError'])
        rtRange = rt_tol_range(negRT, parameters['maxRTDiffAdjFrame'])
        matchesH2 = list(numpy.where(
                (pmz >= mzRange[0]) & (pmz <= mzRange[1]) & (prt >= rtRange[0])
                & (prt <= rtRange[1]))[0])
        # First, look for H2 matches
        if (matchesH2):
            indMatch = __bestMatch__(matchesH2, negMassH2, pmz, negRT, prt,
                                     parameters)
            # Keep the frame with the highest total mean
            if (pmeans[indMatch] > nmeans[i]):
                results = results.append(posData.iloc[indMatch],
                                         ignore_index=True)
                if (parameters['combineIntensities']):
                    results.iloc[-1, firstIndex : lastIndex] = \
                            results.iloc[-1, firstIndex : lastIndex] \
                            + negData.iloc[i, firstIndex : lastIndex]
                    results.iloc[-1, polColIndex] += ' (Combined)'
                else:
                    results.iloc[-1, polColIndex] += ' (Both)'
            else:
                results = results.append(negData.iloc[i], ignore_index=True)
                if (parameters['combineIntensities']):
                    results.iloc[-1, firstIndex : lastIndex] = \
                            results.iloc[-1, firstIndex : lastIndex] \
                            + posData.iloc[indMatch, firstIndex : lastIndex]
                    results.iloc[-1, polColIndex] += ' (Combined)'
                else:
                    results.iloc[-1, polColIndex] += ' (Both)'
            logger.info('Match found: Negative ID %d - Positive ID %d.',
                         negData.iloc[i, 0], posData.iloc[indMatch, 0])
            # Remove match from positive dataframe, avoiding writing
            # the action to the log file
            if (isinstance(posData, LFDataFrame)):
                super(LFDataFrame, posData).drop(indMatch, inplace=True)
            else:
                posData.drop(indMatch, inplace=True)
            posData.reset_index(inplace=True, drop=True)
            pmz = posData[mzCol].values
            prt = posData[rtCol].values
            pmeans = posData['TotalMean'].values
            continue
        # If there are no H2 matches, look for CH4 matches
        negMassCH4 = negMass + METHANE
        mzRange = mz_tol_range(negMassCH4, parameters['mzFixedError'],
                               parameters['mzPPMError'])
        matchesHCH3 = list(numpy.where(
                (pmz >= mzRange[0]) & (pmz <= mzRange[1]) & (prt >= rtRange[0])
                & (prt <= rtRange[1]))[0])
        if (matchesHCH3):
            indMatch = __bestMatch__(matchesHCH3, negMassCH4, pmz, negRT, prt,
                                     parameters)
            # Keep the frame with the highest total mean
            if (pmeans[indMatch] > nmeans[i]):
                results = results.append(posData.iloc[indMatch],
                                         ignore_index=True)
                if (parameters['combineIntensities']):
                    results.iloc[-1, firstIndex : lastIndex] = \
                            results.iloc[-1, firstIndex : lastIndex] \
                            + negData.iloc[i, firstIndex : lastIndex]
                    results.iloc[-1, polColIndex] += ' (Combined)'
                else:
                    results.iloc[-1, polColIndex] += ' (Both)'
            else:
                results = results.append(negData.iloc[i], ignore_index=True)
                if (parameters['combineIntensities']):
                    results.iloc[-1, firstIndex : lastIndex] = \
                            results.iloc[-1, firstIndex : lastIndex] \
                            + posData.iloc[indMatch, firstIndex : lastIndex]
                    results.iloc[-1, polColIndex] += ' (Combined)'
                else:
                    results.iloc[-1, polColIndex] += ' (Both)'
            logger.info('Match found: Negative ID %d - Positive ID %d.',
                         negData.iloc[i, 0], posData.iloc[indMatch, 0])
            # Remove match from positive dataframe, avoiding writing
            # the action to the log file
            if (isinstance(posData, LFDataFrame)):
                super(LFDataFrame, posData).drop(indMatch, inplace=True)
            else:
                posData.drop(indMatch, inplace=True)
            posData.reset_index(inplace=True, drop=True)
            pmz = posData[mzCol].values
            prt = posData[rtCol].values
            pmeans = posData['TotalMean'].values
            continue
        results = results.append(negData.iloc[i], ignore_index=True)
    # Append what remains in the positive dataframe (unmatched positive
    # m/z values)
    results = results.append(posData, ignore_index=True)
    if (pandas.__version__ < '0.23.0'):
        # Fix unexpected column sorting from append() in pandas v0.20.3
        # or newer (solved in v0.23.0 with argument "sort=False")
        results = results.reindex(negCol, axis=1)
    results.drop('TotalMean', axis=1, inplace=True)
    # Sort results by m/z and retention time and create the CSV file
    results.sort_values([mzCol, rtCol], inplace=True, kind='mergesort')
    results.to_csv(os.path.join(dst, 'amalgamated.csv'), index=False)
    # Update progress bar
    print_progress_bar(total, total, prefix='Amalgamator progress:')
    # Write the final information in log file and remove handler
    logger.info('Amalgamator completed. Output dataframe has %d rows.\n',
                 len(results.index))
    handler.close()
    logger.removeHandler(handler)


def __hitScore__(srcMZ, targetMZ, srcRT, targetRT, parameters):
    # type: (float, float, float, float, LFParameters) -> float
    """Return the hit score of the target frame for the given source
    frame.

    Keyword Arguments:
        srcMZ      -- source m/z
        targetMZ   -- target m/z
        srcRT      -- source retention time
        targetRT   -- target retention time
        parameters -- LipidFinder's Amalgamator parameters instance
    """
    mzDelta = mz_delta(srcMZ, parameters['mzFixedError'],
                       parameters['mzPPMError'])
    mzDiff = abs(srcMZ - targetMZ)
    rtDelta = rt_delta(parameters['maxRTDiffAdjFrame'])
    rtDiff = abs(srcRT - targetRT)
    return sqrt(min(mzDiff / mzDelta, 1.0) ** 2 \
                + min(rtDiff / rtDelta, 1.0) ** 2)


def __bestMatch__(matches,   # list
                  negMass,   # float
                  pmz,       # pandas.Series
                  negRT,     # float
                  prt,       # pandas.Series
                  parameters # LFParameters
                  ):
    # type: (...) -> int
    """Return the index of 'matches' with the highest hit score for the
    given negative m/z and retention time.

    Keyword Arguments:
        matches    -- list of indexes of positive frames that match the
                      negative m/z and retention time (RT)
        negMass    -- negative m/z
        pmz        -- all positive m/z
        negRT      -- negative RT
        prt        -- all positive RT
        parameters -- LipidFinder's Amalgamator parameters instance
    """
    maxScore = 0.0
    maxScoreIndex = 0
    for ind in matches:
        score = __hitScore__(negMass, pmz[ind], negRT, prt[ind], parameters)
        if (score > maxScore):
            maxScore = score
            maxScoreIndex = ind
    return maxScoreIndex
