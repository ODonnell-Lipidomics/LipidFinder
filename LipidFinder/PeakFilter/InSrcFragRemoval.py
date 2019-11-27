# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to remove in-source fragments:
    > remove_in_src_frags():
        Remove in-source ion fragments with a retention time tolerance
        of 0.05 minutes.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import InSrcFragRemoval
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> InSrcFragRemoval.remove_in_src_frags(data, parameters)
"""

import numpy
import pandas

from LipidFinder._py3k import range, viewkeys
from LipidFinder._utils import mz_tol_range, rt_tol_range


RT_TOLERANCE = 0.05


def remove_in_src_frags(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Remove in-source ion fragments.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Read the CSV file with the know in-source ion fragments
    if (parameters['polarity'] == 'Negative'):
        fragments = pandas.read_csv(parameters['negIonFragsCSVPath'])
    else:
        fragments = pandas.read_csv(parameters['posIonFragsCSVPath'])
    # Create an array from 'data' with one m/z, retention time (RT) and
    # index per row
    array = numpy.stack((data[parameters['mzCol']].values,
                         data[parameters['rtCol']].values,
                         data.index.values),
                        axis=-1)
    # Get the in-source fragments to be removed
    frags = fragments[fragments['Type'].str.lower() == 'fragment']
    if (frags['MZ'].count != 0):
        # Get the indexes of the features that correspond to common
        # in-source fragments
        fragsIndex = rm_full_frags(array, frags, parameters)
        # Change each index's data type from 'float64' to 'int'
        rmIndexes = array[fragsIndex, 2].astype('int')
        # Remove detected in-source fragments
        array = numpy.delete(array, fragsIndex, axis=0)
        data.drop('In-source fragment removal (fragments)', labels=rmIndexes,
                  inplace=True)
    # Get the neutral losses to detect the in-source fragments
    losses = fragments[fragments['Type'].str.lower() == 'neutral loss']
    if (losses['MZ'].count != 0):
        # Get the indexes of the features that correspond to in-source
        # fragments of common neutral losses
        fragsIndex = rm_neutral_loss_frags(array, losses, parameters)
        # Change each index's data type from 'float64' to 'int'
        rmIndexes = array[fragsIndex, 2].astype('int')
        # Remove detected in-source fragments. If we plan to use 'array'
        # afterwards, we need to remove them from it too:
        #     array = numpy.delete(array, fragsIndex, axis=0)
        data.drop('In-source fragment removal (neutral loss)', labels=rmIndexes,
                  inplace=True)
    # Reset the index of dataframe after the removal
    data.reset_index(inplace=True, drop=True)


def rm_full_frags(array,     # type: numpy.ndarray
                  fragments, # type: pandas.DataFrame
                  parameters # type: LFParameters
                  ):
    # type: (...) -> list[float]
    """Return an index list corresponding to common in-source fragments
    in the given sample array.

    Return the index list of all 'array' features that match the m/z
    values provided in 'fragments' for which there is at least another
    feature above the given m/z cut-off at the same retention time (RT).
    All m/z and RT matching are computed within tolerance.

    Keyword arguments:
        array      -- array with m/z, RT and index of the original
                      dataframe
        fragments  -- in-source fragments to be removed
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Create an array with one in-source fragment m/z cut-off and m/z
    # offset per row
    fragsArray = numpy.stack(
            (fragments['MZ'].values, fragments['MZCutOff'].values),
            axis=-1)
    fragsIndex = []
    for fragMZ, mzCutOff in fragsArray:
        mzRange = mz_tol_range(fragMZ, parameters['mzFixedError'],
                               parameters['mzPPMError'])
        # Get the index of 'array' features that match the in-source
        # fragment m/z value
        mzMatches = numpy.searchsorted(array[:, 0], mzRange)
        if (mzMatches[0] == mzMatches[1]):
            continue
        for index in range(mzMatches[0], mzMatches[1]):
            # To be a match, each feature must have the same RT
            minRT, maxRT = rt_tol_range(array[index, 1], RT_TOLERANCE)
            rtMatches = numpy.where(
                    (array[:, 0] >= mzCutOff) & (array[:, 1] >= minRT)
                    & (array[:, 1] <= maxRT))[0]
            if (len(rtMatches) > 0):
                # Mark the feature as an in-source fragment
                fragsIndex.append(index)
    return fragsIndex


def rm_neutral_loss_frags(array,     # type: numpy.ndarray
                          losses,    # type: pandas.DataFrame
                          parameters # type: LFParameters
                          ):
    # type: (...) -> list[float]
    """Return an index list corresponding to the features in the given
    sample array that have been fragmented.

    Return the index list of all 'array' features that have lost one of
    the m/z in 'losses' and their complete counterpart is present in the
    data. The features to be removed must be higher than the given
    cut-off. All m/z and retention time (RT) matching are computed
    within tolerance.

    Keyword arguments:
        array      -- array with m/z, RT and index of the original
                      dataframe
        losses     -- neutral losses to subtract in order to detect
                      fragmented features
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Create an array with one m/z cut-off and neutral loss m/z per row
    fragsArray = numpy.stack((losses['MZCutOff'].values, losses['MZ'].values),
                             axis=-1)
    # Create a dictionary with cut-off values as keys and their
    # corresponding neutral loss m/z in lists as values
    fragsDict = {}
    for mzCutOff, mzLoss in fragsArray:
        fragsDict.setdefault(mzCutOff, []).append(mzLoss)
    matchIndexSet = set()
    for mzCutOff in viewkeys(fragsDict):
        # Get the index of the first m/z value in 'array' greater than
        # the m/z cut-off
        firstIndex = numpy.searchsorted(array[:, 0], mzCutOff)
        for index in range(firstIndex, len(array)):
            for mzLoss in fragsDict[mzCutOff]:
                # Look for in-source fragments, that is, features that
                # are the result of subtracting the neutral loss to the
                # parent's m/z and elute at the same RT
                fragMZ = array[index, 0] - mzLoss
                mzRange = mz_tol_range(fragMZ, parameters['mzFixedError'],
                                       parameters['mzPPMError'])
                # Get first and last indexes of the features within the
                # m/z range
                mzMatches = numpy.searchsorted(array[:, 0], mzRange)
                if (mzMatches[0] == mzMatches[1]):
                    continue
                # In order to be considered a match, each feature must
                # have the same RT
                minRT, maxRT = rt_tol_range(array[index, 1], RT_TOLERANCE)
                rtMatches = numpy.where(
                        (array[mzMatches[0] : mzMatches[1], 1] >= minRT)
                        & (array[mzMatches[0] : mzMatches[1], 1] <= maxRT)
                        )[0]
                if (len(rtMatches) == 0):
                    continue
                # The resultant indexes are based on the starting index
                # of the search ('mzMatches[0]')
                rtMatches += mzMatches[0]
                # The union of sets will handle any index repetition
                matchIndexSet.update(set(rtMatches))
    return list(matchIndexSet)
