# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to target and remove outliers:
    > remove_outliers():
        Removes outliers from a set of replicates on a row by row basis
        until the relative standard deviation is reduced below the
        established threshold.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import OutlierCorrection
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> OutlierCorrection.remove_outliers(data, parameters)
"""

import numpy
import pandas

from LipidFinder._py3k import range


def remove_outliers(data, parameters, src='samples'):
    # type: (LFDataFrame, LFParameters) -> None
    """Removes outliers from a set of replicates on a row by row basis.

    All sample replicates may be discarded if the relative standard
    deviation (RSD) of the remaining replicates cannot be reduced below
    the established threshold.

    Keyword Arguments:
        data         -- LFDataFrame instance
        parameters   -- LipidFinder's PeakFilter parameters instance
        src          -- columns where to check for outliers: "samples"
                        or "blanks" [default: "samples"]
    """
    if (src not in ['samples', 'blanks']):
        raise ValueError('Unexpected value. Options: samples, blanks')
    # Set the corresponding values regarding the columns to evaluate
    if (src == 'samples'):
        startIndex = parameters['firstSampleIndex'] - 1
        endIndex = startIndex + (parameters['numSamples']
                                 * parameters['numTechReps'])
        repsPerGroup = parameters['numTechReps']
    else:
        startIndex = parameters['firstSampleIndex'] \
                + (parameters['numSamples'] * parameters['numTechReps']) \
                + parameters['numQCReps'] - 1
        endIndex = startIndex + parameters['numSolventReps']
        repsPerGroup = parameters['numSolventReps']
    # Add dummy row to avoid unexpected behavior when using apply(): "In
    # the current implementation, apply calls func twice on the first
    # column/row to decide whether it can take a fast or slow code
    # path."
    tmpData = data.iloc[0, :].to_frame().transpose()
    tmpData = tmpData.append(data, ignore_index=True)
    # Loop through each set of replicates per sample, in each case
    # slicing out and processing 1 sample's replicate
    for firstIndex in range(startIndex, endIndex, repsPerGroup):
        lastIndex = firstIndex + repsPerGroup
        tmpData.iloc[:, firstIndex : lastIndex] = \
                tmpData.iloc[:, firstIndex : lastIndex].apply(
                        __reps_frame__, axis=1, parameters=parameters)
    # Copy to data the new replicates values after removing the first
    # dummy row
    tmpData = tmpData.iloc[1 : ]
    tmpData.index = tmpData.index - 1
    data.iloc[:, startIndex : endIndex] = \
            tmpData.iloc[:, startIndex : endIndex]
    # Drop empty frames (if any)
    data.drop_empty_frames('Empty frames after Outlier Correction', parameters)


def __non_zero_mean__(inArray):
    # type: (numpy.ndarray) -> float
    """Return the mean of non-zero values of an array.

    Keyword Arguments:
        inArray -- intensities of one frame
    """
    return inArray[numpy.nonzero(inArray)[0]].mean()


def __non_zero_std__(inArray):
    # type: (numpy.ndarray) -> float
    """Return the standard deviation of non-zero values of an array.

    Keyword Arguments:
        inArray -- intensities of one frame
    """
    return inArray[numpy.nonzero(inArray)[0]].std()


def __reps_frame__(inArray, parameters):
    # type: (pandas.Series, LFParameters) -> float
    """Remove any value out of the parameters' thresholds.

    Keyword Arguments:
        inArray    -- intensities of one frame
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Get a reference to the value array of the passed series
    intensities = inArray.values
    # By default use the lower relative standard deviation (RSD) cut off
    gapFillRSDCutOff = parameters['intensityRSD'][0]
    # Number of replicates
    repCount = len(intensities)
    # Number of non-zero values
    repValCount = numpy.count_nonzero(intensities)
    if (repCount <= 3):
        dels = 0
    elif (repCount > 5):
        dels = 2
    else:
        dels = 1
    # Check at least half of replicates have values
    while ((2 * repValCount) > repCount):
        curMean = __non_zero_mean__(intensities)
        if (curMean > parameters['intenOutlierCutOff']):
            gapFillRSDCutOff = parameters['intensityRSD'][1]
        # RSD check
        if ((__non_zero_std__(intensities) / curMean * 100) > gapFillRSDCutOff):
            # Can we delete any replicates?
            if (dels > 0):
                # Create copy of series
                temp = intensities.copy()
                # Set maximum deviation replicate to zero in series copy
                remRep = abs(intensities - curMean).argmax()
                temp[remRep] = 0
                # Is the largest mean deviation replicate a clear
                # outlier of the remaining replicates?
                if (abs(intensities[remRep] - __non_zero_mean__(temp))
                    > (3 * __non_zero_std__(temp))):
                    # Remove the outlier replicate
                    intensities[remRep] = 0
                    dels -= 1
                else:
                    # No individual replicate is a clear outlier: delete
                    # the lot
                    numpy.copyto(intensities, 0)
                    break
            else:
                # Cannot reduce the RSD by removing sample replicates:
                # delete the lot
                numpy.copyto(intensities, 0)
                break
        else:
            # RSD is within tolerance
            break
    else:
        # Delete the lot, not enough non-zero replicates
        numpy.copyto(intensities, 0)
    return inArray
