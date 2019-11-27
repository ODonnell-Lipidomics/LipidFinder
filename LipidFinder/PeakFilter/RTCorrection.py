# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods focused on correct any retention time misalignments:
    > correct_retention_time():
        Correct minor misalignments in retention time not detected by the
        pre-processing software.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import RTCorrection
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> RTCorrection.correct_retention_time(data, parameters)
"""

import numpy
import pandas

from LipidFinder._py3k import range


def correct_retention_time(data, parameters, means=False):
    # type: (LFDataFrame, LFParameters, bool) -> None
    """Correct minor misalignments in retention time not detected by the
    pre-processing software.

    These misalignments might have become clearer during PeakFilter's
    pipeline.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
        means      -- perform the correction over mean columns instead
                      of each sample replicate? [default: False]
    """
    # Create groupby object on feature clusters and process each one
    data = data.groupby(['FeatureClusterID']).apply(
            __process_feature__, parameters=parameters, means=means)
    if (not means):
        # Drop empty frames (if any)
        data.drop_empty_frames('Empty frames after Retention Time Correction',
                               parameters)


def __process_feature__(featureCluster, parameters, means):
    # type: (pandas.DataFrame, LFParameters, bool) -> pandas.DataFrame
    """Correct retention time misalignment in the given feature cluster.

    Keyword Arguments:
        featureCluster -- frames with the same feature cluster ID
        parameters     -- LipidFinder's PeakFilter parameters instance
        means          -- perform the correction over mean columns
                          instead of each sample replicate?
    """
    if (len(featureCluster) == 1):
        return featureCluster
    if (means):
        # The sample means for the feature cluster
        tmpData = featureCluster.iloc[:, -parameters['numSamples'] : ].copy()
        # Get the index of frames with at least 1 column with a non-zero
        # intensity
        nonZeroIndices = numpy.where(tmpData.sum(axis=1) > 0)[0]
        if (nonZeroIndices.size > 1):
            # Get array of retention times (RT)
            rtArray = featureCluster[parameters['rtCol']].values
            # Get an array of the time difference to next frame
            rtDiff = numpy.roll(rtArray[nonZeroIndices], -1) \
                    - rtArray[nonZeroIndices]
            # Get the array of intensities for the frames with at least
            # 1 column with a non-zero intensity
            intensity = tmpData.values[nonZeroIndices]
            __process_sample__(intensity, rtDiff, parameters,
                               parameters['numSamples'])
            # Replace old values with the new ones
            tmpData.values[nonZeroIndices] = intensity
            featureCluster.iloc[:, -parameters['numSamples'] : ] = tmpData
    else:
        firstSampleIndex = parameters['firstSampleIndex'] - 1
        lastSampleIndex = firstSampleIndex + (parameters['numSamples']
                                              * parameters['numTechReps'])
        # Get array of RTs
        rtArray = featureCluster[parameters['rtCol']].values
        # Loop through each set of replicates per sample
        for firstIndex in range(firstSampleIndex, lastSampleIndex,
                                parameters['numTechReps']):
            lastIndex = firstIndex + parameters['numTechReps']
            tmpData = featureCluster.iloc[:, firstIndex : lastIndex].copy()
            # Get the index of frames with at least 1 replicate with a
            # non-zero intensity
            nonZeroIndices = numpy.where(tmpData.sum(axis=1) > 0)[0]
            if (nonZeroIndices.size > 1):
                # Get an array of the time difference to next frame
                rtDiff = numpy.roll(rtArray[nonZeroIndices], -1) \
                        - rtArray[nonZeroIndices]
                # Get the array of intensities for the frames with at least
                # 1 replicate with a non-zero intensity
                intensity = tmpData.values[nonZeroIndices]
                __process_sample__(intensity, rtDiff, parameters,
                                   parameters['numTechReps'])
                # Replace old values with the new ones
                tmpData.values[nonZeroIndices] = intensity
                featureCluster.iloc[:, firstIndex : lastIndex] = tmpData
    return featureCluster


def __process_sample__(intensity, rtDiff, parameters, repsPerGroup):
    # type: (numpy.ndarray, numpy.ndarray, LFParameters, int) -> None
    """Correct retention time misalignment in the given sample.

    Keyword Arguments:
        intensity    -- intensity per frame and sample's replicate
        rtDiff       -- time differences between consecutive frames
        parameters   -- LipidFinder's PeakFilter parameters instance
        repsPerGroup -- number of replicates per sample
    """
    while True:
        # Copy 'intensity' array to check later if it has been modified
        oldIntensity = numpy.copy(intensity)
        # Number of frames and replicates in the given feature cluster
        numRows, numCols = intensity.shape
        for rep in range(0, numCols):
            for row in range(0, numRows):
                if (intensity[row][rep] != 0):
                    continue
                # Require at least half non-zero intensity values
                elif ((2 * numpy.count_nonzero(intensity[row]))
                      >= repsPerGroup):
                    # Adjacent frame (row -/+ 1) intensity values
                    adjFrameValues = [0, 0]
                    if ((row > 0) and (intensity[row - 1][rep] != 0)
                        and (rtDiff[row - 1]
                             < parameters['maxRTDiffAdjFrame'])):
                        # The frame above has a non-zero intensity and
                        # is within the allowed retention time (RT)
                        # threshold
                        adjFrameValues[0] = intensity[row - 1][rep]
                    if ((row < (numRows - 1)) and (intensity[row + 1][rep] != 0)
                        and (rtDiff[row] < parameters['maxRTDiffAdjFrame'])):
                        # The frame below has a non-zero intensity and
                        # is within the allowed RT threshold
                        adjFrameValues[1] = intensity[row + 1][rep]
                    if (any(adjFrameValues)):
                        # Save the contiguous frame (if any) where to
                        # swap the intensity values
                        swapIndex = 0
                        # At least one contiguous intensity is greater
                        # than zero. Get mean and standard deviation of
                        # current frame (non-zero values).
                        repMean = intensity[row][numpy.nonzero(
                                intensity[row])[0]].mean()
                        repStdDev = intensity[row][numpy.nonzero(
                                intensity[row])[0]].std()
                        # Calculate the maximum standard deviation
                        stDev = parameters['intensityStDev'] * repStdDev
                        if ((adjFrameValues[0] != 0)
                            and (adjFrameValues[0] >= repMean - stDev)
                            and (adjFrameValues[0] <= repMean + stDev)):
                            if ((2 * numpy.count_nonzero(intensity[row - 1]))
                                < repsPerGroup):
                                swapIndex = -1
                            elif ((2 * numpy.count_nonzero(intensity[row - 1]))
                                  == repsPerGroup):
                                prevFrameMean = intensity[row - 1][
                                        numpy.nonzero(intensity[row - 1])[0]
                                        ].mean()
                                if (repMean >= prevFrameMean):
                                    swapIndex = -1
                        if ((adjFrameValues[1] != 0)
                            and (adjFrameValues[1] >= repMean - stDev)
                            and (adjFrameValues[1] <= repMean + stDev)):
                            # If 'swapIndex' is not 0, swap with the
                            # closest intensity value to the mean of the
                            # current frame
                            if ((swapIndex == 0)
                                or ((swapIndex != 0)
                                    and (abs(repMean - adjFrameValues[1])
                                         < abs(repMean - adjFrameValues[0])))):
                                nextNonZeroReps = numpy.count_nonzero(
                                        intensity[row + 1])
                                if ((2 * nextNonZeroReps) < repsPerGroup):
                                    swapIndex = 1
                                elif ((2 * nextNonZeroReps) == repsPerGroup):
                                    nextFrameMean = intensity[row + 1][
                                            numpy.nonzero(intensity[row + 1])[0]
                                            ].mean()
                                    if (repMean >= nextFrameMean):
                                        swapIndex = 1
                        if (swapIndex != 0):
                            # Swap with the chosen contiguous frame
                            intensity[row][rep] = \
                                    intensity[row + swapIndex][rep]
                            intensity[row + swapIndex][rep] = 0
        # Repeat the process until no more modifications are performed
        if (numpy.array_equal(intensity, oldIntensity)):
            break
