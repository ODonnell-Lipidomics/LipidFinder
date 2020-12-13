# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to target and remove broad contaminants:
    > process_all_features():
        Remove ions that elute across a chromatogram for a particular m/z
        with similar intensities that are likely to be contaminants.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import BroadContaminant
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> BroadContaminant.process_all_features(data, parameters)
"""

import numpy
import pandas


def process_all_features(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Remove ions that elute across a chromatogram for a particular m/z
    with similar intensities that are likely to be contaminants.

    The m/z matches are done within a tolerance.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Add dummy cluster to avoid unexpected behavior when using apply():
    # "In the current implementation, apply calls func twice on the
    # first column/row to decide whether it can take a fast or slow code
    # path."
    firstGroup = data[data['mzClusterID'] == 1].copy()
    firstGroupIndices = firstGroup.index.values
    firstGroup.loc[:, 'mzClusterID'] = 0
    tmpData = pandas.concat([firstGroup, data], ignore_index=True)
    # Get array of retention time column
    rtArray = tmpData[parameters['rtCol']].values
    # Get the sample mean columns as new dataframe
    means = tmpData.iloc[:, -parameters['numSamples'] : ]
    # Add dummy column
    dummyCol = tmpData.iloc[:, -parameters['numSamples']]
    means.insert(0, 'Temp', dummyCol)
    # Create groupby object on Mass Clusters and process each feature
    means = means.groupby(tmpData['mzClusterID']).apply(
            __process_feature__, rtArray=rtArray, parameters=parameters)
    # Drop dummy column
    means.drop('Temp', axis=1, inplace=True)
    # Drop dummy cluster
    means.drop(firstGroupIndices, inplace=True)
    means.reset_index(inplace=True, drop=True)
    # Copy the new sample means to data
    data.iloc[:, -parameters['numSamples'] : ] = means
    # Drop empty frames (if any)
    data.drop_empty_frames('Empty frames after Broad Contaminant Removal',
                           parameters, True)


def __process_feature__(mzCluster, # pandas.DataFrame
                        rtArray,   # numpy.array
                        parameters # LFParameters
                        ):
    # type: (...) -> pandas.DataFrame
    """Process each sample mean of the same mass cluster independently.

    Keyword Arguments:
        mzCluster  -- mass cluster dataframe with sample means
        rtArray    -- array of retention times from source data
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    mzCluster = mzCluster.apply(__process_sample_mean__, rtArray=rtArray[mzCluster.index],
                                parameters=parameters)
    return mzCluster


def __process_sample_mean__(sample, rtArray, parameters):
    # type: (pandas.Series, numpy.array, LFParameters) -> pandas.Series
    """Remove contaminant from the given sample mean intensities.

    Keyword Arguments:
        sample     -- sample mean series
        rtArray    -- array of retention times (RT) from source data
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    if (numpy.count_nonzero(sample.values) > parameters['minNonZeroPoints']):
        intensities = numpy.copy(sample.values)
        # Get the index of every non-zero intensity
        nonZeroIndex = intensities.nonzero()[0]
        # Create two arrays, one with non-zero intensities and the other
        # with their corresponding RT
        nonZeroIntensities = numpy.copy(intensities[nonZeroIndex])
        nonZeroRT = numpy.copy(rtArray[nonZeroIndex])
        # Create an array to store the outlier indices
        outliers = numpy.array([], dtype=int)
        # Create an array to store the just the high outlier indices
        highOutliers = numpy.copy(outliers)
        while ((nonZeroIndex.size - outliers.size)
               > parameters['minNonZeroPoints']):
            # There is a significant number of non-zero intensities not
            # classified as outlier
            if (__get_rsd__(nonZeroIntensities, outliers)
                < parameters['intenRSDCutOff']):
                # The remaining non-outlier values are candidates to be
                # deleted
                if (__get_stdev__(nonZeroRT, outliers) < parameters['rtSDCutOff']):
                    # Group is compact enough: delete all non-outliers
                    tmp = numpy.zeros_like(intensities)
                    # Keep the intensity of 'highOutlier' indices
                    tmp[nonZeroIndex[highOutliers]] = \
                            nonZeroIntensities[highOutliers]
                    intensities = numpy.copy(tmp)
                break
            else:
                newOutliers = __find_outliers__(nonZeroIntensities, outliers,
                                                parameters)
                if (newOutliers.size > 0):
                    # Get every new high outlier
                    newHighOutliers = __find_high_outliers__(
                                nonZeroIntensities, outliers, parameters)
                    highOutliers = numpy.append(highOutliers, newHighOutliers)
                    # Add every new outlier to 'outliers'
                    outliers = numpy.append(outliers, newOutliers)
                else:
                    # There are no more outliers and the results cannot
                    # be improved
                    break
        # Copy the new intensity values to the series
        numpy.copyto(sample.values, intensities)
    return sample


def __find_outliers__(inArray, outliers, parameters):
    # type: (numpy.array, numpy.array, LFParameters) -> numpy.array
    """Return outliers in the input array where values represented by
    input outliers are discarded.

    Keyword Arguments:
        inArray    -- input array
        outliers   -- indices to omit
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Get values of remaining elements in 'inArray' without outliers
    index = __all_except__(inArray, outliers)
    inArrayNoOut = inArray[index]
    mean = numpy.mean(inArrayNoOut)
    delta = parameters['outlierMinDiff'] * numpy.std(inArrayNoOut)
    # Indices of outliers below and above the thresholds
    lowOutliers = numpy.where(inArrayNoOut < (mean - delta))[0]
    highOutliers = numpy.where(inArrayNoOut > (mean + delta))[0]
    # Return an array of indices of both low and high outliers
    return numpy.append(index[lowOutliers], index[highOutliers])


def __find_high_outliers__(inArray, outliers, parameters):
    # type: (numpy.array, numpy.array, LFParameters) -> numpy.array
    """Return high outliers in the input array where values represented
    by input outliers are discarded.

    Keyword Arguments:
        inArray    -- input array
        outliers   -- indices to omit
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Get values of remaining elements in 'inArray' without outliers
    index = __all_except__(inArray, outliers)
    inArrayNoOut = inArray[index]
    mean = numpy.mean(inArrayNoOut)
    delta = parameters['outlierMinDiff'] * numpy.std(inArrayNoOut)
    # Return an array of indices of outliers above the threshold
    return index[numpy.where(inArrayNoOut > (mean + delta))[0]]


def __get_rsd__(inArray, outliers):
    # type: (numpy.array, numpy.array) -> float
    """Return the relative standard deviation of values of input array
    excluding the indices in 'outliers'.

    The returned value has been rounded up to 3 decimal numbers.

    Keyword Arguments:
        inArray  -- input array
        outliers -- indices to omit
    """
    index = __all_except__(inArray, outliers)
    inArrayNoOut = inArray[index]
    return round(numpy.std(inArrayNoOut) / numpy.mean(inArrayNoOut) * 100, 3)


def __get_stdev__(inArray, outliers):
    # type: (numpy.array, numpy.array) -> float
    """Return the standard deviation of values of input array excluding
    the indices in 'outliers'.

    The returned value has been rounded up to 3 decimal numbers.

    Keyword Arguments:
        inArray  -- input array
        outliers -- indices to omit
    """
    index = __all_except__(inArray, outliers)
    return round(numpy.std(inArray[index]), 3)


def __all_except__(inArray, toOmit):
    # type: (numpy.array, numpy.array) -> numpy.array
    """Return an array of indices of the input array without the indices
    to omit.

    Keyword Arguments:
        inArray -- input array
        toOmit  -- indices to omit
    """
    # Create an array of indices from 0 to the size of 'inArray'
    index = numpy.arange(0, inArray.size)
    # Create boolean mask
    mask = numpy.ones(inArray.size, dtype=bool)
    # Set indices contained in 'toOmit' to 0
    mask[toOmit] = 0
    return index[mask]
