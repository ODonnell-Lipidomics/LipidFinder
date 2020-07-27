# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods focused on calculate the sample intensity mean:
    > calculate_sample_means():
        Calculate and add the mean of the intensity of each sample
        replicates.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import SampleMeansCalc
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> SampleMeansCalc.calculate_sample_means(data, parameters)
"""

import re

import numpy

from LipidFinder._py3k import range


def calculate_sample_means(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Calculate and add the mean of the intensity of each sample
    replicates in the input dataframe.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Column index of first sample replicate
    startIndex = parameters['firstSampleIndex'] - 1
    # Column index of last sample replicate
    endIndex = startIndex + (parameters['numSamples']
                             * parameters['numTechReps'])
    if (parameters['numTechReps'] == 1):
        # The mean of a single replicate is the replicate itself, so the
        # mean column will have a copy of the single sample replicate
        for firstIndex in range(startIndex, endIndex):
            colName = data.columns[firstIndex] + '_mean'
            data[colName] = data.iloc[:, firstIndex].astype(float).round(0).astype(int)
    else:
        for firstIndex in range(startIndex, endIndex,
                                parameters['numTechReps']):
            lastIndex = firstIndex + parameters['numTechReps']
            # Create the column name for the mean of the current sample
            colName = re.sub('\d+$', "", data.columns[firstIndex]) + '_mean'
            # Get means (not taking into account zeros) of the sample
            rawMeans = data.iloc[:, firstIndex : lastIndex].apply(
                    lambda x: x.sum() /
                        (x.astype(bool).sum() if (x.astype(bool).sum()) else 1),
                    axis=1)
            # Round to nearest integer, cast to integer and insert
            # sample means into the dataframe
            data[colName] = rawMeans.round(0).astype("int64")
