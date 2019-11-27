# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods focused on creating a summary of the dataset:
    > create_summary():
        Create a summary CSV file containing only the mean sample
        intensity of each feature within the given retention time
        window.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import Summary
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> Summary.create_summary(data, parameters)
"""

import os

import pandas


def create_summary(data, parameters, dst=''):
    """Create a summary CSV file containing only the mean sample
    intensity of each feature within the given retention time window.

    If 'dst' is not an absolute path, the current working directory will
    be used as starting point. If "peakfilter_<polarity>_summary.csv"
    file already exists, it will be overwritten without warning.
    "<polarity>" stands for "positive" or "negative", as stated in the
    parameters.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
        dst        -- destination directory where the file will be saved
                      [default: current working directory]
    """
    # The summary dataframe will contain the ID, m/z, retention time,
    # polarity and samples mean columns. The polarity column is added
    # in case Amalgamator is used afterwards.
    rtCol = parameters['rtCol']
    # Biological sample means are before the isotope annotation
    firstIndex = -parameters['numSamples'] * 2
    lastIndex = firstIndex + parameters['numSamples']
    summaryData = pandas.concat(
            [data.iloc[:, 0], data[parameters['mzCol']], data[rtCol],
             data.iloc[:, firstIndex : lastIndex]],
            axis=1)
    summaryData.insert(3, 'Polarity', parameters['polarity'])
    # Restrict the summary information to the retention time window from
    # the parameters
    rtRange = parameters['rtRange']
    summaryData = summaryData.loc[(rtRange[0] <= summaryData[rtCol])
                                  & (summaryData[rtCol] <= rtRange[1])]
    # Create the CSV file with the summary information in 'dst'
    fileName = 'peakfilter_{0}_summary.csv'.format(
            parameters['polarity'].lower())
    summaryData.to_csv(os.path.join(dst, fileName), index=False)
