# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to filter common mass defects:
    > remove_salt_clusters():
        Remove features identified as salt clusters based on m/z defect.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import MassDefectFilter
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> MassDefectFilter.remove_salt_clusters(data, parameters)
"""

import pandas


def remove_salt_clusters(data, parameters):
    # type: (LFDataFrame, LFParameters) -> pandas.DataFrame
    """Remove features identified as salt clusters based on m/z defect.

    Only take into account frames under the given retention time (RT)
    threshold, excluding m/z values in the inclusion list.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    mzCol = parameters['mzCol']
    # Filter m/z values under the RT threshold
    tmpData = data.loc[data[parameters['rtCol']] <= parameters['rtCutOff'], :]
    # Read the CSV file with the inclusion list
    if (parameters['polarity'] == 'Negative'):
        inclusionList = pandas.read_csv(parameters['negMassDefectCSVPath'])
    else:
        inclusionList = pandas.read_csv(parameters['posMassDefectCSVPath'])
    # Set matching mass delta (in Daltons)
    mzDelta = parameters['mzDelta']
    # Remove the m/z values in the inclusion list to ensure we keep them
    for index, mz in inclusionList['MZ'].iteritems():
        tmpData = tmpData.loc[(tmpData[mzCol] < (mz - mzDelta))
                              | (tmpData[mzCol] > (mz + mzDelta)), :]
    # Get the decimal part of each m/z value
    mzDefectSeries = tmpData[mzCol].apply(
            lambda x: round(x % 1, len(repr(x).split('.')[1])))
    fitMZSeries = tmpData[mzCol].apply(lambda x: 0.00112 * x + 0.01953)
    # Remove rows in 'tmpData' where 'mzDefectSeries' > 'fitMZSeries'
    toRemove = mzDefectSeries.loc[(mzDefectSeries > fitMZSeries)].index.values
    data.drop('Mass defect filter', labels=toRemove, inplace=True)
    # Reset the index of dataframe after the removals
    data.reset_index(inplace=True, drop=True)
