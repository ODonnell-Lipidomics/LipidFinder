# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods focused on the quality control (QC) samples of an
experiment:
    > qc_rsd_ratio():
        Calculate the ratio (%) of quality control (QC) samples with a
        relative standard deviation (RSD) lower than the lower cut off
        to QC samples with RSD lower than the upper cut off.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import QCCalcs
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> QCCalcs.qc_rsd_ratio(data, parameters)
"""

import re


def qc_rsd_ratio(data, parameters):
    # type: (LFDataFrame, LFParameters) -> float
    """Return ratio (%) of quality control (QC) samples with a relative
    standard deviation (RSD) lower than QCRSD's lower part to QC samples
    with RSD lower than QCRSD's upper part.

    The mean and RSD of the set of QC samples for each row is calculated
    and added at the end of the input dataframe.
    "QCRSD" is the key of one of the parameters in "parameters".

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Column index of first QC sample
    firstIndex = parameters['firstSampleIndex'] \
                 + (parameters['numSamples'] * parameters['numTechReps']) - 1
    # Column index of last QC sample
    lastIndex = firstIndex + parameters['numQCReps']
    # Get the common QC samples' column name from the first QC column
    colName = re.sub('\d+$', '', data.iloc[:, firstIndex].name)
    meanColName = colName + '_mean'
    rsdColName = colName + '_RSD'
    # Calculate and insert in the input dataframe the mean and RSD of QC
    # samples (row-wise).
    data[meanColName] = data.iloc[:, firstIndex : lastIndex].mean(axis=1)
    # pandas.DataFrame.std() returns sample standard deviation over
    # requested axis, so we need to change the degrees of freedom (ddof)
    # to 0 to get the population standard deviation.
    data[rsdColName] = data.iloc[:, firstIndex : lastIndex].std(axis=1, ddof=0)\
                       * 100.0 / data[meanColName]
    # Return the ratio of QC samples with RSD less than QCRSD's lower
    # part to QC samples with RSD less than QCRSD's upper part.
    lowerRSDCount = (data[rsdColName] < parameters['QCRSD'][0]).sum()
    upperRSDCount = (data[rsdColName] < parameters['QCRSD'][1]).sum()
    return round(lowerRSDCount / float(upperRSDCount) * 100, 1)
