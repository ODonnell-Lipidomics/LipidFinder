# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods focused on adjust sample replicate intensities that
may be altered during the experiment:
    > remove_solvent_effect():
        Remove the effect of solvent samples on biological samples.

    > remove_low_intensity_frames():
        Discard features (rows) where intensity is below the set
        threshold.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import SolventCalcs
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> SolventCalcs.remove_solvent_effect(data, parameters)
    >>> SolventCalcs.remove_low_intensity_frames(data, parameters)
"""

import numpy
import re

import pandas

from LipidFinder.PeakFilter import OutlierCorrection


def remove_solvent_effect(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Remove the effect of solvent samples on biological samples.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Column index of first solvent sample
    firstIndex = parameters['firstSampleIndex'] \
                 + (parameters['numSamples'] * parameters['numTechReps']) \
                 + parameters['numQCReps'] - 1
    # Column index of last solvent sample
    lastIndex = firstIndex + parameters['numSolventReps']
    if (parameters['numSolventReps'] >= 3):
        # Outlier correction of solvent samples
        OutlierCorrection.remove_outliers(data, parameters, src='blanks')
    # Insert mean colvent intensity column into the dataframe
    solMeanCol = re.sub('\d+$', "", data.columns[firstIndex]) + '_mean'
    # Get means (not taking into account zeros) of solvent samples
    solMeans = data.iloc[:, firstIndex : lastIndex].apply(
        lambda x: x[numpy.where(x>0)[0]].mean(), axis=1).fillna(0)
    # Round to nearest integer and convert means to integer types
    data[solMeanCol] = solMeans
    # Subtracts solvent mean intensity from each sample replicate
    firstIndex = parameters['firstSampleIndex'] - 1
    lastIndex = firstIndex \
                + (parameters['numSamples'] * parameters['numTechReps'])
    # Solvent fold elimination: remove frames where all technical
    # replicates of all samples are less than "solventMinFoldDiff"
    # parameter times the solvent mean
    toRemove = numpy.where(data.iloc[:, firstIndex : lastIndex].max(axis=1)
                           < parameters['solventMinFoldDiff'] * data[solMeanCol]
                          )[0]
    data.drop('Solvent removal', labels=toRemove, inplace=True)
    data.reset_index(inplace=True, drop=True)
    # Remove solvent from all remaining samples
    data.iloc[:, firstIndex : lastIndex] = numpy.maximum(0.0,
            data.iloc[:, firstIndex : lastIndex].sub(data[solMeanCol], axis=0))
    # Drop empty frames (if any)
    data.drop_empty_frames('Empty frames after Solvent removal', parameters)


def remove_low_intensity_frames(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Discard features (rows) where intensity is below the threshold.

    The threshold is set by "intenSignifCutOff" parameter.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    firstIndex = parameters['firstSampleIndex'] - 1
    lastIndex = firstIndex \
                + (parameters['numSamples'] * parameters['numTechReps'])
    # Set all replicate intensities that are less than
    # "intenSignifCutOff" to zero
    temp = data.iloc[:, firstIndex : lastIndex]
    temp[temp < parameters['intenSignifCutOff']] = 0
    data.iloc[:, firstIndex : lastIndex] = temp
    # Remove features where all sample replicates are zero
    data.drop_empty_frames('Background correction', parameters)
