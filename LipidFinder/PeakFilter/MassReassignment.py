# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to reassign a m/z value to every feature in a
cluster:
    > reassign_frame_masses():
        Assign each mass in either a mass cluster or feature cluster to
        the mass of the row containing the highest sample mean
        intensity.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import MassReassignment
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> MassReassignment.reassign_frame_masses(data, parameters)
"""

def reassign_frame_masses(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Assign each mass in either a mass cluster or feature cluster to
    the mass of the row containing the highest sample mean intensity.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    if (parameters['featMassAssignment']):
        data = data.groupby('FeatureClusterID').apply(__max_group_mass__,
                                                      parameters=parameters)
    else:
        data = data.groupby('mzClusterID').apply(__max_group_mass__,
                                                 parameters=parameters)


def __max_group_mass__(groupMass, parameters):
    # type: (pandas.DataFrame, LFParameters) -> pandas.DataFrame
    """Replace the m/z value of each row by the m/z with the highest
    sample mean intensity.

    Keyword Arguments:
        groupMass  -- mass or feature cluster
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    maxIntensityIndex = groupMass.iloc[:, -parameters['numSamples'] : ].max(
            axis=1).idxmax()
    mzCol = parameters['mzCol']
    groupMass.loc[:, mzCol] = groupMass[mzCol][maxIntensityIndex].copy()
    return groupMass
