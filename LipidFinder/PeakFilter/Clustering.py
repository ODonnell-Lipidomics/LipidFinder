# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to cluster features by diverse criteria:
    > cluster_by_mz():
        Cluster m/z artifacts that differ from each other by a mass less
        than the defined tolerance.

    > cluster_by_features():
        Cluster contiguous ions within the same mass cluster where each
        member is separated by a retention time difference of less than
        'maxRTDiffAdjFrame' (in 'parameters').

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import Clustering
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> Clustering.cluster_by_mz(data, parameters)
    >>> Clustering.cluster_by_features(data, parameters)
"""

import numpy
import pandas
from scipy.cluster import hierarchy
from scipy.spatial import distance

from LipidFinder._utils import mz_delta
from LipidFinder._py3k import range


def cluster_by_mz(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Cluster m/z artifacts that differ from each other by a mass less
    than the defined tolerance.

    Hierarchical clustering is employed to group the ions into the most
    appropriate groups. Mass clusters are assigned an arbitrary unique
    integer identifier.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    firstRepIndex = parameters['firstSampleIndex'] - 1
    mzCol = parameters['mzCol']
    # Create a new dataframe with auxiliary information:
    # "mzDiffNextFrame": m/z difference between current and next frames
    # "mzClusterSectionID": cluster section ID given to each m/z
    auxData = pandas.DataFrame(
            {'mzDiffNextFrame': data[mzCol].shift(-1) - data[mzCol]})
    auxData['mzClusterSectionID'] = numpy.nan
    # Calculate the cluster section ID for each m/z
    numRowsData = len(data)
    sectionBegin = 0
    # Minimum amount of m/z that will belong to the same cluster section
    sectionMinSize = 49
    clusterSectionID = 1
    while ((numRowsData - sectionBegin) >= sectionMinSize):
        sectionEnd = sectionBegin + sectionMinSize
        while (sectionEnd < (numRowsData - 1)):
            # If the m/z difference to the next frame is greater than
            # the sum of the m/z delta of the largest mass in the
            # current group and the smallest mass in the next group, we
            # can close this cluster section and start a new one
            currentDelta = mz_delta(data.loc[sectionEnd, mzCol],
                                    parameters['mzFixedError'],
                                    parameters['mzPPMError'])
            nextDelta = mz_delta(data.loc[sectionEnd + 1, mzCol],
                                 parameters['mzFixedError'],
                                 parameters['mzPPMError'])
            if (auxData.iloc[sectionEnd, 0] > (currentDelta + nextDelta)):
                break
            sectionEnd += 1
        sectionEnd += 1
        auxData.iloc[sectionBegin : sectionEnd, 1] = clusterSectionID
        clusterSectionID += 1
        sectionBegin = sectionEnd
    if (sectionBegin < numRowsData):
        # Group the remaining masses in another cluster section
        auxData.iloc[sectionBegin : numRowsData, 1] = clusterSectionID
    else:
        # The last cluster section ID was not used so get the total
        # number of IDs assigned
        clusterSectionID -= 1
    # Add a column to dataframe where the mass cluster IDs will be saved
    data['mzClusterID'] = numpy.nan
    currentMaxClusterID = 0
    for sectionID in range(1, clusterSectionID + 1):
        sectionRows = auxData.iloc[:, 1] == sectionID
        # Copy the masses in the current cluster into a list of single
        # item lists (one per mass)
        vectorMZ = data.loc[sectionRows, mzCol].values.reshape((-1, 1))
        if (len(vectorMZ) == 1):
            # Give the next cluster ID to the item and move to next
            # cluster section
            currentMaxClusterID += 1
            data.loc[sectionRows, 'mzClusterID'] = currentMaxClusterID
        else:
            # Perform hierarchical clustering:
            # Get maximum m/z error in current cluster (based on maximum
            # m/z). This will be the cut off for hierarchical clustering.
            maxMZ = data.loc[sectionRows, mzCol].max()
            currentMaxMZError = 2 * mz_delta(maxMZ, parameters['mzFixedError'],
                                             parameters['mzPPMError'])
            # Calculate distance between every mass in cluster section
            mzDistMatrix = distance.pdist(vectorMZ)
            # Calculate linkage
            mzLinkage = hierarchy.complete(mzDistMatrix)
            # Return a list of flat cluster IDs for each mass, shifting
            # the numbers by the last assigned cluster ID
            mzClusters = hierarchy.fcluster(mzLinkage, currentMaxMZError,
                                            'distance') + currentMaxClusterID
            # Add this information to the dataframe
            data.loc[sectionRows, 'mzClusterID'] = mzClusters
            # Increment the current cluster ID by the number of unique
            # clusters in the current mass section
            currentMaxClusterID += len(set(mzClusters))
    # Renumber Cluster IDs based on their appearance in the dataframe
    clusterIDs = data['mzClusterID'].values
    id = 1
    numRowsData = len(data)
    for index in range(0, numRowsData - 1):
        clusterIDs[index] = id
        if (clusterIDs[index] != clusterIDs[index + 1]):
            id += 1
    clusterIDs[numRowsData - 1] = id


def cluster_by_features(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Cluster contiguous ions within the same mass cluster where each
    member is separated by a retention time difference of less than
    'maxRTDiffAdjFrame' (in 'parameters').

    Feature clusters are identified and each assigned an arbitrary
    unique integer identifier.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    firstRepIndex = parameters['firstSampleIndex'] - 1
    mzCol = parameters['mzCol']
    rtCol = parameters['rtCol']
    # Re-sort dataframe ready for feature clustering
    data.sort_values(by=['mzClusterID', rtCol, mzCol], inplace=True,
                     kind='mergesort')
    # Reset index
    data.reset_index(inplace=True, drop=True)
    # Create a new dataframe with auxiliary information:
    # "TimeDiff": retention time difference between current and next
    #     frames
    auxData = pandas.DataFrame(
            {'TimeDiff': data[rtCol].shift(-1) - data[rtCol]})
    # Assign a feature cluster ID to each cluster of contiguous
    # ions within the same mass cluster where each member is separated
    # by a retention time difference of less than 'maxRTDiffAdjFrame'
    data['FeatureClusterID'] = numpy.nan
    timeDiffs = auxData['TimeDiff'].values
    mzClusterIDs = data['mzClusterID'].values
    featureClusterIDs = data['FeatureClusterID'].values
    id = 1
    numRowsData = len(data)
    for index in range(0, numRowsData - 1):
        featureClusterIDs[index] = id
        if ((mzClusterIDs[index] != mzClusterIDs[index + 1])
            or (timeDiffs[index] > parameters['maxRTDiffAdjFrame'])):
            id += 1
    featureClusterIDs[numRowsData - 1] = id
