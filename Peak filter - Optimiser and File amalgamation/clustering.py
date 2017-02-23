from __future__ import division
from builtins import range
from past.utils import old_div
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as spd


def mzCluster(lfd, saveFrame=True):
    """Mass clusters are groups of m/z artifacts that differ from each other in m/z by a
    mass less than the user defined tolerance. Hierarchical clustering is employed to group
    the ions into the most appropriate groups.

    Mass clusters are assigned an arbitory unique integer identifier.

    Args:
        lfd (LipidFinderData): A LipidFinderData object that has NOT been previously clustered
        saveFrame (bool, optional): Used to indicate if the DataFrame post clustering should be saved in the current LipidFinderData object, for general processing the default of True will be unmodified.

    """
    # sort dataframe accending by mz
    lfd.processedDataFrame.sort_values(by=['MZ'], inplace=True)

    # reset index to allow for index selection
    lfd.processedDataFrame.reset_index(inplace=True, drop=True)

    # report in current frame mass difference to between current frame and
    # next frame (after sorting by mass)
    lfd.processedDataFrame.insert(lfd.firstRepOffset, 'Mass Difference - to next Frame',
                                  lfd.processedDataFrame['MZ'].shift(-1) - lfd.processedDataFrame['MZ'])

    # count total number of records
    totalRecCount = len(lfd.processedDataFrame)
    sectBegin = 0
    sectMin = 49

    currentClusterSectionID = 1

    lfd.processedDataFrame.insert(
        lfd.firstRepOffset + 1, 'mz Cluster Section ID', np.nan)

    while totalRecCount - sectBegin >= sectMin:

        sectEnd = sectBegin + sectMin

        while sectEnd < totalRecCount - 1:

            # If the mass difference to the next frame is greater than the
            # average of the range of the largest mass in the current group and
            # the range of the smallest mass in the next group it must be safe
            # to end the mass section
            if lfd.processedDataFrame['Mass Difference - to next Frame'].ix[sectEnd] > old_div((lfd.mzRange(lfd.processedDataFrame['MZ'].ix[sectEnd]) + lfd.mzRange(lfd.processedDataFrame['MZ'].ix[sectEnd + 1])), 2.0):
                break
            sectEnd += 1

        lfd.processedDataFrame.ix[
            sectBegin:sectEnd, 'mz Cluster Section ID'] = currentClusterSectionID
        currentClusterSectionID += 1
        sectBegin = sectEnd + 1

    lfd.processedDataFrame.ix[sectBegin:totalRecCount,
                              'mz Cluster Section ID'] = currentClusterSectionID

    # How many unique cluster sections (could use max in mz Cluster Section
    # ID' column)
    mzClustSectCount = len(lfd.processedDataFrame[
                           'mz Cluster Section ID'].unique())

    lfd.processedDataFrame.insert(
        lfd.firstRepOffset + 2, 'mz Cluster ID', np.nan)

    # variable for keeping track of cluster ids
    currentMaxClusterID = 0

    for clustSect in range(1, mzClustSectCount + 1):

        # Get max error width in current cluster (based on max mass), this will
        # be used as the cut off for hierachical clustering
        currentMaxMZError = lfd.mzRange(lfd.processedDataFrame['MZ'][
                                        lfd.processedDataFrame['mz Cluster Section ID'] == clustSect].max())

        # Get the masses in the current cluster into the correct format. Each
        # mass is a single item list within a list
        vectorMZ = lfd.processedDataFrame['MZ'][lfd.processedDataFrame[
            'mz Cluster Section ID'] == clustSect].values.reshape((-1, 1))

        # Check if there is only 1 sublist within the list
        if len(vectorMZ) == 1:

            # If 1 item then give it next cluster id and move to next mass
            # section
            lfd.processedDataFrame.ix[lfd.processedDataFrame[
                'mz Cluster Section ID'] == clustSect, 'mz Cluster ID'] = currentMaxClusterID + 1
            currentMaxClusterID += 1
        else:
            # perform hierachical clustering

            # calculate distance between every mass in cluster section
            mzDistMatrix = spd.pdist(vectorMZ)

            # calculate linkage
            mzLinkage = hier.complete(mzDistMatrix)

            # return cluster id at each position to an numpy array, note the
            # addition of the last recorded cluster name prior to this
            # (clustering starts at each time, hence currentMaxClusterID
            # starting at 0
            mzClusters = hier.fcluster(
                mzLinkage, currentMaxMZError, 'distance') + currentMaxClusterID

            # Append these cluster ids to the pandas data frame
            lfd.processedDataFrame.ix[lfd.processedDataFrame[
                'mz Cluster Section ID'] == clustSect, 'mz Cluster ID'] = mzClusters

            # Increment currentMaxClusterID by number of unique clusters in
            # current mass section
            currentMaxClusterID += len(set(mzClusters))

    lfd.processedDataFrame.drop(
        'Mass Difference - to next Frame', inplace=True, axis=1)
    lfd.processedDataFrame.drop('mz Cluster Section ID', inplace=True, axis=1)

    # Cluster ids not in order, keep the ids in order with mz
    mc = lfd.processedDataFrame['mz Cluster ID'].values
    CurClus = 1
    lenTF = len(lfd.processedDataFrame)
    for i in range(lenTF - 1):
        if mc[i + 1] == mc[i]:
            mc[i] = CurClus
        else:
            mc[i] = CurClus
            CurClus += 1
    mc[lenTF - 1] = CurClus

    # Re-sort dataframe ready for feature clustering
    lfd.processedDataFrame.sort_values(
        by=['mz Cluster ID', 'Time', 'MZ'], inplace=True)

    # reset index as we have sorted
    lfd.processedDataFrame.reset_index(inplace=True, drop=True)

    # Increment first rep position by 1
    lfd.incFirstRepOffset()

    if saveFrame:
        lfd.mzClusteredData = lfd.processedDataFrame.copy()


def featureCluster(lfd, saveFrame=True):
    """Feature clusters are groups of contiguous artefacts within the same mass cluster
    where each member is separated by a retention time difference of less than the
    'peakAdjacentFrameMaxRT' parameter. Feature clusters are identified and each assigned
    an arbitory unique integer identifier.

    Args:
        lfd (LipidFinderData): A LipidFinderData object that has NOT been previously feature clustered
        saveFrame (bool, optional): Used to indicate if the DataFrame post clustering should be saved in the current LipidFinderData object, for general processing the default of True will be unmodified.

    """

    lfd.processedDataFrame.insert(lfd.firstRepOffset, 'Time Diff', lfd.processedDataFrame[
                                  'Time'].shift(-1) - lfd.processedDataFrame['Time'])
    lfd.processedDataFrame.insert(
        lfd.firstRepOffset + 1, 'Feature Cluster ID', np.nan)
    td = lfd.processedDataFrame['Time Diff'].values
    mc = lfd.processedDataFrame['mz Cluster ID'].values
    fID = lfd.processedDataFrame['Feature Cluster ID'].values
    CurFeat = 1
    lenTF = len(lfd.processedDataFrame)
    for i in range(lenTF - 1):
        if mc[i + 1] == mc[i] and td[i] <= lfd.peakAdjacentFrameMaxRT:
            fID[i] = CurFeat
        else:
            fID[i] = CurFeat
            CurFeat += 1
    fID[lenTF - 1] = CurFeat

    lfd.processedDataFrame.drop('Time Diff', inplace=True, axis=1)

    # Increment first rep position by 1
    lfd.incFirstRepOffset()

    if saveFrame:
        lfd.featureClusteredData = lfd.processedDataFrame.copy()
