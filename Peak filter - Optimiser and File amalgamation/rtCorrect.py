from builtins import range
import pandas as pd
import numpy as np
import sys

lfd = None
rtMeans = None
startIndex = None
endIndex = None
repsPerGroup = None
allRT = None


def process(inLfd, inRtMeans=False):
    """Corrects minor misalignments in retention time at the not picked up by SIEVE or XCMS
    in the original processing or that may have become more apparant during PeakFilter
    processing.

    Args:
        inLfd (LipidFinderData): The LipidFinderData object holding the 'processedDataFrame' for retention time correction
        inRtMeans (bool, optional): A toggle to indicate if the correction should occur at the replicate level or at the sample means level

    Returns:
        pandas DataFrame: A copy of the inLFD.processedDataFrame that has been RT corrected
    """
    global lfd
    lfd = inLfd

    global rtMeans
    rtMeans = inRtMeans

    # repsPerGroup is the number of distinct groups between start and end index
    global startIndex
    startIndex = lfd.firstRepOffset

    global endIndex
    endIndex = lfd.lastRepOffset

    global repsPerGroup
    # repsPerGroup is different depending on where the module is run
    if rtMeans:
        repsPerGroup = len(lfd.meanColList)
    else:
        repsPerGroup = lfd.numberOfTechReps

    global allRt

    # Add dummy cluster to start (to deal with 1st frame possibly being run
    # twice for optimisation reasons
    firstGroup = lfd.processedDataFrame[
        lfd.processedDataFrame['Feature Cluster ID'] == 1].copy()
    firstGroupIndices = firstGroup.index.values
    firstGroup.loc[:, 'Feature Cluster ID'] = 0
    tempDataFrame = pd.concat(
        [firstGroup, lfd.processedDataFrame], ignore_index=True)

    # Set allRT array after dummy cluster added above
    allRt = tempDataFrame['Time'].values

    # Create groupby object on Feature Cluster ID'
    rdGroup = tempDataFrame.groupby(['Feature Cluster ID'])

    # Apply function to each series feature
    reps = rdGroup.apply(__processFeature__)

    # Drop dummy cluster
    reps.drop(firstGroupIndices, inplace=True)
    reps.reset_index(inplace=True, drop=True)

    return reps


def __processFeature__(featClust):

    if rtMeans:

        # The sample means for the cluster
        temp = featClust.ix[:, lfd.meanColList].copy()

        # Get the local index of all frames in one samples replicates that have
        # at least 1 replicate with a non-zero intensity
        nonZeroInd = np.where(temp.sum(axis=1) > 0)[0]

        # Get the global index of all frames in one samples replicates that
        # have at least 1 replicate with a non-zero intensity in the current
        # cluster
        origNonZeroInd = temp.index.values[nonZeroInd]

        # An array of the time difference to next frame (tdnf)
        tdnf = np.roll(allRt[origNonZeroInd], -1) - allRt[origNonZeroInd]

        # No RT correction possible if <2 frames
        if origNonZeroInd.size > 1:

            # Generate a 2d array of 1 sample's replicate intensities where at
            # least 1 non-zero per frame exists
            intys = temp.values[nonZeroInd]

            # Perform RT correction on 1 sample's replicates in the current
            # cluster
            __processSample__(intys, tdnf)

            # Copy the RT corrected intensities for the replicates back to the
            # 1 sample's DataFrame slice
            temp.values[nonZeroInd] = intys

        # Copy the 1 sample DataFrame slice back to the DataFarme feature
        # cluster
        featClust.ix[:, lfd.meanColList] = temp

    else:

        # Loop through each set of replicates per sample, in each case slicing
        # out and processing 1 sample's replicates
        for sampFirstReps in range(startIndex, endIndex, repsPerGroup):

            # One samples replicates for the cluster
            temp = featClust.ix[
                :, sampFirstReps:sampFirstReps + repsPerGroup].copy()

            # Get the local index of all frames in one samples replicates that
            # have at least 1 replicate with a non-zero intensity
            nonZeroInd = np.where(temp.sum(axis=1) > 0)[0]

            # Get the global index of all frames in one samples replicates that
            # have at least 1 replicate with a non-zero intensity in the
            # current cluster
            origNonZeroInd = temp.index.values[nonZeroInd]

            # An array of the time difference to next frame (tdnf)
            tdnf = np.roll(allRt[origNonZeroInd], -1) - allRt[origNonZeroInd]

            # No RT correction possible if <2 frames
            if origNonZeroInd.size > 1:

                # Generate a 2d array of 1 sample's replicate intensities where
                # at least 1 non-zero per frame exists
                intys = temp.values[nonZeroInd]

                # Perform RT correction on 1 sample's replicates in the current
                # cluster
                __processSample__(intys, tdnf)

                # Copy the RT corrected intensities for the replicates back to
                # the 1 sample's DataFrame slice
                temp.values[nonZeroInd] = intys

            # Copy the 1 sample DataFrame slice back to the DataFarme feature
            # cluster
            featClust.ix[:, sampFirstReps:sampFirstReps + repsPerGroup] = temp

    return featClust

# Recursive!!


def __processSample__(intys, tdnf):

    # Take a copy of 2d array (intys) that represents the current feature
    # cluster replicate values, used later to check if changes have been made
    # to intys
    oldIntys = np.copy(intys)

    # Count replicates
    rows = len(intys)

    # Count of frames in cluster
    reps = intys[0].size

    # Loop through each row of each frame
    for rep in range(reps):
        for row in range(rows):

            # If the value of the current frame current replicate is not 0
            if intys[row][rep] != 0:
                continue

            # Are there at least half the intensity values for the current
            # frame non-zero
            elif (2 * np.count_nonzero(intys[row]) >= repsPerGroup):

                # A list to store the values of the cell above (index 0) and
                # the frame below (index 1) - If they exist then 0 will be
                # overwritten with the intensity value in the corresponding
                # frame
                adjacentFrameValues = [0, 0]

                # Is there is a frame above (not the case if currently in the
                # first frame)
                if row > 0:
                    # The intensity of the replicate frame above is not 0 and
                    # the time difference to next frame above is within the
                    # allowable time
                    if intys[row - 1][rep] != 0 and tdnf[row - 1] < lfd.peakAdjacentFrameMaxRT:
                        adjacentFrameValues[0] = intys[row - 1][rep]

                # Is there is a frame below (not the case if currently in the
                # last frame)
                if row < rows - 1:
                    # The intensity of the replicate frame below is not 0 and
                    # the time difference to next frame below is within the
                    # allowable time
                    if intys[row + 1][rep] != 0 and tdnf[row] < lfd.peakAdjacentFrameMaxRT:
                        adjacentFrameValues[1] = intys[row + 1][rep]

                # Check if any valid intensities found either above or below
                if np.sum(adjacentFrameValues) > 0:

                    # Get mean and standard deviation of current frame non-zero values
                    # Don't have to worry about nan caused by taking mean where
                    # there are 0 nonzero values as we only have >0 nonzero
                    # values
                    repMean = intys[row][np.nonzero(intys[row])[0]].mean()
                    repStdDev = intys[row][np.nonzero(intys[row])[0]].std()

                    # Indicates which frame replicate value to shift into
                    # current current frame
                    swapIndex = 0

                    # Need !=0 as repMean - lfd.rtCorrectStDev * repStdDev
                    # might be <0
                    if (adjacentFrameValues[0] != 0 and adjacentFrameValues[0] >= repMean - lfd.rtCorrectStDev * repStdDev and adjacentFrameValues[0] <= repMean + lfd.rtCorrectStDev * repStdDev):

                        # If frame above has less than half non-zero replicates
                        if (2 * np.count_nonzero(intys[row - 1]) < repsPerGroup):

                            swapIndex = -1

                        # Check if frame above has half the valid replicate
                        # intensities, nb possible for frame above to more than
                        # half non-zero frames owing to moving in from frame
                        # above that, if more than half then cannot move
                        elif 2 * np.count_nonzero(intys[row - 1]) == repsPerGroup:

                            # If the mean intensities in the current frame are
                            # >= mean intensities in frame above then we could
                            # shift the value from the frame above
                            if (repMean >= intys[row - 1][np.nonzero(intys[row - 1])[0]].mean()):

                                swapIndex = -1

                    if (adjacentFrameValues[1] != 0 and adjacentFrameValues[1] >= repMean - lfd.rtCorrectStDev * repStdDev and adjacentFrameValues[1] <= repMean + lfd.rtCorrectStDev * repStdDev):

                        # If the swapIndex is not 0 then there is a value in
                        # the frame above that can move as well as below, move
                        # in the value closest to the mean of the current frame
                        if swapIndex == 0:

                            # If frame below has less than half non-zero
                            # replicates
                            if (2 * np.count_nonzero(intys[row + 1]) < repsPerGroup):

                                swapIndex = 1

                            # Check if frame below has half the valid replicate
                            # intensities
                            elif 2 * np.count_nonzero(intys[row + 1]) == repsPerGroup:

                                # If the mean intensities in the current frame
                                # are >= mean intensities in frame above then
                                # we could shift the value from the frame above
                                if (repMean >= intys[row + 1][np.nonzero(intys[row + 1])[0]].mean()):

                                    swapIndex = 1

                        else:
                            # Need to check which which candidate is closer to
                            # the mean of the rest of the samples in the frame
                            if abs(repMean - adjacentFrameValues[1]) < abs(repMean - adjacentFrameValues[0]):

                                # If frame below has less than half non-zero
                                # replicates
                                if (2 * np.count_nonzero(intys[row + 1]) < repsPerGroup):

                                    swapIndex = 1

                                # Check if frame below has half the valid
                                # replicate intensities
                                elif 2 * np.count_nonzero(intys[row + 1]) == repsPerGroup:

                                    # If the mean intensities in the current
                                    # frame are >= mean intensities in frame
                                    # below then we could shift the value from
                                    # the frame below
                                    if (repMean >= intys[row + 1][np.nonzero(intys[row + 1])[0]].mean()):

                                        swapIndex = 1

                    # Make the changes should a swap be necessary
                    if swapIndex != 0:
                        intys[row][rep] = intys[row + swapIndex][rep]
                        intys[row + swapIndex][rep] = 0
    # Recursive call to repeat until no further changes
    if not np.array_equal(intys, oldIntys):
        __processSample__(intys, tdnf)
