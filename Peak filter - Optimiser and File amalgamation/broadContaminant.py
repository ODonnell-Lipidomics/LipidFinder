import numpy as np
import pandas as pd

broadContsdMult = None
broadContMinPoints = None
broadContRSDCutOff = None
broadContrtSDCutOff = None
allRt = None


def processAllFeatures(lfd):
    """Removes ions that elute across a chromatogram for a particular m/z (within a tolerance)
    with similar intensities that are likely to be contaminants.

    Args:
        lfd (LipidFinderData): A LipidFinderData object that has been previously clustered and peaks found or is of XCMS origin and has been processed by xcmsReformat.py

    """

    # The number of standard deviations from the mean of all other points for
    # current point to be considered an outlier
    global broadContsdMult
    broadContsdMult = lfd.broadContsdMult

    # Minimum number of non-zero points considered representative enough to
    # allow broad contaminant removal to go ahead
    global broadContMinPoints
    broadContMinPoints = lfd.broadContMinPoints

    # The target intensity RSD for the non outliers to be considerd similar
    # enough for removal
    global broadContRSDCutOff
    broadContRSDCutOff = lfd.broadContRSDCutOff

    # The target RT RSD for the non outliers to be considerd similar enough
    # for removal
    global broadContrtSDCutOff
    broadContrtSDCutOff = lfd.broadContrtSDCutOff

    global allRt

    # Add dummy cluster to start (to deal with 1st frame possibly being run
    # twice for optimisation reasons
    firstGroup = lfd.processedDataFrame[
        lfd.processedDataFrame['mz Cluster ID'] == 1].copy()
    firstGroupIndices = firstGroup.index.values
    firstGroup.loc[:, 'mz Cluster ID'] = 0
    tempDataFrame = pd.concat(
        [firstGroup, lfd.processedDataFrame], ignore_index=True)

    # Set allRT array after dummy cluster added above
    allRt = tempDataFrame['Time'].values

    # Select out just the sample means only
    reps = tempDataFrame.ix[:, lfd.meanColList]

    # Add dummy column as first series apply possible run twice
    tempCol = tempDataFrame.ix[:, lfd.meanColList[0]]
    reps.insert(0, 'Temp', tempCol)

    # Create groupby object on 'mz Cluster ID'
    rdGroup = reps.groupby(tempDataFrame['mz Cluster ID'])

    # Apply function to each series feature
    reps = rdGroup.apply(__processFeature__)

    # Drop dummy cluster
    reps.drop(firstGroupIndices, inplace=True)
    reps.reset_index(inplace=True, drop=True)

    # Drop dummy column
    reps.drop('Temp', inplace=True, axis=1)

    # Overwrite data in original processedDataFrame
    lfd.processedDataFrame.ix[:, lfd.meanColList] = reps

    lfd.broadContaminantCorrected = lfd.processedDataFrame.copy()


def __processFeature__(massClust):
    massClust = massClust.apply(__singleRepProcess__)
    return massClust


def __singleRepProcess__(repMassClust):
    # Check there are at least x nonzeros otherwise don't modify
    if np.count_nonzero(repMassClust.values) > broadContMinPoints:
        inty = np.copy(repMassClust.values)
        inty = __process__(inty)
        np.copyto(repMassClust.values, inty)
    return repMassClust

# Find outliers in inArray where values represented by indices in outliers
# are removed before calculations


def __findOutliers__(inArray, outliers):
    # Find indices of remaining elements in inArray with outliers removed
    index = __allExcept__(inArray, outliers)

    # Values in inArray, values at outlier indices removed
    inArrayNoOut = inArray[index]

    # Mean - x * SD
    lowOut = np.mean(inArrayNoOut) - broadContsdMult * np.std(inArrayNoOut)

    # Mean + x * SD
    highOut = np.mean(inArrayNoOut) + broadContsdMult * np.std(inArrayNoOut)

    # Indices of any outliers < Mean - mult x SD
    indLowOuts = np.where(inArrayNoOut < lowOut)[0]

    # Indices of any outliers > Mean + mult x SD
    indHighOuts = np.where(inArrayNoOut > highOut)[0]

    # Return an array of indices of any low or high outliers
    return np.append(index[indLowOuts], index[indHighOuts])

# Find high outliers in inArray where values represented by indices in
# outliers are removed before calculations


def __findHighOutliers__(inArray, outliers):
    # Find indices of remaining elements in inArray with outliers removed
    index = __allExcept__(inArray, outliers)

    # Values in inArray, values at outlier indices removed
    inArrayNoOut = inArray[index]

    # Mean + 2 x SD
    highOut = np.mean(inArrayNoOut) + broadContsdMult * np.std(inArrayNoOut)

    # Return an array of indices of high outliers
    return index[np.where(inArrayNoOut > highOut)[0]]

# Find RSD of inArray with where values represented by indices in outliers
# are removed before calculation


def __findRSD__(inArray, outliers):
    # Find indices of remaining elements in inArray with outliers removed
    index = __allExcept__(inArray, outliers)

    # Values in inArray, values at outlier indices removed
    inArrayNoOut = inArray[index]

    # Return RSD of inArray where indices represented by outlier are not
    # considered in the calculation
    return round(np.std(inArrayNoOut) / np.mean(inArrayNoOut) * 100, 3)

# Find SD of inArray with where values represented by indices in outliers
# are removed before calculation


def __findSD__(inArray, outliers):

    # Find indices of remaining elements in inArray with outliers removed
    index = __allExcept__(inArray, outliers)

    # Values in inArray, values at outlier indices removed
    inArrayNoOut = inArray[index]

    # Return SD of inArray where indices represented by outlier are not
    # considered in the calculation
    return round(np.std(inArrayNoOut), 3)

# Returns an array of indices of inArray with omitted indices removed, an
# inverse omissions when considering inArray


def __allExcept__(inArray, omissions):
    # Create an array of indices from 0 to the sixe of inArray
    index = np.arange(0, inArray.size)

    # Create boolean mask
    mask = np.ones(inArray.size, dtype=bool)

    # Make indices in mask 0 where in omissions
    mask[omissions] = 0

    # Return indices with omissions removed
    return index[mask]


def __process__(intensities):
    # Get an index of all intensities that are not zero
    nonZeroIndex = intensities.nonzero()[0]

    # Create an array of intensities that are not zero (using the nonZeroIndex)
    nonZeroIntensities = np.copy(intensities[nonZeroIndex])

    # Create an array of corresponding retention time values where intensity
    # is is not zero (using the nonZeroIndex)
    nonZeroTime = np.copy(allRt[nonZeroIndex])

    # Create an array to store the all (high and low) outlier indices
    outliers = np.array([], dtype=int)

    # Create an array to store the just the high outlier indices
    highOutliers = np.copy(outliers)

    # While there are a significant number of non zero intensities left that
    # are not outlier
    while nonZeroIndex.size - outliers.size > broadContMinPoints:

        # If the RSD of the non outliers drop below RSD cut off
        if __findRSD__(nonZeroIntensities, outliers) < broadContRSDCutOff:

            # The remaining non-outlier values are now candidates to be deleted
            # If the SD of time values with outliers removed is less than cut
            # off - (All outliers as low outliers will be removed anyway)
            if __findSD__(nonZeroTime, outliers) > broadContrtSDCutOff:

                # The group is compact enough, delete all non outliers
                # Create an array to hold the new intensity values - Those will
                # be populated with outlier values only
                temp = np.zeros_like(intensities)

                # Populate the corresponding highOutlier indices with the
                # values they represent
                temp[nonZeroIndex[highOutliers]
                     ] = nonZeroIntensities[highOutliers]

                # Copy this array and change the reference of intensities
                intensities = np.copy(temp)

            break

        # The RSD of the non outliers is above RSD cut off
        else:
            # Get any new outliers since outlier lsit was last updated
            newOutliers = __findOutliers__(nonZeroIntensities, outliers)

            # If there are no new outliers
            if newOutliers.size > 0:

                # AGet any new high outliers and append to the high outliers
                # list
                highOutliers = np.append(
                    highOutliers, __findHighOutliers__(nonZeroIntensities, outliers))

                # Add newOutlier to outliers list
                outliers = np.append(outliers, newOutliers)

            # There are no more outliers according to the input parameters, the
            # RSD is above cut off and cannot be improved. There are gretaer
            # values than the nonZero value cut off and again this cannot be
            # improved
            else:
                break
    return intensities
