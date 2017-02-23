from builtins import range
import pandas as pd
import numpy as np


def process(lfd, startIndex, endIndex, repsPerGroup):
    """Removes outliers from a set of replicates on a row by row basis. All replicates may
    be discarded if this does not sufficiently reduce RSD of the remaining replicates.

    Args:
        lfd (LipidFinderData): A LipidFinderData object
        startIndex (int): The index of the first replicate of the group
        endIndex (int): The index of the last replicate of the group
        repsPerGroup (TYPE): The number of replicates in the group

    Returns:
        pandas DataFrame: A copy of the lfd.processedDataFrame that has been outlier corrected
    """
    # repsPerGroup is the number of distinct groups between start and end index

    # Add dummy row to start (to deal with 1st row possibly being run twice
    # for optimisation reasons)
    firstRow = lfd.processedDataFrame.ix[0:0].copy()
    tempDataFrame = pd.concat(
        [firstRow, lfd.processedDataFrame], ignore_index=True)

    # Loop through each set of replicates per sample, in each case slicing out
    # and processing 1 sample's replicates
    for column in range(startIndex, endIndex, repsPerGroup):
        tempDataFrame.ix[:, column:column + repsPerGroup] = tempDataFrame.ix[:,
                                                                             column:column + repsPerGroup].apply(__repsFrameData__, axis=1, lfData=lfd)

    # Drop the dummy row and reset index
    correctedProcessedDataFrame = tempDataFrame.drop(tempDataFrame.index[0])
    correctedProcessedDataFrame.reset_index(inplace=True, drop=True)

    return correctedProcessedDataFrame

# Helper function to return the non-zero mean of an array


def __meanNoZero__(inArray):
    return inArray[np.nonzero(inArray)[0]].mean()

# Helper function to return the non-zero population standard deviation of
# an array


def __sdNoZero__(inArray):
    return inArray[np.nonzero(inArray)[0]].std()


def __repsFrameData__(inArray, lfData):

    # Get a reference to the value array of the passed series
    inty = inArray.values

    # By default we use the lower RSD cut off
    gapFillRSDCutOff = lfData.outlierLowIntensityRSD

    # Number of replicates
    repCount = len(inty)

    # Number of non-Nan values
    repValCount = np.count_nonzero(inty)

    # Number of allowed deletions 1 for 4 or 5 for, 2 for 6 or more
    dels = 0 if repCount <= 3 else 2 if repCount > 5 else 1

    # Check at least 50% of replicates have values
    while 2 * repValCount > repCount:

        curMean = __meanNoZero__(inty)

        if curMean > lfData.outlierHighIntensityValue:

            gapFillRSDCutOff = lfData.outlierHighIntensityRSD

        # RSD check
        # If RSD > allowed:
        if (__sdNoZero__(inty) / curMean * 100) > gapFillRSDCutOff:

            # Can we delete any replicates?
            if dels > 0:

                # Create copy of series
                temp = inty.copy()

                # Max rep index with maximum deviation
                remRep = abs(inty - curMean).argmax()

                # Set set max deviation index to np.nan in series copy
                temp[remRep] = 0

                # Is the largest mean deviation replicate a clear outlier of
                # the remaining replicates
                if abs(inty[remRep] - __meanNoZero__(temp)) > 3 * __sdNoZero__(temp):
                    # Remove the offending replicate
                    inty[remRep] = 0

                    dels -= 1
                else:
                    # No individual replicate is a clear outlier therefore
                    # delete the lot
                    np.copyto(inty, 0)
                    break
            else:
                # Cannot reduce the RSD by removing replicates therefore delete
                # the lot
                np.copyto(inty, 0)
                break

        else:
            # The RSD is within tolerance
            break
    else:
        # delete the lot, not enough non-nan replicates
        np.copyto(inty, 0)

    # As inArray.values == inty we return inArray which has been modified all
    # along
    return inArray
