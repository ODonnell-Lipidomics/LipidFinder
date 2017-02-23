from __future__ import division
from past.utils import old_div
import numpy as np
import pandas as pd

peakMinFoldCutOff = None
peakMaxRTWidth = None
peakConcatenateAllFrames = None
allRt = None


def processAllFeatures(lfd, saveFrame=True):
    """For all features in all replicates the intensities of lower intensity features
    (peak frames) within feature clusters are combined into the most intense feature
    (peak centre) where they are part of the same peak. Wide peaks and leading and trailing
    tails that are indicative of contaminants are also removed.

    Args:
        lfd (LipidFinderData): A LipidFinderData object that has been previously feature clustered
        saveFrame (bool, optional): Should the 'processedDataFrame' be archived following solvents calculations.

    """
    global peakMinFoldCutOff
    peakMinFoldCutOff = lfd.peakMinFoldCutOff

    global peakMaxRTWidth
    peakMaxRTWidth = lfd.peakMaxRTWidth

    global peakConcatenateAllFrames
    peakConcatenateAllFrames = lfd.peakConcatenateAllFrames

    global allRt

    # Column index of first sample first replicate
    indFirst = lfd.firstRepOffset

    # Column index of last sample last replicate
    indLast = lfd.lastRepOffset

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

    # Select out just the replicates only
    reps = tempDataFrame.ix[:, indFirst:indLast]

    # Add dummy column as first series apply possible run twice
    tempCol = tempDataFrame.ix[:, indFirst]
    reps.insert(0, 'Temp', tempCol)

    # Create groupby object on Feature Cluster ID'
    rdGroup = reps.groupby(tempDataFrame['Feature Cluster ID'])

    # Apply function to each series feature
    reps = rdGroup.apply(__processFeature__)

    # Drop dummy cluster
    reps.drop(firstGroupIndices, inplace=True)
    reps.reset_index(inplace=True, drop=True)

    # Drop dummy column
    reps.drop('Temp', inplace=True, axis=1)

    # Overwrite data in original processedDataFrame
    lfd.processedDataFrame.ix[:, indFirst:indLast] = reps

    if saveFrame:
        lfd.peakFound = lfd.processedDataFrame.copy()


def __processFeature__(featClust):
    featClust = featClust.apply(__singleRepFeature__)
    return featClust


def __singleRepFeature__(repFeatClust):

    # Count of None nan in repMassGroup - if less than 2 then either single
    # frame peak or no peak so no processing needed
    if np.count_nonzero(repFeatClust.values) > 1:

        rt = allRt[repFeatClust.index.values]

        __featPeakAnalysis__(repFeatClust.values, rt)

    return repFeatClust


def __featPeakAnalysis__(inty, rt):

    # index of start of feature (will =0, but left in for ease of refactoring)
    featureClusterLowestRTIndex = 0

    # index of end of feature
    featureClusterHighestRTIndex = inty.size - 1

    # Create array to record the frame category (consider chnge to fc,
    # currently left in for ease of refactoring)
    mpc = np.empty_like(inty, dtype='<U2')
    mpc.fill('--')

    # Create an array to hold the intensities, each time an intensity is
    # categorisd it is remove thus leaving the only uncatergorised intensities
    # so the highest can be be selected
    intPeakCat = np.copy(inty)

    # Do while there are uncategorised frames within the feature group (i.e.
    # 'Molecule Peak Category' mpc == "--") for none nan intensities.
    # np.count_nonzero works because we have a boolean mask whereby 0 becomes
    # True or 1 and then this is counted as a non-zero
    while sum(mpc == '--') > np.count_nonzero(inty == 0):

        intPeakCat[np.where(mpc != '--')[0]] = -1

        # Find highest intensity uncategorised index and set as the peak centre
        peakCentreIndex = intPeakCat.argmax()

        # set peak start and finish index to be same as peak centre to begin
        # with
        peakLowestRTIndex = peakHighestRTIndex = peakCentreIndex

        # Catergorise the peak centre as PC
        mpc[peakCentreIndex] = 'PC'

        # Retention time at PC, this can change to fairly accomodate wide
        # peaks, but the PC index remains the same
        rtPC = rt[peakCentreIndex]

        # Default is false - wide peak is where a the actual most intense
        # region of a peak is close to the border between 2 frames, so we can
        # allow 2 similar frames to be a peak centre together, see line above
        widePeak = None

        # Indicates whether there are frames to left and right of current peak
        # extremities to consider
        framesLeft = True
        framesRight = True

        # Check left side for peakiness
        # Check there is frame at (peakLowestRTIndex - 1), currently same as
        # (peakCentreIndex - 1)
        if featureClusterLowestRTIndex < peakLowestRTIndex:

            # The frame at (peakLowestRTIndex - 1) intensity is not 0
            if inty[peakLowestRTIndex - 1] != 0:

                # Is the intensity at peakLowestRTIndex greater than
                # peakMinFoldCutOff * (intensity at peakLowestRTIndex - 1) -
                # peakCentreIndex and (peakCentreIndex - 1)
                if peakMinFoldCutOff * inty[peakLowestRTIndex - 1] >= inty[peakLowestRTIndex]:

                    # Potential wide peak, check there is frame at
                    # (peakLowestRTIndex - 2), currently same as
                    # (peakCentreIndex - 2)
                    if featureClusterLowestRTIndex < peakLowestRTIndex - 1:

                        # The frame at (peakLowestRTIndex - 2) intensity is not
                        # 0
                        if inty[peakLowestRTIndex - 2] != 0:

                            # Is the intensity at peakLowestRTIndex greater
                            # than peakMinFoldCutOff * (intensity at
                            # peakLowestRTIndex - 1) - (peakCentreIndex - 1)
                            # and (peakCentreIndex - 2)
                            if peakMinFoldCutOff * inty[peakLowestRTIndex - 2] >= inty[peakLowestRTIndex - 1]:

                                # Categorise as solvent feature and check for further solvents left and right, peakCentreIndex set to 'SF' twice, this avoids running into another feature if PC is located at first of last position
                                # Check for low range of solvent feature
                                mpc = __solventsLowRT__(
                                    mpc, inty, peakCentreIndex)

                                # Check for high range of solvent feature
                                mpc = __solventsHighRT__(
                                    mpc, inty, peakCentreIndex, featureClusterHighestRTIndex)

                                # The feature is fully categorised, exit
                                # current iteration while loop
                                continue

                            else:
                                widePeak = "Left"
                                peakLowestRTIndex -= 2

                        # Frame is already catergorised therefore currently a
                        # 'wide peak' and no further valid frames to left
                        else:
                            framesLeft = False
                            widePeak = "Left"
                            peakLowestRTIndex -= 1

                    # No further frames to the left therefore peak is currently
                    # a 'wide peak'
                    else:
                        framesLeft = False
                        widePeak = "Left"
                        peakLowestRTIndex -= 1

                # The PC is currently a peak
                else:
                    peakLowestRTIndex -= 1
            else:
                framesLeft = False
        else:
            framesLeft = False

        # Check right side for peakiness
        # Check there is frame at (peakLowestRTIndex + 1), currently same as
        # (peakCentreIndex + 1)
        if featureClusterHighestRTIndex > peakHighestRTIndex:

            # The frame at (peakHighestRTIndex + 1) intensity is not 0
            if inty[peakHighestRTIndex + 1] != 0:

                # Is the intensity at peakHighestRTIndex greater than
                # peakMinFoldCutOff * (intensity at peakHighestRTIndex + 1) -
                # peakCentreIndex and (peakCentreIndex + 1)
                if peakMinFoldCutOff * inty[peakHighestRTIndex + 1] >= inty[peakHighestRTIndex]:

                    # If 'wide peak' already then there is no peak (i.e. it's
                    # solvent or end of the peak)
                    if widePeak:

                        # Categorise as solvent feature and check for further solvents left and right, peakCentreIndex set to 'SF' twice, this avoids running into another feature if PC is located at first of last position
                        # Check for low range of solvent feature
                        mpc = __solventsLowRT__(mpc, inty, peakCentreIndex)

                        # Check for high range of solvent feature
                        mpc = __solventsHighRT__(
                            mpc, inty, peakCentreIndex, featureClusterHighestRTIndex)

                        # The feature is fully categorised, exit current
                        # iteration while loop
                        continue

                    # Potential wide peak, check there is frame at
                    # (peakHighestRTIndex + 2), currently same as
                    # (peakCentreIndex + 2)
                    elif featureClusterHighestRTIndex > (peakHighestRTIndex + 1):

                        # The frame at (peakHighestRTIndex + 2) intensity is
                        # not 0
                        if inty[peakHighestRTIndex + 2] != 0:

                            # Is the intensity at peakHighestRTIndex greater
                            # than peakMinFoldCutOff * (intensity at
                            # peakHighestRTIndex + 1) - (peakCentreIndex + 1)
                            # and (peakCentreIndex + 2)
                            if (peakMinFoldCutOff * inty[peakHighestRTIndex + 2]) >= inty[peakHighestRTIndex + 1]:

                                # Categorise as solvent feature and check for further solvents left and right, peakCentreIndex set to 'SF' twice, this avoids running into another feature if PC is located at first of last position
                                # Check for low range of solvent feature
                                mpc = __solventsLowRT__(
                                    mpc, inty, peakCentreIndex)

                                # Check for high range of solvent feature
                                mpc = __solventsHighRT__(
                                    mpc, inty, peakCentreIndex, featureClusterHighestRTIndex)

                                # The feature is fully categorised, exit
                                # current iteration while loop
                                continue

                            # It is a wide peak
                            else:
                                widePeak = "Right"
                                peakHighestRTIndex += 2
                                # Right side is a wide peak and there are
                                # frames on the right and left side of the peak
                                # still to categorise

                        # No further valid frames to right
                        else:
                            framesRight = False
                            peakHighestRTIndex += 1

                    # No further frames to the right therefore peak is
                    # currently a 'wide peak'
                    else:
                        framesRight = False
                        peakHighestRTIndex += 1
                else:
                    peakHighestRTIndex += 1

            else:
                framesRight = False
        else:
            framesRight = False

        # Change rtPC to be centre point of current peak
        if widePeak == "Left":
            rtPC = old_div(
                (rt[peakCentreIndex] + rt[peakCentreIndex - 1]), 2.0)
        elif widePeak == "Right":
            rtPC = old_div(
                (rt[peakCentreIndex] + rt[peakCentreIndex + 1]), 2.0)

        # If we get here then we have peak frames to categorise that we found
        # during the process of checking if the PC is part of a peak or solvent
        # feature. 'PS' will be used for the rare situation where a frame is
        # part of 2 peaks and needs to be shared amongst the 2 peaks
        mpc[peakLowestRTIndex:peakHighestRTIndex +
            1][mpc[peakLowestRTIndex:peakHighestRTIndex + 1] == 'PF'] = 'PS'
        mpc[peakLowestRTIndex:peakHighestRTIndex + 1][(mpc[peakLowestRTIndex:peakHighestRTIndex + 1] != 'PC') & (
            mpc[peakLowestRTIndex:peakHighestRTIndex + 1] != 'PS')] = 'PF'

        # Complete peak left
        if framesLeft:

            # Check there is frame at (peakLowestRTIndex - 1) -
            # peakLowestRTIndex - originally set == peakCentreIndex, but may
            # have changed following checking if 'PC' is a peak
            while featureClusterLowestRTIndex < peakLowestRTIndex:

                # The frame at (peakLowestRTIndex - 1) is not a previously catergorised PF, the intensity is not 0 and will adding the frame keep the peak within the permitted width, remember test each side of peak for time width (it gets 50% of total width each - hence peakMaxRTWidth/2.0)
                # https://docs.python.org/2/tutorial/floatingpoint.html -
                # rounding
                if inty[peakLowestRTIndex - 1] == 0:
                    break

                if round(rtPC - rt[peakLowestRTIndex - 1], 3) <= round(old_div(peakMaxRTWidth, 2.0), 3):

                    # Is the intensity at peakLowestRTIndex greater than
                    # peakMinFoldCutOff * (intensity at peakLowestRTIndex - 1)
                    if peakMinFoldCutOff * inty[peakLowestRTIndex - 1] >= inty[peakLowestRTIndex]:

                        # is it low enough to be solvent
                        if peakMinFoldCutOff * inty[peakLowestRTIndex] >= inty[peakLowestRTIndex - 1]:

                            # It is solvent, so now check for low range of
                            # solvent feature
                            mpc = __solventsLowRT__(
                                mpc, inty, peakLowestRTIndex - 1)

                        # (peakLowestRTIndex - 1) is either too large to be solvent so end of left peak OR any solvent chains are identified
                        break
                    else:
                        peakLowestRTIndex -= 1

                        # Check if frame already 'PF' if so then set as shared
                        # frame ('PS')
                        if mpc[peakLowestRTIndex] == 'PF':
                            mpc[peakLowestRTIndex] = 'PS'
                        # Otherwise it's a new peak frame
                        else:
                            mpc[peakLowestRTIndex] = 'PF'

                # no more frames to add, end of left peak, but may be frames
                # that would be in peak if peak were wider, tail frames, we
                # need to ensure they cannot be assigned as 'PC' themselves
                else:
                    # Is the next frame part of the tail of the last Peak, if
                    # so it's solvent frame. This stops the tails of peaks
                    # being catergorised as 'PC' giving false positives. If
                    # their has been any previous catergorisation then no need
                    # to bother and wouldn't want to overwrite 'PF' from
                    # another frame with 'SF'
                    if mpc[peakLowestRTIndex - 1] == '--' and peakMinFoldCutOff * inty[peakLowestRTIndex] >= inty[peakLowestRTIndex - 1]:
                        # It is solvent, so now check for low range of solvent
                        # feature
                        mpc = __solventsLowRT__(
                            mpc, inty, peakLowestRTIndex - 1)
                    break

        # Complete peak right
        if framesRight:

            # Check there is frame at (peakHighestRTIndex + 1) -
            # peakHighestRTIndex - originally set == peakCentreIndex, but may
            # have changed following checking if 'PC' is a peak
            while featureClusterHighestRTIndex > peakHighestRTIndex:

                # The frame at (peakHighestRTIndex + 1) is not a previously catergorised PF, the intensity is not 0 and will adding the frame keep the peak within the permitted width
                # https://docs.python.org/2/tutorial/floatingpoint.html -
                # rounding
                if inty[peakHighestRTIndex + 1] == 0:
                    break
                if round(rt[peakHighestRTIndex + 1] - rtPC, 3) <= round(old_div(peakMaxRTWidth, 2.0), 3):

                    # Is the intensity at peakHighestRTIndex greater than
                    # peakMinFoldCutOff * (intensity at peakHighestRTIndex + 1)
                    if peakMinFoldCutOff * inty[peakHighestRTIndex + 1] >= inty[peakHighestRTIndex]:

                        # is it solvent
                        if peakMinFoldCutOff * inty[peakHighestRTIndex] >= inty[peakHighestRTIndex + 1]:

                            # Check for high range of solvent feature
                            mpc = __solventsHighRT__(
                                mpc, inty, peakHighestRTIndex + 1, featureClusterHighestRTIndex)

                        # (peakHighestRTIndex + 1) is either too large to be solvent so end of left peak OR any solvent chains are identified
                        break
                    else:
                        peakHighestRTIndex += 1

                        # Check if frame already 'PF' if so then set as shared
                        # frame ('PS')
                        if mpc[peakHighestRTIndex] == 'PF':

                            mpc[peakHighestRTIndex] = 'PS'

                        # Otherwise it's a new peak frame
                        else:
                            mpc[peakHighestRTIndex] = 'PF'

                # no more frames to add, end of right peak, but may be frames
                # that would be in peak if peak were wider, tail frames, we
                # need to ensure they cannot be assigned as 'PC' themselves
                else:
                    # Is the next frame part of the tail of the last Peak, if
                    # so it's solvent frame. This stops the tails of peaks
                    # being catergorised as 'PC' giving false positives. If
                    # their has been any previous catergorisation then no need
                    # to bother and wouldn't want to overwrite 'PF' from
                    # another frame with 'SF'
                    if mpc[peakHighestRTIndex + 1] == '--' and peakMinFoldCutOff * inty[peakHighestRTIndex] >= inty[peakHighestRTIndex + 1]:
                        # It is solvent, so now check for low range of solvent
                        # feature
                        mpc = __solventsHighRT__(
                            mpc, inty, peakHighestRTIndex + 1, featureClusterHighestRTIndex)
                    break

        # Determine peak concatenation type
        if peakConcatenateAllFrames:

            # The 'PC' intensity is set to be the sum of all frames in the peak
            inty[peakCentreIndex] = inty[
                peakLowestRTIndex:peakHighestRTIndex + 1].sum()

        else:

            # The 'PC' intensity is set to be the sum of the 'PC' and the most
            # intense 'PF', if there is only 'PC' and this code runs at this
            # point a nan will return and cause an error. In theory this cannot
            # happen here since a Peak here must consist of at least 1 'PF' and
            # the 'PC'
            if sum(mpc[peakLowestRTIndex:peakHighestRTIndex + 1] == 'PF') > 0:
                inty[peakCentreIndex] += inty[peakLowestRTIndex:peakHighestRTIndex +
                                              1][np.where(mpc[peakLowestRTIndex:peakHighestRTIndex + 1] == 'PF')[0]].max()

    # Set all non PC frames to 0 prior to return
    inty[np.where(mpc != 'PC')[0]] = 0
    # return inty


def __solventsLowRT__(mpc, inty, peakLowestRTIndex):

    mpc[peakLowestRTIndex] = 'SF'

    while 0 < peakLowestRTIndex:

        # Is new frame already categorised, if 'PF', 'PC' or 'PS'? we wouldn't
        # overwrite and if SF no point in continuing AND does not have 0
        # intensity
        if mpc[peakLowestRTIndex - 1] == '--' and inty[peakLowestRTIndex - 1] != 0:

            # Is the intensity at (peakLowestRTIndex - 1) less than
            # peakMinFoldCutOff * (intensity at peakLowestRTIndex)
            if peakMinFoldCutOff * inty[peakLowestRTIndex] >= inty[peakLowestRTIndex - 1]:

                peakLowestRTIndex -= 1

                mpc[peakLowestRTIndex] = 'SF'

                continue

            # (peakLowestRTIndex - 1) intensity is too large to be solvent, end of left peak
            else:
                break

        # either the new frame is already categorised or it's intensity is 0
        # meaning the end of the feature to the left
        else:
            break

    return mpc


def __solventsHighRT__(mpc, inty, peakHighestRTIndex, featureClusterHighestRTIndex):

    mpc[peakHighestRTIndex] = 'SF'

    while featureClusterHighestRTIndex > peakHighestRTIndex:

        # Is new frame already categorised, if 'PF', 'PC' or 'PS'? we wouldn't
        # overwrite and if SF no point in continuing AND does not have 0
        # intensity
        if mpc[peakHighestRTIndex + 1] == '--' and inty[peakHighestRTIndex + 1] != 0:

            # Is the intensity at (peakHighestRTIndex + 1) less than
            # peakMinFoldCutOff * (intensity at peakHighestRTIndex)
            if peakMinFoldCutOff * inty[peakHighestRTIndex] >= inty[peakHighestRTIndex + 1]:

                peakHighestRTIndex += 1

                mpc[peakHighestRTIndex] = 'SF'

                continue

            # (peakLowestRTIndex + 1) intensity is too large to be solvent, end of left peak
            else:
                break

        # either the new frame is already categorised or it's intensity is 0
        # meaning the end of the feature to the left
        else:
            break
    return mpc
