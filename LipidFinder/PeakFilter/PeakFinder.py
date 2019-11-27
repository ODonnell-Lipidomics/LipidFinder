# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to merge all features of each feature cluster
into one:
    > process_features():
        Combine the lower intensity features into the most intense
        feature (peak centre) within each feature clusters. Wide peaks
        and leading and trailing tails that are indicative of
        contaminants are also removed.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import PeakFinder
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> PeakFinder.process_features(data, parameters)
"""

import numpy
import pandas

from LipidFinder.PeakFilter import Clustering


def process_features(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Combine the intensities of lower intensity features into the most
    intense feature within feature clusters.

    For all features in all replicates the intensities of lower
    intensity features (peak frames) within feature clusters are
    combined into the most intense feature (peak centre) where they are
    part of the same peak. Wide peaks and leading and trailing tails
    that are indicative of contaminants are also removed.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    firstSampleIndex = parameters['firstSampleIndex'] - 1
    lastSampleIndex = firstSampleIndex \
                      + (parameters['numSamples'] * parameters['numTechReps'])
    # Perform m/z and feature clustering
    Clustering.cluster_by_mz(data, parameters)
    Clustering.cluster_by_features(data, parameters)
    featureIDCol = data.columns.values[-1]
    # Add dummy cluster to avoid unexpected behavior when using apply():
    # "In the current implementation, apply calls func twice on the
    # first column/row to decide whether it can take a fast or slow code
    # path."
    firstGroup = data[data[featureIDCol] == 1].copy()
    firstGroupIndices = firstGroup.index.values
    firstGroup.loc[:, featureIDCol] = 0
    tmpData = pandas.concat([firstGroup, data], ignore_index=True)
    # Get array of retention time column
    rtArray = tmpData[parameters['rtCol']].values
    # Select out just the replicates only
    replicates = tmpData.iloc[:, firstSampleIndex : lastSampleIndex]
    # Add dummy column as first Series to avoid unexpected behavior when
    # using apply()
    tempCol = tmpData.iloc[:, firstSampleIndex]
    replicates.insert(0, 'DummyColumn', tempCol)
    # Create groupby object on "Feature Clusters" and process each
    # feature
    replicates = replicates.groupby(tmpData[featureIDCol]).apply(
            __process_feature__, parameters=parameters, rtArray=rtArray)
    # Drop dummy cluster and column
    replicates.drop(firstGroupIndices, inplace=True)
    replicates.drop('DummyColumn', axis=1, inplace=True)
    replicates.reset_index(inplace=True, drop=True)
    # Copy the new samples intensities to data
    data.iloc[:, firstSampleIndex : lastSampleIndex] = replicates
    # Drop empty frames (if any)
    data.drop_empty_frames('Empty frames after Peak Finder', parameters)


def __process_feature__(featureCluster, parameters, rtArray):
    # type: (pandas.DataFrame, LFParameters, numpy.array)
    #       -> pandas.DataFrame
    """Process the feature cluster.

    Keyword Arguments:
        featureCluster -- feature cluster dataframe
        parameters     -- LipidFinder's PeakFilter parameters instance
        rtArray        -- array of retention times from source data
    """
    return featureCluster.apply(__single_rep_feature__, parameters=parameters,
                                rtArray=rtArray)


def __single_rep_feature__(repFeature, parameters, rtArray):
    # type: (pandas.Series, LFParameters, numpy.array) -> pandas.Series
    """Process the each sample replicate of the feature.

    Keyword Arguments:
        repFeature -- sample replicate intensities
        parameters -- LipidFinder's PeakFilter parameters instance
        rtArray    -- array of retention times from source data
    """
    # Count of None/NaN in 'repFeature': if less than 2 then there is
    # either a single frame peak or no peak, so no processing is needed
    if (numpy.count_nonzero(repFeature.values) > 1):
        repRT = rtArray[repFeature.index.values]
        __feat_peak_analysis__(parameters, repFeature.values, repRT)
    return repFeature


def __feat_peak_analysis__(parameters, intensities, repRT):
    # type: (LFParameters, numpy.ndarray, numpy.array) -> None
    """Analyse feature peak.

    Keyword Arguments:
        parameters  -- LipidFinder's PeakFilter parameters instance
        intensities -- array of feature peak intensities
        repRT       -- array of retention times of the sample replicate
    """
    # Index of start of feature (left in to ease code refactoring)
    lowestIndex = 0
    # Index of end of feature
    highestIndex = intensities.size - 1
    # Create array to record the frame category in little-endian
    # unicode-8 (left in to ease code refactoring)
    peakCategory = numpy.empty_like(intensities, dtype='<U2')
    peakCategory.fill('--')
    # Create an array to hold the intensities. Each time an intensity
    # is categorised it is removed, leaving only uncatergorised
    # intensities so the highest can be be selected.
    intensityPeakCat = numpy.copy(intensities)
    # Do while there are uncategorised frames within the feature group.
    # numpy.count_nonzero() counts the number of intensities equal to 0.
    while (sum(peakCategory == '--') > numpy.count_nonzero(intensities == 0)):
        intensityPeakCat[numpy.where(peakCategory != '--')[0]] = -1
        # Get the index of the highest intensity uncategorised and set
        # it as the peak centre
        peakCentreIndex = intensityPeakCat.argmax()
        # Initialize peak start and finish index to peak centre
        peakLowestIndex = peakCentreIndex
        peakHighestIndex = peakCentreIndex
        # Catergorise the peak centre as "PC"
        peakCategory[peakCentreIndex] = 'PC'
        # Retention time at PC (this can change to fairly accommodate
        # wide peaks, but the PC index remains the same)
        rtPC = repRT[peakCentreIndex]
        # Wide peak is where the actual most intense region of a peak is
        # close to the border between two frames, so we can allow two
        # similar frames to be a peak centre together (default: False)
        widePeak = None
        # Indicates whether there are frames to left and right of
        # current peak edges to consider
        framesLeft = True
        framesRight = True
        # Check left side for peakiness
        if ((lowestIndex < peakLowestIndex)
            and (intensities[peakLowestIndex - 1] != 0)):
            decIndex = peakLowestIndex - 1
            if ((parameters['peakMinFoldDiff'] * intensities[decIndex])
                >= intensities[peakLowestIndex]):
                # Potential wide peak
                # Check if there are more frames
                if ((lowestIndex < decIndex)
                    and (intensities[peakLowestIndex - 2] != 0)):
                    dec2Index = peakLowestIndex - 2
                    # Is the intensity at 'peakLowestIndex' greater
                    # than parameters["peakMinFoldDiff"] * (intensity at
                    # 'peakLowestIndex' - 1) - ('peakCentreIndex' - 1)
                    # and ('peakCentreIndex' - 2)
                    if ((parameters['peakMinFoldDiff'] * intensities[dec2Index])
                        >= intensities[decIndex]):
                        # Categorise as solvent feature and check for
                        # further solvents left and right.
                        # 'peakCentreIndex' set to "SF" twice: this
                        # avoids running into another feature if PC is
                        # located at first of last position.
                        # Check for low range of solvent feature
                        peakCategory = __solvents_low_rt__(
                                parameters, peakCategory, intensities,
                                peakCentreIndex, lowestIndex)
                        # Check for high range of solvent feature
                        peakCategory = __solvents_high_rt__(
                                parameters, peakCategory, intensities,
                                peakCentreIndex, highestIndex)
                        # The feature is fully categorised: exit current
                        # iteration
                        continue
                    else:
                        widePeak = "Left"
                        peakLowestIndex -= 2
                else:
                    # No further left frames or the left frame has been
                    # categorised, hence the peak is a wide peak
                    framesLeft = False
                    widePeak = "Left"
                    peakLowestIndex -= 1
            else:
                # The PC is currently a peak
                peakLowestIndex -= 1
        else:
            framesLeft = False
        # Check right side for peakiness
        if ((highestIndex > peakHighestIndex)
            and (intensities[peakHighestIndex + 1] != 0)):
            incIndex = peakHighestIndex + 1
            if ((parameters['peakMinFoldDiff'] * intensities[incIndex])
                >= intensities[peakHighestIndex]):
                if (widePeak):
                    # It is a solvent or the end of the peak. Categorise
                    # as solvent feature and check for further solvents
                    # left and right. 'peakCentreIndex' set to "SF"
                    # twice: this avoids running into another feature if
                    # PC is located at first of last position.
                    # Check for low range of solvent feature
                    peakCategory = __solvents_low_rt__(
                            parameters, peakCategory, intensities,
                            peakCentreIndex, lowestIndex)
                    # Check for high range of solvent feature
                    peakCategory = __solvents_high_rt__(
                            parameters, peakCategory, intensities,
                            peakCentreIndex, highestIndex)
                    # The feature is fully categorised: exit current
                    # iteration
                    continue
                # Potential wide peak
                elif ((highestIndex > incIndex)
                      and (intensities[peakHighestIndex + 2] != 0)):
                    inc2Index = peakHighestIndex + 2
                    # Is the intensity at 'peakHighestIndex' greater
                    # than parameters["peakMinFoldDiff"] * (intensity at
                    # 'peakHighestIndex' + 1) - ('peakCentreIndex' + 1)
                    # and ('peakCentreIndex' + 2)
                    if ((parameters['peakMinFoldDiff'] * intensities[inc2Index])
                        >= intensities[incIndex]):
                        # Categorise as solvent feature and check for
                        # further solvents left and right.
                        # 'peakCentreIndex' set to "SF" twice: this
                        # avoids running into another feature if PC is
                        # located at first of last position.
                        # Check for low range of solvent feature
                        peakCategory = __solvents_low_rt__(
                                parameters, peakCategory, intensities,
                                    peakCentreIndex, lowestIndex)
                        # Check for high range of solvent feature
                        peakCategory = __solvents_high_rt__(
                                parameters, peakCategory, intensities,
                                    peakCentreIndex, highestIndex)
                        # The feature is fully categorised: exit
                        # current iteration
                        continue
                    else:
                        # Right side is a wide peak and there are still
                        # frames on the right and left side of the peak
                        # to categorise
                        widePeak = "Right"
                        peakHighestIndex += 2
                else:
                    # No further right frames or the right frame has
                    # been categorised, hence the peak is a wide peak
                    framesRight = False
                    peakHighestIndex += 1
            else:
                peakHighestIndex += 1
        else:
            framesRight = False
        # Change 'rtPC' to be centre point of current peak
        if (widePeak == "Left"):
            rtPC = (repRT[peakCentreIndex] + repRT[peakCentreIndex - 1]) / 2.0
        elif (widePeak == "Right"):
            rtPC = (repRT[peakCentreIndex] + repRT[peakCentreIndex + 1]) / 2.0
        # There are peak frames to categorise that we found during the
        # process of checking if the PC is part of a peak or solvent
        # feature. "PS" will be used for the rare situation where a
        # frame is part of 2 peaks and needs to be shared amongst both.
        nextIndex = peakHighestIndex + 1
        peakCategory[peakLowestIndex : nextIndex][
                peakCategory[peakLowestIndex : nextIndex] == 'PF'
                ] = 'PS'
        peakCategory[peakLowestIndex : nextIndex][
                (peakCategory[peakLowestIndex : nextIndex] != 'PC')
                & (peakCategory[peakLowestIndex : nextIndex] != 'PS')
                ] = 'PF'
        # Complete left part of the peak
        if (framesLeft):
            # Check there is frame at ('peakLowestIndex' - 1) originally
            # set as 'peakCentreIndex' that may have changed after
            # checking if "PC" is a peak
            while (lowestIndex < peakLowestIndex):
                prevIndex = peakLowestIndex - 1
                if (intensities[prevIndex] == 0):
                    break
                # The frame at ('peakLowestIndex' - 1) is not a
                # previously catergorised "PF": the intensity is not 0
                # and adding the frame will keep the peak within the
                # allowed width. Test each side of peak for time width
                # limits.
                if (round(rtPC - repRT[prevIndex], 3)
                    <= round(parameters['peakMaxRTWidth'] / 2.0, 3)):
                    if ((parameters['peakMinFoldDiff'] * intensities[prevIndex])
                        >= intensities[peakLowestIndex]):
                        if ((parameters['peakMinFoldDiff']
                             * intensities[peakLowestIndex])
                            >= intensities[prevIndex]):
                            # It is a solvent: check for the low range
                            # of the feature
                            peakCategory = __solvents_low_rt__(
                                    parameters, peakCategory, intensities,
                                    prevIndex, lowestIndex)
                        # ('peakLowestIndex' - 1) is either too large to
                        # be a solvent or at least one solvent chain has
                        # been identified: end of left peak
                        break
                    else:
                        peakLowestIndex -= 1
                        if (peakCategory[peakLowestIndex] == 'PF'):
                            # Set as shared frame ("PS")
                            peakCategory[peakLowestIndex] = 'PS'
                        else:
                            # New peak frame
                            peakCategory[peakLowestIndex] = 'PF'
                else:
                    # End of left peak, but there may be frames that
                    # would be in the peak if it was wider (tail
                    # frames). Is the next frame part of the tail of the
                    # last peak? If so, it is a solvent frame. This
                    # stops tails of peaks being categorised as "PC",
                    # avoiding false positives.
                    if ((peakCategory[prevIndex] == '--')
                        and ((parameters['peakMinFoldDiff']
                              * intensities[peakLowestIndex])
                             >= intensities[prevIndex])):
                        # It is a solvent: check for the low range of
                        # the feature
                        peakCategory = __solvents_low_rt__(
                                parameters, peakCategory, intensities,
                                prevIndex, lowestIndex)
                    break
        # Complete right side of the peak
        if (framesRight):
            # Check there is frame at ('peakHighestIndex' - 1) originally
            # set as 'peakCentreIndex' that may have changed after
            # checking if "PC" is a peak
            while (highestIndex > peakHighestIndex):
                nextIndex = peakHighestIndex + 1
                if (intensities[nextIndex] == 0):
                    break
                # The frame at ('peakHighestIndex' + 1) is not a
                # previously catergorised "PF": the intensity is not 0
                # and adding the frame will keep the peak within the
                # allowed width
                if (round(repRT[nextIndex] - rtPC, 3)
                    <= round(parameters['peakMaxRTWidth'] / 2.0, 3)):
                    if ((parameters['peakMinFoldDiff'] * intensities[nextIndex])
                        >= intensities[peakHighestIndex]):
                        if ((parameters['peakMinFoldDiff']
                             * intensities[peakHighestIndex])
                            >= intensities[nextIndex]):
                            # It is a solvent: check for the high range
                            # of the feature
                            peakCategory = __solvents_high_rt__(
                                    parameters, peakCategory, intensities,
                                    nextIndex, highestIndex)
                        # ('peakHighestIndex' - 1) is either too large
                        # to be a solvent or at least one solvent chain
                        # has been identified: end of right peak
                        break
                    else:
                        peakHighestIndex += 1
                        if (peakCategory[peakHighestIndex] == 'PF'):
                            # Set as shared frame ("PS")
                            peakCategory[peakHighestIndex] = 'PS'
                        else:
                            # New peak frame
                            peakCategory[peakHighestIndex] = 'PF'
                else:
                    # End of right peak, but there may be frames that
                    # would be in the peak if it was wider (tail
                    # frames). Is the next frame part of the tail of the
                    # last peak? If so, it is a solvent frame. This
                    # stops tails of peaks being categorised as "PC",
                    # avoiding false positives.
                    if ((peakCategory[nextIndex] == '--')
                        and ((parameters['peakMinFoldDiff']
                              * intensities[peakHighestIndex])
                             >= intensities[nextIndex])):
                        # It is a solvent: check for the high range of
                        # the feature
                        peakCategory = __solvents_high_rt__(
                                parameters, peakCategory, intensities,
                                nextIndex, highestIndex)
                    break
        # Determine peak concatenation type
        nextIndex = peakHighestIndex + 1
        if (parameters['concatAllFrames']):
            # "PC" intensity is set to be the sum of all frames in peak
            intensities[peakCentreIndex] = \
                    intensities[peakLowestIndex : nextIndex].sum()
        else:
            # "PC" intensity is set to be the sum of "PC" and the most
            # intense "PF". If there is only "PC", NaN will be returned
            # and raise an exception. This should not happen since a
            # peak at this point must have at least one "PF" and the
            # "PC".
            if (sum(peakCategory[peakLowestIndex : nextIndex] == 'PF') > 0):
                pfsArray = numpy.where(
                        peakCategory[peakLowestIndex : nextIndex] == 'PF')[0]
                intensities[peakCentreIndex] += \
                        intensities[peakLowestIndex : nextIndex][pfsArray].max()
    # Set all non "PC" frames to 0
    intensities[numpy.where(peakCategory != 'PC')[0]] = 0


def __solvents_low_rt__(parameters,      # LFParameters
                        peakCategory,    # numpy.array
                        intensities,     # numpy.array
                        peakLowestIndex, # int
                        lowestIndex      # int
                        ):
    # type: (...) -> numpy.array
    """Check for the low solvent range of the given feature.

    Keyword Arguments:
        parameters      -- LipidFinder's PeakFilter parameters instance
        peakCategory    -- frame category array
        intensities     -- array of feature peak intensities
        peakLowestIndex -- lowest index of the peak
        lowestIndex     -- lowest index of the feature
    """
    peakCategory[peakLowestIndex] = 'SF'
    while (lowestIndex < peakLowestIndex):
        prevIndex = peakLowestIndex - 1
        if ((peakCategory[prevIndex] == '--')
            and (intensities[prevIndex] != 0)):
            if ((parameters['peakMinFoldDiff'] * intensities[peakLowestIndex])
                >= intensities[prevIndex]):
                peakLowestIndex -= 1
                peakCategory[peakLowestIndex] = 'SF'
                continue
            else:
                # ('peakLowestIndex' - 1) intensity is too large to be a
                # solvent: end of left peak
                break
        else:
            # Either the new frame is already categorised or its
            # intensity is 0: end of the left part of the feature
            break
    return peakCategory


def __solvents_high_rt__(parameters,       # LFParameters
                         peakCategory,     # numpy.array
                         intensities,      # numpy.array
                         peakHighestIndex, # int
                         highestIndex      # int
                         ):
    # type: (...) -> numpy.array
    """Check for the high solvent range of the given feature.

    Keyword Arguments:
        parameters       -- LipidFinder's PeakFilter parameters instance
        peakCategory     -- frame category array
        intensities      -- array of feature peak intensities
        peakHighestIndex -- highest index of the peak
        highestIndex     -- highest index of the feature
    """
    peakCategory[peakHighestIndex] = 'SF'
    while (highestIndex > peakHighestIndex):
        nextIndex = peakHighestIndex + 1
        if ((peakCategory[nextIndex] == '--')
            and (intensities[nextIndex] != 0)):
            if ((parameters['peakMinFoldDiff'] * intensities[peakHighestIndex])
                >= intensities[nextIndex]):
                peakHighestIndex += 1
                peakCategory[peakHighestIndex] = 'SF'
                continue
            else:
                # ('peakLowestIndex' + 1) intensity is too large to be a
                # solvent: end of right peak
                break
        else:
            # Either the new frame is already categorised or its
            # intensity is 0: end of the right part of the feature
            break
    return peakCategory
