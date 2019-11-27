# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to remove diverse types of contaminants:
    > remove_contaminants():
        Removes straight m/z contaminants included in the contaminants
        CSV file from input data.

    > remove_adducts():
        Retain only the highest intensity adduct from the given list of
        pairs.

    > remove_stacks():
        Detect lipid and contaminant stacks and delete all ions present
        (in lipid stacks the parent is retained). A stack is a series of
        ions differing in m/z by a user-defined fixed mass shift.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import ContaminantRemoval
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> ContaminantRemoval.remove_contaminants(data, parameters)
    >>> ContaminantRemoval.remove_adducts(data, parameters)
    >>> ContaminantRemoval.remove_stacks(data, parameters)
"""

import numpy
import pandas

from LipidFinder._utils import mz_delta, mz_tol_range, rt_tol_range


# Set minimum number of features that integrate a lipid stack
MIN_LIPID_STACK = 4
# Set minimum number of features that integrate a contaminant stack
MIN_CONTAM_STACK = 4


def remove_contaminants(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Remove straight m/z contaminants included in the contaminants CSV
    file from input data.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    mzCol = parameters['mzCol']
    # Read the CSV file with the contaminants information
    if (parameters['polarity'] == 'Negative'):
        contaminants = pandas.read_csv(parameters['negContaminantsCSVPath'])
    else:
        contaminants = pandas.read_csv(parameters['posContaminantsCSVPath'])
    # Remove every frame that matches with a known contaminant
    toRemove = []
    for index, mz in contaminants['MZ'].iteritems():
        minMZ, maxMZ = mz_tol_range(mz, parameters['mzFixedError'],
                                    parameters['mzPPMError'])
        toRemove.extend(data[(minMZ <= data[mzCol])
                             & (data[mzCol] <= maxMZ)].index.tolist())
    if (toRemove):
        # Remove duplicates in the list to avoid errors
        toRemove = list(set(toRemove))
        data.drop('Contaminants removal', labels=toRemove, inplace=True)
        # Reset the index of dataframe after the update
        data.reset_index(inplace=True, drop=True)


def remove_adducts(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Retain only the highest intensity adduct from those found in
    the input data.

    Adduct pairs for the appropriate mode are imported from the adducts
    CSV file. An offset is generated for each given pair based upon
    their mass difference. The m/z column in 'data' is searched for
    pairs differing by this offset with the same retention time (within
    tolerance).

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    firstSampleIndex = parameters['firstSampleIndex'] - 1
    lastSampleIndex = firstSampleIndex \
                      + (parameters['numSamples'] * parameters['numTechReps'])
    # Read the CSV file with the adducts information
    if (parameters['polarity'] == 'Negative'):
        adducts = pandas.read_csv(parameters['negAdductsCSVPath'])
        adductsPairs = parameters['negAdductsPairs']
    else:
        adducts = pandas.read_csv(parameters['posAdductsCSVPath'])
        adductsPairs = parameters['posAdductsPairs']
    mzArray = data[parameters['mzCol']].values
    rtArray = data[parameters['rtCol']].values
    # Get the replicates intensities as a new dataframe
    replicates = data.iloc[:, firstSampleIndex : lastSampleIndex]
    replicates = replicates.apply(
            __rep_adduct_removal__, adductsPairs=adductsPairs,
            adducts=adducts, mzArray=mzArray, rtArray=rtArray,
            parameters=parameters)
    # Overwrite data in original dataframe
    data.iloc[:, firstSampleIndex : lastSampleIndex] = replicates
    # Remove rows with every sample intensity equal to 0
    toRemove = replicates[(replicates == 0).all(axis=1)].index.tolist()
    data.drop('Adducts removal', labels=toRemove, inplace=True)
    # Reset the index of dataframe after the removals
    data.reset_index(inplace=True, drop=True)


def __rep_adduct_removal__(replicate,    # pandas.Series
                           adductsPairs, # list
                           adducts,      # pandas.DataFrame
                           mzArray,      # numpy.array
                           rtArray,      # numpy.array
                           parameters    # LFParameters
                           ):
    # type: (...) -> pandas.Series
    """Detect pairs of adducts in the given sample replicate and set to
    zero the lowest intensity of each pair.

    The m/z and retention time (RT) matches are done within a tolerance.

    Keyword Arguments:
        replicate    -- replicate's intensities
        pairs        -- list of paired adducts
        adducts      -- adducts information
        mzArray      -- sample replicate's m/z values
        rtArray      -- sample replicate's rt values
        parameters   -- LipidFinder's PeakFilter parameters instance
    """
    # Get the index of all intensities that are not zero
    nonZeroIndices = replicate.values.nonzero()[0]
    # Create an array of intensities that are not zero, and the arrays
    # with their corresponding m/z and RT values
    nzIntensities = numpy.copy(replicate.values[nonZeroIndices])
    nzMZ = numpy.copy(mzArray[nonZeroIndices])
    nzRT = numpy.copy(rtArray[nonZeroIndices])
    # Create an array to hold a reference to the first adduct found
    # (default: "")
    adductTags = numpy.empty_like(nonZeroIndices, dtype=str)
    adductTags.fill('')
    for pair in adductsPairs:
        # Get the adducts information of the pair to create the lambda
        # function to calculate the offset of the given mass
        pairInfo = adducts.loc[adducts.iloc[:, 0].isin(pair)]
        # abs() function makes the pairs order insensitive
        get_offset = lambda x: \
                abs(x - (pairInfo.iloc[1, 1] * (x - pairInfo.iloc[0, 2]) /
                         pairInfo.iloc[0, 1] + pairInfo.iloc[1, 2]))
        index = 0
        while (index < (nonZeroIndices.size - 1)):
            # Get first adduct tag of the source index
            tag = adductTags[index]
            # If source frame has been previously tagged and tag is not
            # equal to lower mass species then this frame cannot be
            # considered further so consider next frame
            if (tag and (tag != pair[0])):
                index += 1
                continue
            # Get the m/z and RT of the current frame
            mz = nzMZ[index]
            rt = nzRT[index]
            # The adduct mass for the source mass with current pair
            adductMZ = mz + get_offset(mz)
            # The limits of a mass that could be an adduct of the source
            # mass, including the error tolerance
            minAdductMZ, maxAdductMZ = mz_tol_range(
                    adductMZ, parameters['mzFixedError'],
                    parameters['mzPPMError'])
            # Get the tolerance range for the RT
            minRT, maxRT = rt_tol_range(rt, parameters['maxRTDiffAdjFrame'])
            # Create a dataframe of potential adducts by m/z and RT
            potentialAdducts = numpy.where(
                    (nzMZ >= minAdductMZ) & (nzMZ <= maxAdductMZ)
                    & (nzRT >= minRT) & (nzRT <= maxRT))[0]
            if (potentialAdducts.size > 0):
                # Get the index of the adduct with the closest RT to the
                # subject RT
                adductIndex = potentialAdducts[numpy.absolute(
                        nzRT[potentialAdducts] - rt).argmin()]
                # If the intensity of the adduct is greater than the
                # intensity of the source, set the latter to 0.
                # Otherwise, set the former to 0.
                if (nzIntensities[adductIndex] > nzIntensities[index]):
                    # Go ahead if the adduct has not been tagged yet
                    if (not adductTags[adductIndex]):
                        # Record adduct species
                        adductTags[adductIndex] = pairInfo.iloc[1, 0]
                        if (parameters['adductAddition']):
                            replicate.values[nonZeroIndices[adductIndex]] += \
                                    nzIntensities[index]
                        # Set source intensity to 0
                        replicate.values[nonZeroIndices[index]] = 0
                        # Remove source from each array
                        nonZeroIndices = numpy.delete(nonZeroIndices, index)
                        nzIntensities = numpy.delete(nzIntensities, index)
                        nzMZ = numpy.delete(nzMZ, index)
                        nzRT = numpy.delete(nzRT, index)
                        adductTags = numpy.delete(adductTags, index)
                        # 'index' will now point to the next frame
                        continue
                else:
                    if (not tag):
                        # Record adduct species
                        adductTags[index] = pairInfo.iloc[0, 0]
                    if (parameters['adductAddition']):
                        replicate.values[nonZeroIndices[index]] += \
                                nzIntensities[adductIndex]
                    # Set adduct intensity to 0
                    replicate.values[nonZeroIndices[adductIndex]] = 0
                    # Remove adduct from each array
                    nonZeroIndices = numpy.delete(nonZeroIndices, adductIndex)
                    nzIntensities = numpy.delete(nzIntensities, adductIndex)
                    nzMZ = numpy.delete(nzMZ, adductIndex)
                    nzRT = numpy.delete(nzRT, adductIndex)
                    adductTags = numpy.delete(adductTags, adductIndex)
            index += 1
    return replicate


def remove_stacks(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Detect lipid and contaminant stacks and delete all ions present
    (in lipid stacks the parent is retained).

    A stack is a series of ions differing in m/z by a user-defined fixed
    mass shift. Lipid stacks elute at same retention time (RT) whilst
    contaminant stacks increase their RT as overall m/z increases.
    Firstly, the m/z is checked for a lipid stack and if a stack is
    present, all ions except the parent are deleted and the next m/z is
    checked. If no lipid stack is found then the m/z is checked for
    contaminant stacks. If found, the whole stack including the parent
    is removed. The list of lipid and contaminant stack mass differences
    is imported from the stacks CSV file.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    mzCol = parameters['mzCol']
    rtCol = parameters['rtCol']
    firstSample = parameters['firstSampleIndex'] - 1
    lastSample = firstSample \
                 + (parameters['numSamples'] * parameters['numTechReps'])
    # Read the CSV file with the stacks information
    stacks = pandas.read_csv(parameters['stacksCSVPath'])
    # Separate lipid and contaminant stacks
    lipidStacksMZ = stacks.loc[stacks['Category'] == 'Lipid', 'MZ'].values
    contStacksMZ = stacks.loc[stacks['Category'] == 'Contaminant', 'MZ'].values
    # Build an array with m/z, RT and index to track the frames removed
    # as part of a stack
    array = numpy.stack((data[mzCol].values, data[rtCol].values,
                             data.index.values), axis=-1)
    # Start the loop to find every stack in the dataset
    parentIndex = 0
    toRemove = {'lipid': [], 'contam': []}
    while (parentIndex < (len(array) - 1)):
        parentMZ, parentRT = array[parentIndex, 0:2]
        # Lipid stack removal where m/z and RT have to be an exact match
        # (within tolerance)
        minRT, maxRT = rt_tol_range(parentRT, parameters['maxRTDiffAdjFrame'])
        # Stacks can be of only one type: if a lipid stack is found the
        # parent m/z will not be analysed as part of a contaminant stack
        for stackMZ in lipidStacksMZ:
            stackDiff = 0
            gapCount = 0
            stackList = []
            while (gapCount <= parameters['maxStackGap']):
                # Calculate the expected m/z of the next feature
                stackDiff += stackMZ
                minMZ, maxMZ = mz_tol_range(parentMZ + stackDiff,
                                            parameters['mzFixedError'],
                                            parameters['mzPPMError'])
                matches = numpy.where(
                        (minMZ <= array[:, 0]) & (array[:, 0] <= maxMZ)
                        & (minRT <= array[:, 1]) & (array[:, 1] <= maxRT))[0]
                if (len(matches) == 0):
                    gapCount += 1
                else:
                    # Select the feature with the closest RT to parent
                    if (len(matches) == 1):
                        stackIndex = matches[0]
                    else:
                        stackIndex = matches[numpy.absolute(
                                array[matches, 1] - parentRT).argmin()]
                    # Add the frame as member of the stack
                    stackList.append(stackIndex)
            if (len(stackList) >= MIN_LIPID_STACK):
                # Mark frames to be removed as part of the lipid stack
                indexes = array[stackList, 2].tolist()
                toRemove['lipid'].extend(indexes)
                # Remove frames in stack from array
                array = numpy.delete(array, stackList, axis=0)
                if (parameters['lipidStackAddition']):
                    # Add stack intensities to parent frame
                    data.iloc[parentIndex, firstSample : lastSample] += \
                            data.iloc[indexes, firstSample : lastSample].sum(
                                    axis=0)
                break
        else:
            # No lipid stacks where found: search for contamination
            # stacks where m/z has to be an exact match (within
            # tolerance) and the RT between every two consecutive
            # features is the same (in accordance with the m/z
            # difference)
            for stackMZ in contStacksMZ:
                stackDiff = 0
                stackList = []
                # The RT gap is unknown until the stack has two
                # features, so calculate the minimum RT of next feature
                minRT = rt_tol_range(parentRT,
                                     parameters['maxRTDiffAdjFrame'])[1]
                # Get every possible stack within the maximum gap
                # distance and keep the largest one
                for gapCount in range(0, parameters['maxStackGap']):
                    # Calculate the expected m/z of next stack feature
                    stackDiff += stackMZ
                    minMZ, maxMZ = mz_tol_range(parentMZ + stackDiff,
                                                parameters['mzFixedError'],
                                                parameters['mzPPMError'])
                    matches = numpy.where(
                            (minMZ <= array[:, 0]) & (array[:, 0] <= maxMZ)
                            & (minRT < array[:, 1]))[0]
                    # Explore every possible stack (different RT gap)
                    # and keep the largest one
                    for i in range(0, len(matches)):
                        nextMZ, nextRT = array[matches[i], 0:2]
                        rtGap = (nextRT - parentRT) / (gapCount + 1)
                        matchStackList = _collect_stack(
                                array, matches[i], rtGap, stackMZ, parameters)
                        if (len(matchStackList) > len(stackList)):
                            stackList = matchStackList
                if ((len(stackList) + 1) >= MIN_CONTAM_STACK):
                    # Mark frames to be removed as part of the
                    # contaminant stack (parent included)
                    stackList.append(parentIndex)
                    indexes = array[stackList, 2].tolist()
                    toRemove['contam'].extend(indexes)
                    # Remove frames in stack from array
                    array = numpy.delete(array, stackList, axis=0)
                    # Adjust index after parent has been removed
                    parentIndex -= 1
                    break
        parentIndex += 1
    # Remove lipid and/or contaminant stack features
    if (toRemove['lipid'] or toRemove['contam']):
        data.drop('Stacks removal (lipid)', labels=toRemove['lipid'],
                  inplace=True)
        data.drop('Stacks removal (contaminant)', labels=toRemove['contam'],
                  inplace=True)
        # Reset the index of dataframe after the update
        data.reset_index(inplace=True, drop=True)


def _collect_stack(array, index, rtGap, stackMZ, parameters):
    # type: (numpy.ndarray, int, float, float, LFParameters) -> list
    """Get every feature that matches the stack m/z difference and
    retention time gap from the previous feature to shape the
    contaminant stack.

    Keyword Arguments:
        array      -- array with every feature m/z, retention time (RT)
                      and index
        index      -- index of the previous feature in 'array'
        rtGap      -- RT difference between consecutive features
        stackMZ    -- contaminant m/z difference
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    nextMZ = array[index, 0]
    lastHitRT = array[index, 1]
    rtDiff = 0
    stackList = [index]
    gapCount = 0
    while (gapCount <= parameters['maxStackGap']):
        # Calculate the expected m/z and RT of next stack feature
        nextMZ += stackMZ
        minMZ, maxMZ = mz_tol_range(nextMZ,
                                    parameters['mzFixedError'],
                                    parameters['mzPPMError'])
        rtDiff += rtGap
        expectedRT = lastHitRT + rtDiff
        minRT, maxRT = rt_tol_range(expectedRT, parameters['maxRTDiffAdjFrame'])
        matches = numpy.where(
                (minMZ <= array[:, 0]) & (array[:, 0] <= maxMZ)
                & (minRT <= array[:, 1]) & (array[:, 1] <= maxRT))[0]
        if (len(matches) == 0):
            gapCount += 1
        else:
            # Select the frame with the closest RT to the expected one
            if (len(matches) == 1):
                stackIndex = matches[0]
            else:
                stackIndex = matches[numpy.absolute(
                        array[matches, 1] - expectedRT).argmin()]
            # Add the frame as member of the stack
            stackList.append(stackIndex)
            # Reset the information to calculate the next RT
            lastHitRT = array[stackIndex, 1]
            rtDiff = 0
            # Reset the number of gaps
            gapCount = 0
    return (stackList)
