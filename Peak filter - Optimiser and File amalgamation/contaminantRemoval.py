from __future__ import division
from past.utils import old_div
import pandas as pd
import numpy as np


def contaminantsRemoval(lfd):
    """Removes straight m/z contaminants for the appropriate ionisation mode when a match
    (within a small tolerance is found in the contaminants.csv file

    Args:
        lfd (LipidFinderData): A LipidFinderData object that has been previously clustered
        and peaks found or is of XCMS origin and has been processed by xcmsReformat.py

    """
    # Column index of first sample first replicate
    indFirst = lfd.firstRepOffset

    # Column index of last sample last replicate
    indLast = lfd.lastRepOffset

    fullContaminantFrame = pd.read_table('contaminants.csv', sep=',')
    if lfd.filePolarityMode == 'P':
        filtContaminantFrame = fullContaminantFrame[
            fullContaminantFrame['Polarity'] != 'Neg']
    else:
        filtContaminantFrame = fullContaminantFrame[
            fullContaminantFrame['Polarity'] != 'Pos']
    cmz = filtContaminantFrame['MZ'].values

    # Make a copy of lfd.processedDataFrame
    tempDataFrame = lfd.processedDataFrame.copy()

    mz = tempDataFrame['MZ'].values

    # Select out just the replicates only
    reps = tempDataFrame.ix[:, indFirst:indLast]

    # Add dummy column as first series apply possible run twice
    tempCol = tempDataFrame.ix[:, indFirst]
    reps.insert(0, 'Temp', tempCol)

    reps = reps.apply(__singleRepContaminantRemoval__, mz=mz, cmz=cmz,
                      lowerMZLimit=lfd.lowerMZLimit, upperMZLimit=lfd.upperMZLimit)

    # Drop dummy column
    reps.drop('Temp', inplace=True, axis=1)

    # Overwrite data in original processedDataFrame
    lfd.processedDataFrame.ix[:, indFirst:indLast] = reps

    lfd.contaminantsRemoved = lfd.processedDataFrame.copy()


def __singleRepContaminantRemoval__(replicate, mz, cmz, lowerMZLimit, upperMZLimit):
    for mass in cmz:
        contIndList = np.where(
            ((mz >= lowerMZLimit(mass)) & (mz <= upperMZLimit(mass))))[0]
        replicate.ix[contIndList] = 0
    return replicate


def adductsRemoval(lfd):
    """Adduct pairs for the appropriate mode are imported from adducts.csv. For each pair
    an offset is generated based upon the mass difference between them. The m/z list in the
    'processedDataFrame' is searched for pairs differening in m/z by this mass and having
    the same retention time. Where found the highest intensity adduct is retained.

    Args:
        lfd (LipidFinderData): A LipidFinderData object that has been previously clustered
        and peaks found or is of XCMS origin and has been processed by xcmsReformat.py

    """
    adductFrame = pd.read_table('adducts.csv', sep=',', index_col=0)

    ao = adductFrame['Adduct Offset'].values
    mo = adductFrame['Mass Offset'].values

    # Column index of first sample first replicate
    indFirst = lfd.firstRepOffset

    # Column index of last sample last replicate
    indLast = lfd.lastRepOffset

    if lfd.filePolarityMode == 'N':
        # Get negative mode adduct pairs relating to index in adducts.csv file
        adductList = lfd.negativeModeAdductList

    # Must be 'P' if not 'N'
    else:
        # Get positive mode adduct pairs relating to index in adducts.csv file
        adductList = lfd.positiveModeAdductList

    # Make a copy of lfd.processedDataFrame
    tempDataFrame = lfd.processedDataFrame.copy()

    mz = tempDataFrame['MZ'].values
    rt = tempDataFrame['Time'].values

    # Select out just the replicates only
    reps = tempDataFrame.ix[:, indFirst:indLast]

    # Add dummy column as first series apply possible run twice
    tempCol = tempDataFrame.ix[:, indFirst]
    reps.insert(0, 'Temp', tempCol)

    reps = reps.apply(__singleRepAdductRemoval__, adductList=adductList, ao=ao, mo=mo, mz=mz, rt=rt, lowerMZLimit=lfd.lowerMZLimit,
                      upperMZLimit=lfd.upperMZLimit, mzRange=lfd.mzRange, lowerRTLimit=lfd.lowerRTLimit, upperRTLimit=lfd.upperRTLimit, adductAddition=lfd.adductAddition)

    # Drop dummy column
    reps.drop('Temp', inplace=True, axis=1)

    # Overwrite data in original processedDataFrame
    lfd.processedDataFrame.ix[:, indFirst:indLast] = reps

    lfd.adductsRemoved = lfd.processedDataFrame.copy()


def __singleRepAdductRemoval__(replicate, adductList, ao, mo, mz, rt, lowerMZLimit, upperMZLimit, mzRange, lowerRTLimit, upperRTLimit, adductAddition):

    # Get an index of all intensities that are not zero
    nonZeroIndex = replicate.values.nonzero()[0]

    # Create an array of intensities that are not zero (using the nonZeroIndex)
    nzInty = np.copy(replicate.values[nonZeroIndex])

    # Create an array of corresponding mz values where intensity is is not
    # zero (using the nonZeroIndex)
    nzMZ = np.copy(mz[nonZeroIndex])

    # Create an array of corresponding retention time values where intensity
    # is is not zero (using the nonZeroIndex)
    nzRT = np.copy(rt[nonZeroIndex])

    # Create an array to hold a reference to the first adduct found, default -
    # 1. This will be populated by the index of the adduct found from the
    # adducts.csv file
    fa = np.zeros_like(nonZeroIndex) - 1

    for adduct in adductList:
        massDiff = ao[adduct[1]] - ao[adduct[0]]
        lowerMassOffset = mo[adduct[0]]
        higherMassOffset = mo[adduct[1]]

        # Using a while loop instead of a for loop as the number of iterations
        # will alter as the loops progress because of deletions from tempFrame
        # as adducts are found. No need to check the last frame (hence the -1)
        sourceInd = 0
        while sourceInd < nonZeroIndex.size - 1:

            # Get First Adduct tag of the source index
            sourceTag = fa[sourceInd]

            # If source frame has been previously tagged and tag is not equal
            # to lower mass species then this frame cannot be considered
            # further so consider next frame
            if (sourceTag != -1 and sourceTag != adduct[0]):
                sourceInd += 1
                continue

            # The mass of the current frame
            sourceMass = nzMZ[sourceInd]

            # The adduct mass for the source mass with current adduct. The abs
            # ensures correct offset whatever the order in the csv file,
            # otherwise offsets would need to be sorted ascending
            adductMass = sourceMass + \
                abs(((higherMassOffset * sourceMass) -
                     (lowerMassOffset * sourceMass) + massDiff))

            # The lower limit of a mass that could be an adduct of the source
            # mass, both adduct mass and source mass could be in error, this
            # represents the limit of this error
            # lower mass limit for adduct - max low error on source mass
            addLowMassLim = lowerMZLimit(
                adductMass) - (old_div(mzRange(sourceMass), 2.0))

            # The upper limit of a mass that could be an adduct of the source
            # mass, both adduct mass and source mass could be in error, this
            # represents the limit of this error
            # upper mass limit for adduct + max high error on source mass
            addHighMassLim = upperMZLimit(
                adductMass) - (old_div(mzRange(sourceMass), 2.0))

            # Create a dataframe of potential adducts by mass ready to be
            # tested by retention time
            potAdduct = np.where(
                ((nzMZ >= addLowMassLim) & (nzMZ <= addHighMassLim)))[0]

            if potAdduct.size > 0:

                # The retention time of the mass
                massRet = nzRT[sourceInd]

                # Create a dataframe with candidate adducts, the maximum error
                # for the subject retention time and adduct retention time is
                # adjacentFrameMaxTimeDist/2 each (the list size will rarely be
                # above 1 but it can happen. In most cases it will be 0)
                candAdduct = potAdduct[np.where((nzRT[potAdduct] >= lowerRTLimit(
                    massRet)) & (nzRT[potAdduct] <= upperRTLimit(massRet)))[0]]

                # Check if there are any candidate adducts by retention time
                if candAdduct.size > 0:

                    # The index of the adduct with the closest retention time
                    # to the subject retention time (usually only 1 anyway)
                    indAdduct = candAdduct[np.absolute(
                        nzRT[candAdduct] - massRet).argmin()]

                    # If the intensity of the adduct is greater than the source
                    # then set the source frame intensity to 0, otherwise set
                    # the adduct frame intensity to 0
                    if nzInty[indAdduct] > nzInty[sourceInd]:

                        # Make sure there is no tag already on this frame
                        # (!=-1), if there is do nothing
                        if fa[indAdduct] == -1:

                            # Record adduct species in 'First Adduct' field in
                            # adduct frame
                            fa[indAdduct] = adduct[1]

                            # Add intensity of frame to drop to frame to keep
                            # if toggled
                            if adductAddition:
                                # Update replicate
                                replicate.values[nonZeroIndex[
                                    indAdduct]] += nzInty[sourceInd]

                            # Set source index replicate intensity to 0
                            replicate.values[nonZeroIndex[sourceInd]] = 0

                            # Remove the new intensity 0 index from each array
                            nonZeroIndex = np.delete(nonZeroIndex, sourceInd)
                            nzInty = np.delete(nzInty, sourceInd)
                            nzMZ = np.delete(nzMZ, sourceInd)
                            nzRT = np.delete(nzRT, sourceInd)
                            fa = np.delete(fa, sourceInd)

                            # To prevent sourceInd being increment later before
                            # loop back
                            continue

                    else:

                        # Check if 'First Adduct' tag is already populated, if
                        # not then this adduct species is the first adduct
                        if sourceTag == -1:

                            # Record adduct species in 'First Adduct' field in
                            # source frame
                            fa[sourceInd] = adduct[0]

                        # Add intensity of frame to drop to frame to keep if
                        # toggled
                        if adductAddition:
                            # Update replicate
                            replicate.values[nonZeroIndex[
                                sourceInd]] += nzInty[indAdduct]

                        # Set source index replicate intensity to 0
                        replicate.values[nonZeroIndex[indAdduct]] = 0

                        # Remove the new intensity 0 index from each array
                        nonZeroIndex = np.delete(nonZeroIndex, indAdduct)
                        nzInty = np.delete(nzInty, indAdduct)
                        nzMZ = np.delete(nzMZ, indAdduct)
                        nzRT = np.delete(nzRT, indAdduct)
                        fa = np.delete(fa, indAdduct)

            sourceInd += 1

    return replicate


def stacksRemoval(lfd):
    """The file stacks.csv contains a list of lipid stack and contaminant stack mass
    differences. A stack is a series of ions differing in m/z by a user defined fixed mass
    multiple. Lipid stacks elute at same RT with the RT of contaminant stacks increasing in
    RT as overall m/z increases. Firstly, the m/z is checked for a lipid stack - if a stack
    is present all ions except the parent are deleted and the next m/z parent is checked.
    If no lipid stack is found then the m/z is checked for contaminant stacks, if found
    the whole stack incliding the parent is removed.

    Args:
        lfd (LipidFinderData): A LipidFinderData object that has been previously clustered
        and peaks found or is of XCMS origin and has been processed by xcmsReformat.py

    """
    stackFrame = pd.read_table('stacks.csv', sep=',')

    lipidStacks = stackFrame[stackFrame['Category'] == 'Lipid']
    lMZ = lipidStacks['MZ'].values
    contaminantStacks = stackFrame[stackFrame['Category'] == 'Contaminant']
    cMZ = contaminantStacks['MZ'].values

    # Column index of first sample first replicate
    indFirst = lfd.firstRepOffset

    # Column index of last sample last replicate
    indLast = lfd.lastRepOffset

    # Make a copy of lfd.processedDataFrame
    tempDataFrame = lfd.processedDataFrame.copy()

    mz = tempDataFrame['MZ'].values
    rt = tempDataFrame['Time'].values

    # Select out just the replicates only
    reps = tempDataFrame.ix[:, indFirst:indLast]

    # Add dummy column as first series apply possible run twice
    tempCol = tempDataFrame.ix[:, indFirst]
    reps.insert(0, 'Temp', tempCol)

    reps = reps.apply(__singleRepStackRemoval__, lMZ=lMZ, cMZ=cMZ, mz=mz, rt=rt, maxStackGap=lfd.maxStackGap, lipidStackAddition=lfd.lipidStackAddition, lowerMZLimit=lfd.lowerMZLimit,
                      upperMZLimit=lfd.upperMZLimit, mzRange=lfd.mzRange, lowerRTLimit=lfd.lowerRTLimit, upperRTLimit=lfd.upperRTLimit, rtTolMultipler=lfd.rtTolMultipler)

    # Drop dummy column
    reps.drop('Temp', inplace=True, axis=1)

    # Overwrite data in original processedDataFrame
    lfd.processedDataFrame.ix[:, indFirst:indLast] = reps

    lfd.stacksRemoved = lfd.processedDataFrame.copy()


def __singleRepStackRemoval__(replicate, lMZ, cMZ, mz, rt, maxStackGap, lipidStackAddition, lowerMZLimit, upperMZLimit, mzRange, lowerRTLimit, upperRTLimit, rtTolMultipler):

    # Get an index of all intensities that are not zero
    nonZeroIndex = replicate.values.nonzero()[0]

    # Create an array of intensities that are not zero (using the nonZeroIndex)
    nzInty = np.copy(replicate.values[nonZeroIndex])

    # Create an array of corresponding mz values where intensity is is not
    # zero (using the nonZeroIndex)
    nzMZ = np.copy(mz[nonZeroIndex])

    # Create an array of corresponding retention time values where intensity
    # is is not zero (using the nonZeroIndex)
    nzRT = np.copy(rt[nonZeroIndex])

    # Using a while loop instead of a for loop as the number of iterations
    # will alter as the loops progress because of deletions from tempFrame as
    # adducts are found. No need to check the last frame (hence the -1)
    parentInd = 0
    while parentInd < (nonZeroIndex.size - 1):
        parentFrameMass = nzMZ[parentInd]
        parentFrameRT = nzRT[parentInd]

        # Shows a stack has been found, serves dual purpose - If lipid stack
        # found, stops pointless continuation search further lipid stack and
        # contamination search, if found in contamination stack search stops
        # further search through possible contaminants - It can only be one
        # type of contaminant
        stackFlag = False

        # Lipid stack removal where mz and RT has to be identical (within tolerance)
        # For else construct used here, if for loop gets to end without a break
        # owing to a stack being found then code within else is executed
        for stackMass in lMZ:
            stackMassMult = stackMass
            maxGapCount = 0
            while maxGapCount <= maxStackGap:
                curStackMass = parentFrameMass + stackMassMult

                stackIndices = list(np.where(((nzMZ >= lowerMZLimit(curStackMass)) & (nzMZ <= upperMZLimit(curStackMass))) & (
                    (nzRT >= lowerRTLimit(parentFrameRT, rtTolMultipler)) & (nzRT <= upperRTLimit(parentFrameRT, rtTolMultipler))))[0])

                if not stackIndices:
                    stackMassMult += stackMass
                    maxGapCount += 1
                else:
                    stackInd = stackIndices[0]

                    if lipidStackAddition:
                        # Add intensity of stack index to be set to 0 to parent
                        # replicate
                        replicate.values[nonZeroIndex[
                            parentInd]] += nzInty[stackInd]

                    # Set replicate value to be 0
                    replicate.values[nonZeroIndex[stackInd]]

                    # Remove the new intensity 0 index from each array
                    nonZeroIndex = np.delete(nonZeroIndex, stackInd)
                    nzInty = np.delete(nzInty, stackInd)
                    nzMZ = np.delete(nzMZ, stackInd)
                    nzRT = np.delete(nzRT, stackInd)

                    stackMassMult += stackMass
                    stackFlag = True

            if stackFlag:
                break

        else:
            # Stack removal where closest RT above last mz stack multiple hit
            # RT
            for stackMass in cMZ:
                stackMassMult = stackMass
                maxGapCount = 0
                lastHitRT = parentFrameRT

                while maxGapCount <= maxStackGap:

                    curStackMass = parentFrameMass + stackMassMult
                    stackIndList = list(np.where(((nzMZ >= lowerMZLimit(curStackMass)) & (
                        nzMZ <= upperMZLimit(curStackMass))) & (nzRT > lowerRTLimit(lastHitRT, rtTolMultipler)))[0])
                    if not stackIndList:
                        maxGapCount += 1
                    else:
                        lowRTInd = nzRT[stackIndList].argmin()
                        lastHitRT = nzRT[stackIndList[lowRTInd]]
                        lastHitRTInd = stackIndList[lowRTInd]

                        # Set replicate value to be 0
                        replicate.values[nonZeroIndex[lastHitRTInd]] = 0

                        # Remove the new intensity 0 index from each array
                        nonZeroIndex = np.delete(nonZeroIndex, lastHitRTInd)
                        nzInty = np.delete(nzInty, lastHitRTInd)
                        nzMZ = np.delete(nzMZ, lastHitRTInd)
                        nzRT = np.delete(nzRT, lastHitRTInd)

                        stackFlag = True
                        maxGapCount = 0

                    stackMassMult += stackMass

                if stackFlag:

                    # Set replicate value to be 0
                    replicate.values[nonZeroIndex[parentInd]] = 0

                    # Remove the new intensity 0 index from each array
                    nonZeroIndex = np.delete(nonZeroIndex, parentInd)
                    nzInty = np.delete(nzInty, parentInd)
                    nzMZ = np.delete(nzMZ, parentInd)
                    nzRT = np.delete(nzRT, parentInd)

                    parentInd -= 1
                    break  # no point in looking further we have already removed a contaminant stack - cannot be 2 things at once therefore move to next mass
        parentInd += 1
    return replicate
