from __future__ import print_function
from __future__ import division
from builtins import input
from builtins import str
from past.utils import old_div
import time
import pandas as pd
import numpy as np
from math import sqrt
import os
import LipidFinderData
import timer as t
import sys

# Proton mass
hmass = 1.007276

# Electron mass = 0.00054858 - http://fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator/
# Carbon 12 mass = 12 exactly

# [CH3]-
# 1 * C12 + 3 * proton + 2 * electron {12 + (3 * 1.007276) + (2 * 0.00054858)}
methmass = 15.022925

# Global
lfd = None


def main():
    """The summary positive and negative file for the experiment are combined with lower
    intensity redundant matching ions present in both files removed. The columns headings
    in each file are checked for matches.

    Two output files are created, a results file (the input to WebSearch) and a matches
    file to allow match checking.
    """
    global lfd
    lfd = LipidFinderData.LipidFinderData()

    lfd.setNumSamples()

    lfd.setSourcePath(3)

    lfd.setDestFolder(3)

    lfd.setFinalFileName(3)

    negData = __setData__("Negative file name:-")
    posData = __setData__("Positive file name:-")

    # Combine intensities flag
    combIntFlag = __combInt__()

    # Create a custom timer object, this will be the start time
    timer = t.Timer()

    # Lambda  Round up to next integer, but leave as float (rint), convert nan to 0 (nan_to_num) (some indices will produce nan when mean is taken as there will be 0 nonzeros)
    # get the indices of nonzero values in the array (nonzero) and slice the
    # array to leave only nonzero values (x(nonzero)), take the mean of these
    # non zero values (mean) and then convert to integer
    negData['Total Mean'] = negData.ix[:, 4:].apply(lambda x: np.rint(
        np.nan_to_num(x[np.nonzero(x)[0]].mean())).astype(int), axis=1)
    posData['Total Mean'] = posData.ix[:, 4:].apply(lambda x: np.rint(
        np.nan_to_num(x[np.nonzero(x)[0]].mean())).astype(int), axis=1)

    # Get column names to ensure files have same samples
    negCols = list(negData.columns.values)
    posCols = list(posData.columns.values)

    # To allow the same column names in a different order as only taking
    # mean of all samples
    negCols.sort()
    posCols.sort()

    # check column names match in pos and neg file
    if not negCols == posCols or not set(['MZ', 'Time']).issubset(negCols):

        sys.exit(
            'Column names in the pos and neg files do not match or there is no \'MZ\' or no \'Time\' column')

    negData.reset_index(inplace=True)
    posData.reset_index(inplace=True)

    nind = negData.index.values
    nmz = negData['MZ'].values
    nrt = negData['Time'].values
    nmeans = negData['Total Mean'].values

    negCol = list(negData.columns.values)
    posCol = list(posData.columns.values)

    # Empty results dataframe
    results = pd.DataFrame(columns=negCol)

    # Empty matches dataframe
    matchesOut = pd.DataFrame(columns=['match_Index'] + negCol)

    # Mass differences
    H2 = 2 * hmass
    HCH3 = hmass + methmass

    # Loop through indices in negative file
    for i in nind:
        negMass = nmz[i]
        negRT = nrt[i]
        pmz = posData['MZ'].values
        prt = posData['Time'].values
        pmeans = posData['Total Mean'].values

        negMassH2 = negMass + H2
        negMassHCH3 = negMass + HCH3

        matchesH2 = list(np.where(((pmz >= lfd.lowerMZLimit(negMassH2)) & (pmz <= lfd.upperMZLimit(
            negMassH2))) & ((prt >= lfd.lowerRTLimit(negRT)) & (prt <= lfd.upperRTLimit(negRT))))[0])

        # Look for 2*H+ match first, if none look for [H+] + [CH3+ +2e-] match
        if matchesH2:
            indMatch = __bestMatch__(matchesH2, negMassH2, pmz, negRT, prt)
            if pmeans[indMatch] > nmeans[i]:

                results = __addPosFrame__(results, posData, indMatch)

                if combIntFlag:
                    # Last record in results has lower matched ion intensity
                    # added

                    results.iloc[-1, 5:] = __calcRowUpdate__(
                        results.iloc[-1, 5:], negData.ix[i, 5:])

                    results.iloc[-1, 3] = 'POSC'

                matchesOut = __addMatches__(
                    matchesOut, negData, posData, i, indMatch, negCol, posCol)

                # Remove pos match from pos file as it's been considered
                posData.drop(indMatch, inplace=True)
                posData.reset_index(inplace=True, drop=True)
                pmz = posData['MZ'].values
                prt = posData['Time'].values
                pmeans = posData['Total Mean'].values
                continue

            else:

                results = __addNegFrame__(results, negData, i)

                if combIntFlag:
                    # Last record in results has lower matched ion intensity
                    # added

                    results.iloc[-1, 5:] = __calcRowUpdate__(
                        results.iloc[-1, 5:], posData.ix[indMatch, 5:])

                    results.iloc[-1, 3] = 'NEGC'

                matchesOut = __addMatches__(
                    matchesOut, negData, posData, i, indMatch, negCol, posCol)

                # Remove pos match from pos file as it's been considered
                posData.drop(indMatch, inplace=True)
                posData.reset_index(inplace=True, drop=True)
                pmz = posData['MZ'].values
                prt = posData['Time'].values
                pmeans = posData['Total Mean'].values
                continue

        matchesHCH3 = list(np.where(((pmz >= lfd.lowerMZLimit(negMassHCH3)) & (pmz <= lfd.upperMZLimit(
            negMassHCH3))) & ((prt >= lfd.lowerRTLimit(negRT)) & (prt <= lfd.upperRTLimit(negRT))))[0])

        if matchesHCH3:
            indMatch = __bestMatch__(matchesHCH3, negMassHCH3, pmz, negRT, prt)
            if pmeans[indMatch] > nmeans[i]:

                results = __addPosFrame__(results, posData, indMatch)

                if combIntFlag:
                    # Last record in results has lower matched ion intensity
                    # added

                    results.iloc[-1, 5:] = __calcRowUpdate__(
                        results.iloc[-1, 5:], negData.ix[i, 5:])

                    results.iloc[-1, 3] = 'POSC'

                matchesOut = __addMatches__(
                    matchesOut, negData, posData, i, indMatch, negCol, posCol)

                # Remove pos match from pos file as it's been considered
                posData.drop(indMatch, inplace=True)
                posData.reset_index(inplace=True, drop=True)
                pmz = posData['MZ'].values
                prt = posData['Time'].values
                pmeans = posData['Total Mean'].values
                continue
            else:

                results = __addNegFrame__(results, negData, i)

                if combIntFlag:
                    # Last record in results has lower matched ion intensity
                    # added

                    results.iloc[-1, 5:] = __calcRowUpdate__(
                        results.iloc[-1, 5:], posData.ix[indMatch, 5:])

                    results.iloc[-1, 3] = 'NEGC'

                matchesOut = __addMatches__(
                    matchesOut, negData, posData, i, indMatch, negCol, posCol)

                # Remove pos match from pos file as it's been considered
                posData.drop(indMatch, inplace=True)
                posData.reset_index(inplace=True, drop=True)
                pmz = posData['MZ'].values
                prt = posData['Time'].values
                pmeans = posData['Total Mean'].values
                continue

        results = __addNegFrame__(results, negData, i)

    # Append whatever remains in positive file as there were no matches in
    # negative file
    results = results.append(posData, ignore_index=True)

    results.drop(['index', 'Total Mean'], axis=1, inplace=True)

    results.to_csv(lfd.destPath + os.sep + lfd.finalFileName + '_results' +
                   lfd.timeStamp + '.csv', sep=',', float_format='%.6f', index=False)

    matchesOut = matchesOut.rename(columns={'index': 'original_file_index'})

    matchesOut.to_csv(lfd.destPath + os.sep + lfd.finalFileName + '_matches' +
                      lfd.timeStamp + '.csv', sep=',', float_format='%.6f', index=False)

    timer.mark()
    print("\n\n    *** Positive and negative files amalgamated in %.3f secs ***\n" % timer.total())

# Prompt user to create input file name for negative file and return dataframe


def __setData__(text=""):
    while True:
        try:
            file = lfd.sourcePath + os.sep + input(text + "\n\n    ")
            data = pd.read_table(file, sep=',', index_col=0)
            return data
        except IOError as e:
            print("\n    " + str(e) + "!")


def __combInt__():
    Flag = "N"
    while True:
        Flag = input(
            "Do you want to combine intensities for ions of the same molecule found in positive and negative mode? (use 'Y' for Yes and 'N' for No) [%s]: " % Flag) or Flag
        upperFlag = Flag.upper()
        if upperFlag == 'Y':
            return True
        elif upperFlag == 'N':
            return False
        else:
            print("\n   'Y' or 'N' required!")


def __hitScore__(sourceMass, targetMass, sourceRT, targetRT):
    mzRange = sourceMass - lfd.lowerMZLimit(sourceMass)
    mzDiff = abs(sourceMass - targetMass)
    rtRange = sourceRT - lfd.lowerRTLimit(sourceRT)
    rtDiff = abs(sourceRT - targetRT)
    return sqrt(min(old_div(mzDiff, mzRange), 1.0)**2.0 + min(old_div(rtDiff, rtRange), 1.0)**2)


def __bestMatch__(matches, negMass, pmz, negRT, prt):
    maxScore = 0.0
    maxScoreIndex = 0
    for ind in matches:
        curScore = __hitScore__(negMass, pmz[ind], negRT, prt[ind])
        if curScore > maxScore:
            maxScore = curScore
            maxScoreIndex = ind
    return maxScoreIndex


def __addPosFrame__(results, posData, posInd):
    # Pos frame as a series
    posFrame = posData.ix[posInd]

    # Append posFrame frame to results
    return results.append(posFrame, ignore_index=True)


def __addNegFrame__(results, negData, negInd):
    # Neg frame as a series
    negFrame = negData.ix[negInd]

    # Append negFrame frame to results
    return results.append(negFrame, ignore_index=True)


def __addMatches__(matchesOut, negData, posData, negInd, posInd, negCol, posCol):

    # Neg frame as a series
    negFrame = negData.ix[negInd]

    # Pos frame as a series
    posFrame = posData.ix[posInd]

    # Get current match number
    matchNo = old_div((2 + len(matchesOut)), 2)

    # Create negFrame series with matchNo at the beginning for matches file
    negFrameInd = negFrame.append(pd.Series(matchNo, ['match_Index'])).reindex(
        index=['match_Index'] + negCol)

    # Create posFrame series with negFrame index at the beginning for matches
    # file
    posFrameInd = posFrame.append(pd.Series(matchNo, ['match_Index'])).reindex(
        index=['match_Index'] + posCol)

    return matchesOut.append([negFrameInd, posFrameInd], ignore_index=True)

# Helper function to get update array when summing positive and negative ions


def __calcRowUpdate__(df1_row, df2_row):
    a = df1_row.values
    b = df2_row.values

    return a + b

if __name__ == "__main__":
    main()
