from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import input
import pandas as pd
import os
import numpy as np
import time
import pickle
import sys

import LipidFinderData
import solventCalcs as sol
import clustering as clus
import peakFinder as pf
import outlierCorrect as oc
import timer as t
import HillClimber as HC


def main():
    """Prepare the inputs to and run the hill HillClimber class to optimise the parameter
    inputs for PeakFilter.

    """
    lfd = LipidFinderData.LipidFinderData()

    lfd.checkDataParameters()

    lfd.setSourcePath()

    lfd.setNumberOfFiles()

    lfd.importDataFrame('S')

    # Get target file
    tDF = __getTargetFile__(lfd.sourcePath)

    # Get target replicate column name (the replicate the target file was
    # built from)
    tRep = __getTargetReplicate__(lfd)

    # Prompt for number of HC runs as random restarts to help avoid the best
    # parameter set being sub optimal
    hc_run_limit = __get_HCRun_limit__()

    # Prompt for maximum 1opt cycles to complete (to avoid running for too
    # long)
    max_cycle_1opt_count = __get_max_cycle_1opt_count__()

    # Prompt for whether to display progress on screen
    debug = __get_debug__()

    # Create a custom timer object, this will be the start time
    timer = t.Timer()

    # Solvent Calculations
    # Perform mean and RSD on solvent samples if present
    # Gap fill for bad data
    # Delete frames where all replicates for each sample are not 3x > solvent mean
    # Remove solvent mean intensity from remaining samples replicates' intensities
    # Check for 0 solvent samples, in which case no solvent samples present so
    # step not required
    if lfd.numberOfSolventReps > 0 and lfd.removeSolvent:

        sol.performSolCalcs(lfd, False)

        sol.removeLowIntesityFrames(lfd, False)

        timer.mark()

        print("\n\n***Solvent outlier correction and mean calculations completed in %.3f secs ***" %
              timer.individual())

    else:
        print("\n\n*** No Solvent samples in data ***")

    timer.mark()

    # Create hill climb object
    hillClimb = HC.HillClimber(
        lfd, tDF, tRep, timer, hc_run_limit, max_cycle_1opt_count, debug)

    # Perform hillclimbing with multiple restarts
    hillClimb.hc_mr()


def __getTargetFile__(sourcePath):
    while True:
        targetFileName = sourcePath + os.sep + \
            input(
                "\nEnter target file name:-\n\n    ")  # 'HC_Target_cd1_NPP.csv'
        try:
            # Read file into DataFrame and return
            return pd.read_table(targetFileName, sep=',')
        except IOError as e:
            print("\n    " + str(e) + "!")
            # Prompt user again for same file number


def __getTargetReplicate__(lfd):
    while True:
        targetRep = input(
            "\nWhat is the replicate name used to generate the target file:- ")  # 'II_Cd1'

        if targetRep in lfd.fullRawData.columns.values:
            return targetRep
        else:
            print("\n    Cannot find column \'%s\' in the source file(s)!" % targetRep)


def __get_HCRun_limit__():
    rl = 1
    while True:
        try:
            rl = int(input(
                "\nHow many hill climbing iterations to perform? [%s]:- " % rl) or rl)
            if rl > 0:
                return rl
            else:
                print("   Integer > 0 needed!")
        except ValueError:
            print("   Integer > 0 needed!")


def __get_max_cycle_1opt_count__():
    max_cycle_1opt_count = 5
    while True:
        try:
            max_cycle_1opt_count = int(input(
                "\nWhat is the maximum move set iterations allowed? [%s]:- " % max_cycle_1opt_count) or max_cycle_1opt_count)
            if max_cycle_1opt_count > 0:
                return max_cycle_1opt_count
            else:
                print("   Integer > 0 needed!")
        except ValueError:
            print("   Integer > 0 needed!")


def __get_debug__():
    debug = 'Y'
    while True:
        debug_raw = input(
            "\nPrint progress on screen? use 'Y' for yes and 'N' for no) [%s]: " % debug) or debug
        upper_debug_raw = debug_raw.upper()
        if upper_debug_raw in ['Y', 'N']:
            if upper_debug_raw == 'Y':
                return True
            else:
                return False
        else:
            print("\n'Y' or 'N' required!")


if __name__ == "__main__":
    main()
