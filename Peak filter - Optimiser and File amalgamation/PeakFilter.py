from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import input
from past.utils import old_div
import pandas as pd
import os
import numpy as np

import LipidFinderData
import qcCalcs as qc
import solventCalcs as sol
import clustering as clus
import peakFinder as pf
import contaminantRemoval as cr
import rtCorrect as rtc
import outlierCorrect as oc
import sampleMeansCalc as smc
import reassignMass as rm
import broadContaminant as bc
import summary as sum
import timer as t
import time
import pickle
import copy as cpy


def main():
    """A process that takes a raw SIEVE or XCMS file and systematically takes the file
    through 16 stages of clean up to remove contaminants and redundant artifacts. The process
    begins by instantiating a LipidFinderData object. Various prompts are given to the user
    to confirm parameters drawn from the 'parameters.csv' file; the location, name and
    quantity of source files. The source file(s) are combined and stored within the
    LipidFinderData object as the current 'processedDataFrame' process runs through the
    following 16 stages saving a csv file after each stage providing a full audit trail.
    Output files are given a meaningful name according to their stage name combined with
    the time and date.

	Step 1 - Raw data importation.
    Step 2 - QC Sample Calculations and Reporting.
    Step 3 - Solvent removal.
    Step 4 - Low intensity removal.
    Step 5 - m/z clustering.
    Step 6 - Feature set clustering.
    Step 7 - Feature peak analysis.
    Step 8 - Mass contaminant removal.
    Step 9 - Adduct ion removal.
    Step 10 - Stack removal.
    Step 11 - Replicate retention time correction.
    Step 12 - Outlier correction.
    Step 13 - Sample mean calculation.
    Step 14 - Mean retention time correction.
    Step 15 - Mass reassignment.
    Step 16 - Broad retention time contaminant removal.

    """

    lfd = LipidFinderData.LipidFinderData()

    file_origin = lfd.checkSourceFileOrigin()

    lfd.checkDataParameters()

    lfd.setSourcePath()

    lfd.setDestFolder()

    lfd.setFinalFileName()

    lfd.setNumberOfFiles()

    lfd.importDataFrame(file_origin)

    lfd.fullRawData.to_csv(lfd.destPath + os.sep + '1-Combined_raw_data' +
                           lfd.fileSuffix + '.csv', sep=',', index=False)

    # Create a custom timer object, this will be the start time
    timer = t.Timer()

    # QC Sample Calculations and Reporting
    # Perform mean and RSD on QC samples if required
    # Check for 0 QC samples, in which case no QC samples present so step not
    # required
    if lfd.numberOfQCReps > 0:

        qcPropRSDLowerToUpper = qc.performQCCalcs(lfd)

        print("\n\n    *** Percentage %d%% QC-RSD samples to %d%% QC-RSD samples = %.1f%% ***" %
              (lfd.QCLowRSD, lfd.QCHighRSD, qcPropRSDLowerToUpper))

        lfd.qcData.to_csv(lfd.destPath + os.sep + '2-QC_data' +
                          lfd.fileSuffix + '.csv', sep=',', index=False)

        timer.mark()

        print("\n\n    *** QC Calculations completed in %.3f secs ***" %
              timer.individual())

    else:
        print("\n\n    *** No QC samples in data ***")

    # Solvent Calculations
    # Perform mean and RSD on solvent samples if present
    # Outlier correct
    # Delete frames where all replicates for each sample are not x * > solvent mean
    # Remove solvent mean intensity from remaining samples replicates' intensities
    # Check for 0 solvent samples, in which case no solvent samples present so
    # step not required
    if lfd.numberOfSolventReps > 0 and lfd.removeSolvent:

        sol.performSolCalcs(lfd)

        lfd.solventData.to_csv(lfd.destPath + os.sep + '3-Solvent_Data' +
                               lfd.fileSuffix + '.csv', sep=',', index=False)

        timer.mark()

        print("\n\n    *** Solvent calculations and processing completed in %.3f secs ***" %
              timer.individual())

    else:
        print("\n\n    *** No Solvent samples in data ***")

    sol.removeLowIntesityFrames(lfd)

    lfd.lowIntData.to_csv(lfd.destPath + os.sep + '4-Low_intensities_removed' +
                          lfd.fileSuffix + '.csv', sep=',', index=False)

    timer.mark()
    print("\n\n    *** Intensity cut off completed in %.3f secs ***" %
          timer.individual())

    # Only run if SIEVE input file
    if file_origin == 'S':

        clus.mzCluster(lfd)

        # output file after mz clustering but before finding feature clusters
        lfd.mzClusteredData.to_csv(
            lfd.destPath + os.sep + '5-Mass_clustered_data' + lfd.fileSuffix + '.csv', sep=',', index=False)

        timer.mark()
        print("\n\n    *** mz clustering completed in %.3f secs ***" %
              timer.individual())

        clus.featureCluster(lfd)

        # output file after feature clustering but before peak finding
        lfd.featureClusteredData.to_csv(
            lfd.destPath + os.sep + '6-Feature_clustered_data' + lfd.fileSuffix + '.csv', sep=',', index=False)

        timer.mark()
        print("\n\n    *** Feature clustering completed in %.3f secs ***" %
              timer.individual())

        pf.processAllFeatures(lfd)

        # output file after peak finding
        lfd.peakFound.to_csv(lfd.destPath + os.sep + '7-Peak_found_data' +
                             lfd.fileSuffix + '.csv', sep=',', index=False)

        timer.mark()
        print("\n\n    *** Peak finding completed in %.3f secs ***" %
              timer.individual())

    elif file_origin == 'X':

        # XCMS only
        # Create 'Feature Cluster ID' and 'mz Cluster ID' columns - each row is
        # already a feature.
        lfd_copy = cpy.deepcopy(lfd)

        clus.mzCluster(lfd_copy)
        timer.mark()

        clus.featureCluster(lfd_copy)
        timer.mark()

        lfd.processedDataFrame.sort_values(by=['MZ'], inplace=True)
        lfd.processedDataFrame.reset_index(inplace=True, drop=True)

        copy_mz_cluster_id = lfd_copy.processedDataFrame['mz Cluster ID']
        copy_feature_cluster_id = lfd_copy.processedDataFrame[
            'Feature Cluster ID']

        lfd.processedDataFrame.insert(
            lfd.firstRepOffset, 'mz Cluster ID', copy_mz_cluster_id)
        lfd.incFirstRepOffset()
        lfd.processedDataFrame.insert(
            lfd.firstRepOffset, 'Feature Cluster ID', copy_feature_cluster_id)
        lfd.incFirstRepOffset()

    if lfd.removeContaminants:
        cr.contaminantsRemoval(lfd)

        # output file after contaminant removal
        lfd.contaminantsRemoved.to_csv(
            lfd.destPath + os.sep + '8-Mass_contaminant_removed_data' + lfd.fileSuffix + '.csv', sep=',', index=False)

        timer.mark()
        print("\n\n    *** Mass contaminant removal completed in %.3f secs ***" %
              timer.individual())

    if lfd.removeAdducts:
        cr.adductsRemoval(lfd)

        # output file after contaminant removal
        lfd.adductsRemoved.to_csv(
            lfd.destPath + os.sep + '9-Adducts_removed_data' + lfd.fileSuffix + '.csv', sep=',', index=False)

        timer.mark()
        print("\n\n    *** Adduct removal completed in %.3f secs ***" %
              timer.individual())

    if lfd.removeStacks:
        cr.stacksRemoval(lfd)

        # output file after contaminant removal
        lfd.stacksRemoved.to_csv(
            lfd.destPath + os.sep + '10-Stacks_removed_data' + lfd.fileSuffix + '.csv', sep=',', index=False)

        timer.mark()
        print("\n\n    *** Stack removal completed in %.3f secs ***" %
              timer.individual())

    # RT correct each set of sample replicates to allow for Sieve alignment
    # errors if replicates > 1
    if lfd.numberOfTechReps > 1:
        lfd.processedDataFrame = rtc.process(lfd)
        lfd.rtCorrected = lfd.processedDataFrame.copy()

        # output file after rt correction
        lfd.rtCorrected.to_csv(lfd.destPath + os.sep + '11-Replicate_Retention_time_corrected_data' +
                               lfd.fileSuffix + '.csv', sep=',', index=False)

        timer.mark()
        print("\n\n    *** Replicate retention time correction completed in %.3f secs ***" %
              timer.individual())

    # Outlier correct sample replicates
    lfd.processedDataFrame = oc.process(
        lfd, lfd.firstRepOffset, lfd.lastRepOffset, lfd.numberOfTechReps)
    lfd.outlierCorrected = lfd.processedDataFrame.copy()

    # output file after outlier correction
    lfd.outlierCorrected.to_csv(
        lfd.destPath + os.sep + '12-Outlier_corrected_data' + lfd.fileSuffix + '.csv', sep=',', index=False)

    timer.mark()
    print("\n\n    *** Replicate outlier correction completed in %.3f secs ***" %
          timer.individual())

    smc.processSampleMeans(lfd)

    # output file after contaminant removal
    lfd.MeanedData.to_csv(lfd.destPath + os.sep + '13-Meaned_data' +
                          lfd.fileSuffix + '.csv', sep=',', index=False)

    timer.mark()
    print("\n\n    *** Adding sample means completed in %.3f secs ***" %
          timer.individual())

    if lfd.rtCorrectMeans:
        # RT correct the means to each other to allow for Sieve alignment
        # errors
        lfd.processedDataFrame = rtc.process(lfd, True)
        lfd.rtMeansCorrected = lfd.processedDataFrame.copy()

        # output file after mean rt correction
        lfd.rtMeansCorrected.to_csv(
            lfd.destPath + os.sep + '14-Mean_Retention_time_corrected_data' + lfd.fileSuffix + '.csv', sep=',', index=False)

    rm.reassignFrameMasses(lfd)

    # output file after mass reassignment
    lfd.reassignedFrameMasses.to_csv(
        lfd.destPath + os.sep + '15-Mass_reassigned_data' + lfd.fileSuffix + '.csv', sep=',', index=False)

    timer.mark()
    print("\n\n    *** Reassigning masses completed in %.3f secs ***" %
          timer.individual())

    bc.processAllFeatures(lfd)

    # output file after contaminant removal
    lfd.broadContaminantCorrected.to_csv(
        lfd.destPath + os.sep + '16-Broad_contaminant_removed_data' + lfd.fileSuffix + '.csv', sep=',', index=False)
    timer.mark()
    print("\n\n    *** Broad contaminant removal completed in %.3f secs ***" %
          timer.individual())

    sum.createSummaryDataFrame(lfd)
    # output files after amalgamation processing
    # Leave index in summary, used by amalgamator to allow comparison between
    # files
    lfd.summaryDataFrame.to_csv(
        lfd.destPath + os.sep + lfd.finalFileName + lfd.fileSuffix + '_summary.csv', sep=',')
    lfd.processedDataFrame.to_csv(
        lfd.destPath + os.sep + lfd.finalFileName + lfd.fileSuffix + '_details.csv', sep=',', index=False)

    timer.mark()
    print("\n\n    *** Retention time cut-off and final output files created completed in %.3f secs ***" %
          timer.individual())

    __pickleData__(lfd)

    print("\n=================================================================================================================\n\n    *** LipidFinder completed in %.3f secs ***" % timer.total())


def __pickleData__(lfd):
    pickleSize = np.rint(old_div(len(pickle.dumps(lfd)), 1e6)).astype(int)
    while True:
        pickleFlag = input(
            "\n    Do you want to pickle the LipidFinderData object (%dMB) for future analysis? [N]: " % pickleSize) or 'N'
        if pickleFlag.upper() == 'N':
            break
        elif pickleFlag.upper() == 'Y':
            pickle.dump(lfd, open(lfd.destPath + '/allData' +
                                  lfd.fileSuffix + '.pkl', 'wb'))
            break
        else:
            print("\n    Y or N required!")

if __name__ == "__main__":
    main()
