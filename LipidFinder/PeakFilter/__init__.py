# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Filter contaminants and redundant artifacts from a LC/MS data
pre-processed by XCMS or another pre-processing tool.

The peak filter process is composed by 19 stages (some of them optional)
designed to remove contaminants and redundant artifacts from the input
dataset. Every frame filtered and calculation computed (e.g. False
Discovery Rate) is recorded in a log file by the corresponding module,
so the user can track what is happening with the data at every step of
the process.

Step  1 - QC Samples Calculations and Reporting.
Step  2 - Solvent removal.
Step  3 - Low intensity removal.
Step  4 - Mass clustering.
Step  5 - Feature set clustering.
Step  6 - Feature peak analysis.
Step  7 - In-source ion fragmentation removal.
Step  8 - Mass contaminant removal.
Step  9 - Adduct ion removal.
Step 10 - Stack removal.
Step 11 - Replicate retention time correction.
Step 12 - Outlier correction.
Step 13 - Sample mean calculation.
Step 14 - Mean retention time correction.
Step 15 - Mass reassignment.
Step 16 - Broad retention time contaminant removal.
Step 17 - Isotope removal.
Step 18 - Salt cluster removal.
Step 19 - False Discovery Rate.

Finally, a summary file (CSV format) is created with the most relevant
information of the data after being processed: id, m/z, retention time,
polarity and samples mean. Additionally, the complete filtered data can
be saved in a CSV file too.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> import PeakFilter
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('xcms_output.csv', parameters)
    >>> PeakFilter.peak_filter(data, parameters)
"""

import logging
import os
import warnings

import pandas

from LipidFinder.PeakFilter import BroadContaminant
from LipidFinder.PeakFilter import Clustering
from LipidFinder.PeakFilter import ContaminantRemoval
from LipidFinder.PeakFilter import FalseDiscoveryRate
from LipidFinder.PeakFilter import InSrcFragRemoval
from LipidFinder.PeakFilter import Deisotoping
from LipidFinder.PeakFilter import MassDefectFilter
from LipidFinder.PeakFilter import MassReassignment
from LipidFinder.PeakFilter import OutlierCorrection
from LipidFinder.PeakFilter import PeakFinder
from LipidFinder.PeakFilter import QCCalcs
from LipidFinder.PeakFilter import RTCorrection
from LipidFinder.PeakFilter import SampleMeansCalc
from LipidFinder.PeakFilter import SolventCalcs
from LipidFinder.PeakFilter import Summary
from LipidFinder._utils import print_progress_bar


# Ignore Future Warnings from pandas library
warnings.simplefilter(action='ignore', category=FutureWarning)
# Progress bar increment per step
INCREMENT = 100.0 / 18

def _update_status(data, stepDst, verbose, stepNum):
    # type: (LFDataFrame, str, bool, int) -> None
    """Create CSV file from 'data' in 'stepDst', update progress bar and
    return incremented step number.

    Keyword Arguments:
        data    -- LFDataFrame instance
        stepDst -- destination directory for CSV 'data' file
        verbose -- create CSV 'data' file?
        stepNum -- step number completed
    """
    # Update progress bar
    print_progress_bar(INCREMENT * stepNum, 100, prefix='PeakFilter progress:')
    if (verbose):
        # Create a CSV file with the whole processed dataframe
        outFileName = 'peakfilter_step_{:02d}.csv'.format(stepNum)
        data.to_csv(os.path.join(stepDst, outFileName), index=False)
    stepNum += 1
    return stepNum


def peak_filter(data, parameters, dst='', verbose=False):
    # type: (LFDataFrame, LFParameters, str, bool) -> None
    """Filter contaminants and redundant artifacts from a LC/MS data
    pre-processed by XCMS or another pre-processing tool.

    If 'dst' is not an absolute path, the current working directory will
    be used as starting point. If either "peakfilter_<polarity>.csv" or
    "peakfilter_<polarity>_summary.csv" files already exist, they will
    be overwritten. "<polarity>" stands for "positive" or "negative", as
    stated in the parameters.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
        dst        -- destination directory where the log file, the
                      processed data CSV file and the summary CSV file
                      will be saved [default: current working directory]
        verbose    -- create folder inside 'dst' where the intermediate
                      results will be saved in CSV files
    """
    # Start progress bar
    print_progress_bar(0, 100, prefix='PeakFilter progress:')
    # Set the log file where the information about the steps performed
    # is saved
    logFilePath = 'peakfilter.log'
    if (dst):
        logFilePath = os.path.join(dst, logFilePath)
    # Create logger and its file handler
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(logFilePath)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(asctime)s] %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    # Write initial information in log file
    logger.info('Starting PeakFilter. Input dataframe ("%s") has %d rows.',
                 data.src, len(data.index))
    # Prepare the folder structure to store the intermediate files
    stepDst = os.path.join(dst, 'step_by_step')
    if (verbose and not os.path.isdir(stepDst)):
        os.makedirs(stepDst)
    stepNum = 1
    # QC Sample Calculations
    if (parameters['numQCReps'] > 0):
        # Perform mean and RSD on QC samples
        qcRatio = QCCalcs.qc_rsd_ratio(data, parameters)
        # Write report in log file
        logger.info(("QC Sample Calculations completed. %.1f%% samples between"
                      " %d%% and %d%% QC-RSD"), qcRatio, parameters["QCRSD"][0],
                     parameters["QCRSD"][1])
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Solvent Calculations
    if ((parameters['numSolventReps'] > 0) and parameters['removeSolvents']):
        # Perform mean and RSD on solvent samples, perform the outlier
        # correction, remove frames where all technical replicates of
        # all samples are less than the 'solventMinFoldDiff' times the
        # solvent mean, and remove solvent mean intensity from remaining
        # intensities of samples replicates
        SolventCalcs.remove_solvent_effect(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Background correction: remove low intensity frames
    SolventCalcs.remove_low_intensity_frames(data, parameters)
    if (parameters['preprocSoftware'] == 'XCMS'):
        # Get m/z clusters required by 'MassReassignment' and
        # 'BroadContaminant' modules
        Clustering.cluster_by_mz(data, parameters)
        # Create the "FeatureClusterID" column that will be used by
        # 'RTCorrection' step. In XCMS, each row is already a feature.
        data['FeatureClusterID'] = range(1, len(data) + 1)
    else:
        # Perform peak finding for any other pre-processing software
        PeakFinder.process_features(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # In-source ion fragment removal
    if (parameters['removeIonFrags']):
        InSrcFragRemoval.remove_in_src_frags(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Contaminant removal
    if (parameters['removeContaminants']):
        ContaminantRemoval.remove_contaminants(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Adduct removal
    if (parameters['removeAdducts']):
        ContaminantRemoval.remove_adducts(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Stack removal
    if (parameters['removeStacks']):
        ContaminantRemoval.remove_stacks(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Retention time correction of each set of sample replicates to fix
    # other pre-processing tool's likely alignment errors
    if ((parameters['numTechReps'] > 1)
        and (parameters['preprocSoftware'] == 'Other')):
        RTCorrection.correct_retention_time(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Remove outliers from sample replicates
    OutlierCorrection.remove_outliers(data, parameters, src='samples')
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Calculate and add the mean of each sample's replicates
    SampleMeansCalc.calculate_sample_means(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Retention time correction to the means of the sample replicates
    if (parameters['correctRTMeans']):
        RTCorrection.correct_retention_time(data, parameters, True)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Assign each m/z in either a mass or feature cluster to the m/z of
    # the row containing the highest sample mean intensity
    MassReassignment.reassign_frame_masses(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Remove ions with similar intensities for the same m/z that are
    # likely to be contaminants
    BroadContaminant.process_all_features(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Isotope removal
    Deisotoping.remove_isotopes(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Mass defect filter: remove salt clusters
    if (parameters['filterMassDefect']):
        MassDefectFilter.remove_salt_clusters(data, parameters)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Calculate the False Discovery Rate
    if (parameters['calculateFDR']):
        try:
            fdrValue = FalseDiscoveryRate.get_fdr(data, parameters)
            message = ("False Discovery Rate for selected data and parameters: "
                       "{0:.2%}").format(fdrValue)
        except ValueError as e:
            message = 'ValueError: ' + e.args[0]
        except Exception as oe:
            message = 'OtherError: ' + oe.args[0]
        logger.info(message)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Create summary CSV file from the processed dataframe
    Summary.create_summary(data, parameters, dst)
    stepNum = _update_status(data, stepDst, verbose, stepNum)
    # Create a CSV file with the whole processed dataframe
    data['Polarity'] = parameters['polarity']
    outFileName = 'peakfilter_{0}.csv'.format(parameters['polarity'].lower())
    data.to_csv(os.path.join(dst, outFileName), index=False)
    # Update progress bar
    print_progress_bar(100, 100, prefix='PeakFilter progress:')
    # Print False Discovery Rate message
    if (parameters['calculateFDR']):
        print(message)
    # Write the final information in log file and close handler
    logger.info('PeakFilter completed. Output dataframe has %d rows.\n',
                 len(data.index))
    handler.close()
    logger.removeHandler(handler)
