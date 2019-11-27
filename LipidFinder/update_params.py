#!/usr/bin/env python

# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Transform the old parameters CSV file for PeakFilter and Amalgamator
(LipidFinder v1.0) to the new parameters JSON files for the same
modules. Files "peakfilter.json" and "amalgamator.json" will be created
in the given folder or in the current working directory (default).

The resultant JSON files will be incomplete (some new parameters have
been introduced in the latest release) and will raise an error when used
as argument for their corresponding module. They must be used as
argument for "config_params.py" first to fill in the missing parameters
and generate a complete version.
"""

import argparse
import ast
import json
import os

import pandas

from LipidFinder._utils import normalise_path
from LipidFinder._py3k import viewitems


# New-to-old parameter name equivalence dict
PARAMS_DICT = {'firstSampleIndex': 'firstRepOffset',
               'numSamples': 'numberOfSamples',
               'numTechReps': 'numberOfTechReps',
               'numQCReps': 'numberOfQCReps',
               'numSolventReps': 'numberOfSolventReps',
               'polarity': 'filePolarityMode',
               'QCRSD': ['QCLowRSD', 'QCHighRSD'],
               'removeSolvents': 'removeSolvent',
               'solventMinFoldDiff': 'solventFoldCutOff',
               'intenSignifCutOff': 'intensitySignificanceCutOff',
               'mzFixedError': 'mzFixedError',
               'mzPPMError': 'mzSizeErrorPPM',
               'peakMaxRTWidth': 'peakMaxRTWidth',
               'peakMinFoldDiff': 'peakMinFoldCutOff',
               'maxRTDiffAdjFrame': 'peakAdjacentFrameMaxRT',
               'concatAllFrames': 'peakConcatenateAllFrames',
               'removeContaminants': 'removeContaminant',
               'removeAdducts': 'removeAdduct',
               'adductAddition': 'adductAddition',
               'removeStacks': 'removeStack',
               'maxStackGap': 'maxStackGap',
               'lipidStackAddition': 'lipidStackAddition',
               'intenOutlierCutOff': 'outlierHighIntensityValue',
               'intensityRSD': ['outlierLowIntensityRSD',
                                'outlierHighIntensityRSD'],
               'featMassAssignment': 'featureLevelMassAssignment',
               'negAdductsPairs': 'negativeModeAdductPairs',
               'posAdductsPairs': 'positiveModeAdductPairs',
               'outlierMinDiff': 'broadContsdMult',
               'minNonZeroPoints': 'broadContminPoints',
               'intenRSDCutOff': 'broadContRSDCutOff',
               'rtSDCutOff': 'broadContrtSDCutOff',
               'rtRange': ['retentionTimeLowCutOff',
                           'retentionTimeHighCutOff'],
               'intensityStDev': 'rtCorrectStDev',
               'correctRTMeans': 'rtCorrectMeans'
              }


def _adduct_rename(adducts, index):
    # type: (pandas.DataFrame, int) -> str
    """Return LipidFinder v2.0 name for the selected adduct: some
    adducts have been renamed to follow the standard nomenclature.

    Keyword parameters:
        adducts -- LipidFinder v1.0 adducts dataframe
        index   -- row index
    """
    name = adducts.iloc[index, 1]
    if (name == 'M+AcO-'):
        name = 'M+OAc'
    elif (name == 'M+Cl-'):
        name = 'M+Cl'
    return name


def main():
    # Create the argument parser and parse the arguments
    parser = argparse.ArgumentParser(
            description=(
                    "Convert LipidFinder's old parameters CSV file to the new "
                    "JSON format: 'peakfilter.json' and 'amalgamator.json' "
                    "files will be created."))
    parser.add_argument('-p', '--params', required=True, type=str,
                        help="LipidFinder v1.0 parameters CSV file")
    parser.add_argument('-a', '--adducts', metavar='FILE', required=True,
                        type=str, help="LipidFinder v1.0 adducts CSV file")
    parser.add_argument('-o', '--output', metavar='DIR', type=str,
                        help="folder where the parameter files will be saved")
    parser.add_argument('--version', action='version',
                        version="LipidFinder 2.0")
    args = parser.parse_args()
    # Check if the output directory exists. If not, create it.
    dst = args.output if (args.output) else '.'
    dst = normalise_path(dst)
    if (not os.path.isdir(dst)):
        os.makedirs(dst)
    # Read old parameters CSV file and generate a dict
    oldParams = pandas.read_csv(args.params, index_col='Parameter name',
                                usecols=['Parameter name', 'Current value'])
    oldParams = oldParams.to_dict()['Current value']
    # Use the name equivalence dict to create the new parameters dict
    pfParams = {}
    for key, value in viewitems(PARAMS_DICT):
        if (isinstance(value, list)):
            pfParams[key] = [ast.literal_eval(oldParams[x]) for x in value]
        else:
            # literal_eval() raises an exception when parsing strings
            try:
                pfParams[key] = ast.literal_eval(oldParams[value])
            except ValueError:
                if (oldParams[value] == 'N'):
                    pfParams[key] = 'Negative'
                elif (oldParams[value] == 'P'):
                    pfParams[key] = 'Positive'
                else:
                    # String containing a boolean value
                    pfParams[key] = oldParams[value] == 'TRUE'
    # Correct starting index in "firstSampleIndex" from 0 to 1
    pfParams['firstSampleIndex'] += 1
    # Ensure float type parameters have float values
    floatParams = ['solventMinFoldDiff', 'mzFixedError', 'mzPPMError',
                   'peakMaxRTWidth', 'peakMinFoldDiff',
                   'maxRTDiffAdjFrame', 'intensityStDev']
    for key in floatParams:
        pfParams[key] = float(pfParams[key])
    pfParams['rtRange'][0] = float(pfParams['rtRange'][0])
    pfParams['rtRange'][1] = float(pfParams['rtRange'][1])
    # Replace indexing by actual adduct names
    adducts = pandas.read_csv(args.adducts, index_col='index_do_not_alter')
    for pair in pfParams['negAdductsPairs']:
        pair[0] = _adduct_rename(adducts, pair[0])
        pair[1] = _adduct_rename(adducts, pair[1])
    for pair in pfParams['posAdductsPairs']:
        pair[0] = _adduct_rename(adducts, pair[0])
        pair[1] = _adduct_rename(adducts, pair[1])
    # "rtTolMultipler" parameter removed in version 2.0, so apply it
    # directly to "maxRTDiffAdjFrame" parameter as it was done
    pfParams['maxRTDiffAdjFrame'] = \
            (pfParams['maxRTDiffAdjFrame']
            * ast.literal_eval(oldParams['rtTolMultipler']))
    # Create the new PeakFilter parameters JSON file
    with open(os.path.join(dst, 'peakfilter.json'), 'w') as paramsFile:
        json.dump(pfParams, paramsFile, indent=4)
    # Set Amalgamator parameters based on PeakFilter's ones
    amalParams = {'numSamples': pfParams['numSamples'],
                  'firstSampleIndex': pfParams['firstSampleIndex'],
                  'mzFixedError': pfParams['mzFixedError'],
                  'mzPPMError': pfParams['mzPPMError'],
                  'maxRTDiffAdjFrame': pfParams['maxRTDiffAdjFrame']}
    # Create the new Amalgamator parameters JSON file
    with open(os.path.join(dst, 'amalgamator.json'), 'w') as paramsFile:
        json.dump(amalParams, paramsFile, indent=4)

if (__name__ == "__main__"):
    main()
