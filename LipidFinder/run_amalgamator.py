#!/usr/bin/env python

# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Read negative and positive input CSV/TSV/XLS/XLSX files and the
parameters JSON file, create the output folder and launch LipidFinder's
Amalgamator.
"""

import argparse
from datetime import datetime
import os

from LipidFinder import Amalgamator
from LipidFinder.Configuration import LFParameters
from LipidFinder.LFDataFrame import LFDataFrame
from LipidFinder._utils import normalise_path


def main():
    # Create the argument parser and parse the arguments
    parser = argparse.ArgumentParser(
            description="Run LipidFinder's Amalgamator.")
    parser.add_argument('-neg', '--negative', metavar='FILE', type=str,
                        required=True, help="negative data file")
    parser.add_argument('-pos', '--positive', metavar='FILE', type=str,
                        required=True, help="positive data file")
    parser.add_argument('-o', '--output', metavar='DIR', type=str,
                        help="folder where the output files will be stored")
    parser.add_argument('-p', '--params', metavar='FILE', type=str,
                        required=True, help="parameters JSON file")
    parser.add_argument('--timestamp', action='store_true',
                        help="add a timestamp to the output folder's name")
    parser.add_argument('--version', action='version',
                        version="LipidFinder 2.0")
    args = parser.parse_args()
    # Load parameters and input data
    parameters = LFParameters(module='amalgamator', src=args.params)
    negData = LFDataFrame(args.negative, parameters)
    posData = LFDataFrame(args.positive, parameters)
    # Check if the output directory exists. If not, create it.
    dst = args.output if (args.output) else ''
    if (args.timestamp):
        timestamp = datetime.utcnow().strftime('%Y%m%d_%H%M%S')
        dst += '_{0}'.format(timestamp) if (dst) else timestamp
    dst = normalise_path(dst)
    if (not os.path.isdir(dst)):
        os.makedirs(dst)
    # Run Amalgamator
    Amalgamator.amalgamate_data(negData, posData, parameters, dst)

if (__name__ == '__main__'):
    main()
