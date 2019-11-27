#!/usr/bin/env python

# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Configure the parameters for the selected LipidFinder's module."""

import argparse

from LipidFinder.Configuration.LFParametersCLI import LFParametersCLI


def main():
    # Create the argument parser and parse the arguments
    parser = argparse.ArgumentParser(
            description="Run LipidFinder's parameter configuration.")
    parser.add_argument('-m', '--module', required=True,
                        choices=['peakfilter', 'amalgamator', 'mssearch'],
                        help="LipidFinder's module to configure")
    parser.add_argument('-p', '--params', metavar='FILE', type=str,
                        help="parameters JSON file")
    parser.add_argument('--version', action='version',
                        version="LipidFinder 2.0")
    args = parser.parse_args()
    # Run the parameter configuration through a command line
    LFParametersCLI(module=args.module, src=args.params)

if (__name__ == '__main__'):
    main()
