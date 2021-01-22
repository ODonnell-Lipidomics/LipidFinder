# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to calculate the False Discovery Rate:
    > get_fdr():
        Return the False Discovery Rate (FDR) of the dataset following a
        target-decoy strategy.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import FalseDiscoveryRate
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> FalseDiscoveryRate.get_fdr(data, parameters)
"""

import os
import requests

import pandas
from requests_toolbelt.multipart.encoder import MultipartEncoder

from LipidFinder._py3k import StringIO, range


LIPIDMAPS_URL = 'https://www.lipidmaps.org/tools/ms/py_bulk_search.php'
# Maximum number of m/z values to send at once to LIPID MAPS
BATCH_SIZE = 150


def get_fdr(data, parameters):
    # type: (LFDataFrame, LFParameters) -> float
    """Return the False Discovery Rate (FDR) of the dataset following a
    target-decoy strategy.

    The value is calculated based on the number of m/z values of 'data'
    found in the COMP_DB database from LIPID MAPS, and the number of m/z
    values of 'data' found in a decoy database, created adding 0.5 Da to
    every m/z in COMP_DB (a very rare lipid mass defect). FDR is equal
    to the number of decoy hits divided by the number of target hits.

    Keyword arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Get the list of unique m/z values from 'data'
    mzList = data[parameters['mzCol']].unique().tolist()
    # Set the target adducts
    if (parameters['polarity'] == 'Positive'):
        targetAdducts = ("M+H,M+H-H2O,M+2H,M+3H,M+4H,M+NH4,M+Ag,M+Na,M+2Na,M+K,"
                         "M+2K,M+Li,M+2Li")
    else:
        targetAdducts = 'M-H,M-CH3,M-2H,M-3H,M-4H,M.F,M.HF2,M.Cl,M.OAc,M.HCOO'
    # Get the number of matches in batches to balance the number of
    # requests and the amount of information requested
    numTargetHits = 0
    numDecoyHits = 0
    for start in range(0, len(mzList), BATCH_SIZE):
        mzBatch = mzList[start : start + BATCH_SIZE]
        # Get a string with one m/z per line (text file alike)
        mzStr = os.linesep.join(map(str, mzBatch))
        numTargetHits += _get_num_matches('COMP_DB', mzStr, targetAdducts)
        numDecoyHits += _get_num_matches('COMP_DB_5', mzStr, targetAdducts)
    # Raise an exception if there are no matches in the target database
    if (numTargetHits == 0):
        raise ValueError(("No matches found in the target database. The FDR "
                          "cannot be computed."))
    # FDR = numDecoyHits / numTargetHits
    return float(numDecoyHits) / numTargetHits


def _get_num_matches(db, mzStr, adducts, tolerance=0.001):
    # type: (str, str, str, float) -> int
    """Return the number of matches for the selected database and
    parameters.

    Keyword Arguments:
        db        -- LIPID MAPS' database
        mzStr     -- string with one m/z per line (text file alike)
        adducts   -- list of adducts separated by commas
        tolerance -- mass tolerance in Daltons (+/- to each m/z)
                     [default: 0.001]
    """
    # Create the target data package with the query
    mpData = MultipartEncoder(
            fields={'CHOICE': db, 'ion': adducts,
                    'file': ('file', StringIO(mzStr), 'text/plain'),
                    'tol': str(tolerance), 'sort': 'DELTA'})
    # Request the table containing the matches from LIPID MAPS
    try:
        response = requests.post(LIPIDMAPS_URL, data=mpData,
                                 headers={'Content-Type': mpData.content_type})
    except:
        raise Exception(("Connection error with the database. Please check your"
                         " network and try again after a few minutes."))
    # Process the response to create a dataframe and count the number of
    # hits (number of m/z that got at least one match in the database)
    if (len(response.text) == 0):
        return 0
    else:
        matches = pandas.read_csv(StringIO(response.text), sep='\t',
                                  engine='python', index_col=False)
        if (matches.empty):
            return 0
        else:
            return len(matches['Input Mass'].unique())
