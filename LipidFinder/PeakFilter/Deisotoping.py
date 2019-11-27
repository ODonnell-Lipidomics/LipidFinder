# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to find and remove isotopes:
    > remove_isotopes():
        Remove isotopes of parent analytes.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from PeakFilter import Deisotoping
    >>> parameters = LFParameters('peakfilter', 'parameters.json')
    >>> data = LFDataFrame('dataset.csv', parameters)
    >>> Deisotoping.remove_isotopes(data, parameters)
"""

import numpy

from LipidFinder._py3k import range, viewvalues, viewitems
from LipidFinder._utils import mz_tol_range, rt_tol_range


# Difference between C13 (13.003354838 u) and C12 (12 u) masses
ISO_OFFSET = 1.003354838


def remove_isotopes(data, parameters):
    # type: (LFDataFrame, LFParameters) -> None
    """Remove isotopes of parent analytes.

    Keyword Arguments:
        data       -- LFDataFrame instance
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    mzCol = parameters['mzCol']
    rtCol = parameters['rtCol']
    # Calculate the location of sample columns based on the current
    # state of the dataframe (before adding isotope annotation)
    firstSampleCol = len(data.columns) - parameters['numSamples']
    lastSampleCol = len(data.columns)
    for i in range(firstSampleCol, lastSampleCol):
        # Create an array from 'data' with m/z, retention time, the
        # samples' intensity mean and index per row
        array = numpy.stack((data[mzCol].values, data[rtCol].values,
                             data.iloc[:, i], data.iloc[:, 0].values), axis=-1)
        tagArray = _detect_sample_isotopes(array, parameters)
        # Set the intensity of the sample detected isotopes to 0
        colName = data.columns[i]
        isoColName = colName + '_isotopes'
        data.insert(len(data.columns), isoColName, tagArray)
        if (parameters['removeIsotopes']):
            data.loc[data[isoColName].str.contains('M\+'), colName] = 0.0
    if (parameters['removeIsotopes']):
        # Drop empty frames, i.e. isotope frames found in every sample
        data.drop_empty_frames(
                'Isotope removal (isotopes found in every sample)', parameters,
                True)


def _detect_sample_isotopes(array, parameters):
    """Return an array with the tagged parents and their corresponding
    isotopes in the same order as in the given sample array.

    Keyword Arguments:
        array      -- array with m/z, retention time (RT), sample's
                      intensity mean and index of the original dataframe
        parameters -- LipidFinder's PeakFilter parameters instance
    """
    # Get the corresponding symbol for the polarity of the data (+ or -)
    polSign = '+' if (parameters['polarity'].lower() == 'positive') else '-'
    # Create an array of empty strings that will contain the tagged
    # parents and their corresponding isotopes
    tagArray = numpy.full(len(array), '', dtype=object)
    # Loop over each m/z to search for isotopes
    isotopesIndex = set()
    for index in range(0, len(array)):
        # Skip if frame has already been identified as an isotope
        if (array[index, 3] in isotopesIndex):
            continue
        for isoPeak in range(1, parameters['numIsotopes'] + 1):
            parentMZ = array[index, 0]
            tagID = int(array[index, 3])
            # Get the first and last indexes of the frames that are
            # within the first isotope m/z range for the current analyte
            isotopeMZ = parentMZ + ISO_OFFSET * isoPeak
            minMZ, maxMZ = mz_tol_range(isotopeMZ, parameters['mzFixedError'],
                                        parameters['mzPPMError'])
            mzMatches = numpy.searchsorted(array[:, 0], [minMZ, maxMZ])
            if (mzMatches[0] == mzMatches[1]):
                # Have not found any analyte with an isotope-like m/z
                if (isoPeak == 1):
                    # The first isotope must exists to search for others
                    break
                else:
                    continue
            # Filter m/z matches with the same RT as the parent
            parentRT = array[index, 1]
            minRT, maxRT = rt_tol_range(parentRT,
                                        parameters['maxRTDiffAdjFrame'])
            rtMatches = numpy.where(
                    (array[mzMatches[0] : mzMatches[1], 1] >= minRT)
                    & (array[mzMatches[0] : mzMatches[1], 1] <= maxRT))[0]
            if (len(rtMatches) == 0):
                # No candidates are within the same RT
                if (isoPeak == 1):
                    # The first isotope must exists to search for others
                    break
                else:
                    continue
            # Resultant indexes are based on the previous search
            rtMatches += mzMatches[0]
            # Filter the candidate isotopes by intensity
            parentInten = array[index, 2]
            # The intensity range coefficients vary depending on the
            # isotope number
            if (isoPeak == 1):
                # Get an estimated maximum number of C in the molecule
                numC = round(parentMZ / 12)
                # Calculate isotopic distribution based on polynomial
                # expansion
                baseIntensity = parentInten * (numC ** 1.3) * 0.002
                minIntensity = baseIntensity * parameters['isoIntensityCoef'][0]
                maxIntensity = baseIntensity * parameters['isoIntensityCoef'][1]
            elif (isoPeak == 2):
                # Get an estimated maximum number of C in the molecule
                numC = round(parentMZ / 12)
                # Calculate isotopic distribution based on polynomial
                # expansion
                baseIntensity = parentInten * (numC ** 1.7) * 0.0001
                minIntensity = baseIntensity * parameters['isoIntensityCoef'][0]
                maxIntensity = baseIntensity * parameters['isoIntensityCoef'][1]
            else:
                # Calculate isotopic distribution with the same formula
                # as CAMERA (from XCMS)
                minIntensity = parentInten * float('1e-{0}'.format(isoPeak + 2))
                maxIntensity = parentInten * 2
            isotopes = numpy.where((array[rtMatches, 2] >= minIntensity)
                                   & (array[rtMatches, 2] <= maxIntensity))[0]
            if (len(isotopes) == 0):
                # No candidates have an intensity within expected range
                if (isoPeak == 1):
                    # The first isotope must exists to search for others
                    break
                else:
                    continue
            # Resultant indexes are based on the previous search
            isotopes += rtMatches[0]
            # Tag the analyte as isotope and save its index to avoid
            # checking it as parent of other analytes
            tagArray[isotopes] = '[{0}][M+{1}]{2}'.format(tagID, isoPeak,
                                                          polSign)
            isotopesIndex.update(array[isotopes, 3])
        else:
            # Tag the analyte as parent
            tagArray[index] = '[{0}][M]{1}'.format(tagID, polSign)
    return tagArray
