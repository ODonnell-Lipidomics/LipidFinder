# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Represent a DataFrame to be processed with LipidFinder's workflow."""

import glob
import logging
import os

import pandas


class LFDataFrame(pandas.core.frame.DataFrame):
    """A LFDataFrame object stores a dataframe to be used as input data
    in LipidFinder.

    The input data file(s) must comply with the following requirements:
      - The format must be: CSV, TSV, XLS or XLSX. For the last two the
        user can also specify the sheet to be read (or the list of
        sheets if a folder is given as 'src').
      - The first column contains an identifier for each row that is
        unique throughout every file.
      - There is one column named as "mzCol" parameter and another one
        as "rtCol" parameter.
      - Starting from the column index in "firstSampleIndex" parameter,
        every intensity column must follow. For instance, for 2 samples
        with 2 technical replicates, 1 quality control sample and 2
        solvents, the columns would be as follows:
            sample11 , sample12 , sample21 , sample22 , QC1 , sol1, sol2
        Ensure that samples with multiple technical replicates are given
        names in the format name1, name2, etc. such that each name is
        unique for each column. Replicates should be suffixed 1, 2, etc.

    Attributes:
        src  (Public[str])
            Source path where the data was loaded from.
        _resolution  (Private[int])
            Number of digits after the radix point in floats.

    Examples:
        LFDataFrame objects can be created in two different ways:
            >>> from Configuration import LFParameters
            >>> from LFDataFrame import LFDataFrame
            >>> params = LFParameters(module='peakfilter')
            >>> csvData = LFDataFrame('input_data.csv', params)
            >>> xlsData = LFDataFrame('input_data.xls', params, sheet=2)
            >>> folderData = LFDataFrame('/home/user/data/', params)

        After loading the required set of parameters, the data can be
        loaded from a single file ('csvData' and 'xlsData' examples) or
        from multiple files located in the same folder ('folderData'
        example). The latter is meant to be used to merge multiple files
        split by time ranges that represent a single run. The first and
        last retention time (RT) minutes of every file are trimmed as
        they are considered unreliable (except for the first and last
        minutes of the first and last files, respectively). The method
        supports overlap (after trimming), and the frames retained will
        be those from the file with the most frames for each overlapping
        minute.

        The number of decimal places to keep from the input m/z column
        can be changed assigning a value to 'resolution' variable. It
        has been predefined to 6, a standard value in high-resolution
        liquid-chromatography coupled to mass-spectrometry.
    """

    def __init__(self, src, parameters, resolution=6, sheet=0):
        # type: (str, LFParameters, int, object) -> LFDataFrame
        """Constructor of the class LFDataFrame.

        Keyword Arguments:
            src        -- source path where to load the data from
            parameters -- LipidFinder's parameters instance (can be for
                          any module)
            resolution -- number of decimal places to keep from m/z
                          column [default: 6]
            sheet      -- sheet number or list of sheet numbers to read
                          when input file(s) have XLS or XLSX extension
                          (zero-indexed position) [default: 0]
        """
        rtCol = parameters['rtCol']
        if (not os.path.isdir(src)):
            data = self._read_file(src, parameters, sheet)
        else:
            # Create a list of the input files in the source folder (in
            # alphabetical order)
            fileList = sorted(glob.iglob(os.path.join(src, '*.*')))
            if (len(fileList) == 0):
                raise FileNotFoundError("No files found in '{0}'".format(src))
            data = self._read_file(fileList[0], parameters, sheet[0])
            if (len(fileList) > 1):
                # Sort first dataframe by RT
                data.sort_values([rtCol], inplace=True, kind='mergesort')
                # Append "minute" column to the dataframe with the
                # integer part of the float values of its RT column
                timeCol = 'minute'
                data = data.assign(minute=data[rtCol].astype(int))
                # Since it is the first file, remove the frames
                # corresponding to the last minute
                data = data[data[timeCol] != data.iloc[-1][timeCol]]
                for index, filePath in enumerate(fileList[1:], start=1):
                    chunk = self._read_file(filePath, parameters, sheet[index])
                    # Sort next chunk dataframe by RT
                    chunk.sort_values([rtCol], inplace=True, kind='mergesort')
                    # Append "minute" column to the dataframe with the
                    # integer part of the float values of its RT column
                    chunk = chunk.assign(minute=chunk[rtCol].astype(int))
                    # Remove the frames of the first minute
                    chunk = chunk[chunk[timeCol] != chunk.iloc[0][timeCol]]
                    if (index < (len(fileList) - 1)):
                        # Since it is not the last file, remove the
                        # frames corresponding to the last minute
                        chunk = chunk[chunk[timeCol] != chunk.iloc[-1][timeCol]]
                    # Create a dataframe with the number of frames per
                    # minute for both the dataframe and the next chunk
                    overlap = pandas.DataFrame(
                            {'data': data.groupby(timeCol).size(),
                             'chunk': chunk.groupby(timeCol).size()}
                            ).fillna(0)
                    # Keep the minutes where the number of frames in the
                    # next chunk is higher than in the current dataframe
                    overlap = overlap[overlap['chunk'] > overlap['data']]
                    minutesToReplace = overlap.index.tolist()
                    if (minutesToReplace):
                        # Remove the dataframe frames to be replaced
                        data = data[~data[timeCol].isin(minutesToReplace)]
                        # Append chunk frames preserving the column
                        # order of the main dataframe
                        data = data.append(
                                chunk[chunk[timeCol].isin(minutesToReplace)],
                                ignore_index=True
                                )[data.columns.tolist()]
                # Drop "minute" column as it will be no longer necessary
                data.drop(timeCol, axis=1, inplace=True)
        # Rename first column if no name was given in the input file(s)
        data.rename(columns={'Unnamed: 0': 'id'}, inplace=True)
        # Sort dataframe by m/z and RT, and reset the indexing
        mzCol = parameters['mzCol']
        data.sort_values([mzCol, rtCol], inplace=True, kind='mergesort')
        data.reset_index(drop=True, inplace=True)
        # Adjust m/z column values to the machine's maximum float
        # resolution
        data[mzCol] = data[mzCol].apply(round, ndigits=resolution)
        super(LFDataFrame, self).__init__(data=data)
        self.src = src
        self._resolution = resolution

    def drop_empty_frames(self, module, parameters, means=False):
        # type: (str, LFParameters, bool) -> None
        """Remove empty frames from the dataframe and reset the index.

        An empty frame is a row for which every sample replicate or
        sample mean has a zero intensity.

        Keyword Arguments:
            module     -- module name to write in the logging file
            parameters -- LipidFinder's parameters instance (can be for
                          any module)
            means      -- check sample means instead of each sample
                          replicate? [default: False]
        """
        if (means):
            meanColIndexes = [i for i, col in enumerate(self.columns)
                                  if col.endswith('_mean')]
            if (parameters['numSolventReps'] > 0):
                # The first mean column is for the solvents
                firstIndex = meanColIndexes[1]
            else:
                firstIndex = meanColIndexes[0]
            lastIndex = meanColIndexes[-1]
        else:
            firstIndex = parameters['firstSampleIndex'] - 1
            lastIndex = firstIndex \
                        + (parameters['numSamples'] * parameters['numTechReps'])
        # Get the indices of all empty frames
        emptyFrames = self.iloc[:, firstIndex : lastIndex].eq(0).all(axis=1)
        indices = self[emptyFrames].index.tolist()
        if (indices):
            # Drop empty frames and reset the index
            self.drop(module, labels=indices, axis=0, inplace=True)
            self.reset_index(drop=True, inplace=True)

    def drop(self, module, **kwargs):
        # type: (str, ...) -> LFDataFrame
        """Wrapper of pandas.DataFrame.drop() with logging report.

        The report will be updated only if the labels correspond to
        rows, i.e. kwargs['axis'] == 0 (default value).

        Keyword Arguments:
            module  -- module name to write in the logging file
            *kwargs -- arguments to pass to pandas.DataFrame.drop()
        """
        # Create logger to print message to the log file
        logger = logging.getLogger(module)
        logger.setLevel(logging.INFO)
        if ((len(kwargs['labels']) > 0) and (kwargs.get('axis', 0) == 0)):
            idCol = self.columns[0]
            idList = [str(x) for x in sorted(self.loc[kwargs['labels'], idCol])]
            logger.info('%s: removed %d rows. IDs: %s', module, len(idList),
                        ','.join(idList))
        return super(LFDataFrame, self).drop(**kwargs)

    @staticmethod
    def _read_file(src, parameters, sheet):
        # type: (str, LFParameters, int) -> pandas.core.frame.DataFrame
        """Return a dataframe with the same content as the source file,
        but with retention time in minutes.

        The read function will be configured based on the file's
        extension. Accepted extensions: CSV, TSV, XLS, XLSX.

        Keyword Arguments:
            src        -- source file path
            parameters -- LipidFinder's parameters instance (can be for
                          any module)
            sheet      -- sheet number to read when the input file has
                          XLS or XLSX extension (zero-indexed position)
        """
        extension = os.path.splitext(src)[1].lower()[1:]
        # Load file based on its extension
        if (extension == 'csv'):
            data = pandas.read_csv(src, float_precision='high')
        elif (extension == 'tsv'):
            data = pandas.read_csv(src, sep='\t', float_precision='high')
        elif (extension in ['xls', 'xlsx']):
            data = pandas.read_excel(src, sheet_name=sheet)
        else:
            raise IOError(("Unknown file extension '{0}'. Expected: csv, tsv, "
                           "xls, xlsx").format(extension))
        if (('timeUnit' in parameters) and
            (parameters['timeUnit'] == 'Seconds')):
            rtCol = parameters['rtCol']
            data[rtCol] = data[rtCol].apply(lambda x: round(x / 60.0, 2))
        return data
