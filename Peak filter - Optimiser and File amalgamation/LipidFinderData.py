from __future__ import print_function
from builtins import str
from builtins import input
from builtins import range
from builtins import object
import time
import numpy as np
import pandas as pd
import ast
import os
import re
import xcmsReformat


class LipidFinderData(object):
    """A class to store all parameters, manage user input, input files, DataFrames from all
    processing stages and m/z and retention time tolerances.

    Attributes:
        adductAddition (bool): Adduct addition toggle. Allows the user to specify whether they wish to add the intensity of any adducts identified to the intensity of the primary mass  (True: Add adduct intensity to primary mass intensity; False: Do NOT add adduct intensity)
        broadContMinPoints (int): Minimum number of non-zero points to allow broad contaminant removal to go ahead
        broadContRSDCutOff (int): The target intensity RSD for the non outliers to be considerd similar enough for removal
        broadContrtSDCutOff (int): The target RT RSD for the non outliers to be considerd similar enough for removal
        broadContsdMult (int): Number of standard deviations a point need to differ from the mean of the others in the chromatogram to be considerd an outlier
        columnType (str): The type of column used for LC  (PO: Polar; NP: Non-polar)
        currentPath (str): Will be set by default to be the current working directory.
        destFolder (str): The name of the destination folder, can also be path from currentPath and name of destination folder.
        destPath (str): The full path to the detination folder including the destination folder.
        featureLevelMassAssignment (bool): This toggle allows the user to specify whether they wish to assign the masses of every frame to the highest intesity mass in the feature set, default is assignment at the mass group level  (TRUE: Re-assign every mass in a feature group to the mass of the highest intensity frame within the feature group; FALSE: Re-assign every mass in a mass cluster to the mass of the highest average intensity frame within the mass cluster)
        filePolarityMode (str): File Polarity Mode (P: Positive mode files; N: Negative mode files)
        fileSuffix (str): A unique per run suffix to be added to every output file composed of the column type, the polarity and the time and date the run began
        finalFileName (str): The name assigned to the final file before the suffix is added.
        firstQCRepOffset (int): The index of the first QC sample replicate column.
        firstRepOffset (int): The index of the first sample, first replicate column.
        firstSolRepOffset (int): The index of the first solvent sample replicate column.
        fullRawData (pandas DataFrame): All the input files are initially held in this DataFrame.
        intensitySignificanceCutOff (int): The level at which the intensity of a sample reading is significant enough to process.
        lastQCRepOffset (int): The index of the last QC sample replicate column.
        lastRepOffset (int): The index of the last sample, last replicate column.
        lastSolRepOffset (int): The index of the first solvent sample replicate column.
        lipidStackAddition (bool): Lipid stack addition toggle. Allows the user to specify whether they wish to add the intensity of any lipid stacks ions identified to the intensity of the primary mass  (True: Add adduct intensity to primary mass intensity; False: Do NOT add adduct intensity)
        listRawFileDataFrame (list): A list of DataFrames that make up the input files.
        maxStackGap (int): The maximum number missing values before a stack search will terminate for a particualr lipid or contaminant.
        mzFixedError (float): The fixed error allowable when observing a mass.
        mzSizeErrorPPM (float): The mass size dependant PPM error to add to the fixed error.
        negativeModeAdductList (list): The pairs of negative adducts relating to the index in the adducts.csv file
        numberOfFiles (int): The number of csv files that make up the input.
        numberOfQCReps (int): The number of QC replicates in the input file(s).
        numberOfSamples (int): The number of samples in the experiment.
        numberOfSolventReps (int): The number of solvent replicates in the input file(s).
        numberOfTechReps (int): The number of replicates of each sample.
        outlierHighIntensityRSD (int): The RSD to use when the mean average intensities are higher than the replicate mean intensity cut off (repMeanOutCorCutOff).
        outlierHighIntensityValue (int): The cut off point of the replicate means of a sample between using the lower RSD cut off (outlierLowIntensityRSD) and the higher RSD cut off (outlierHighIntensityRSD).
        outlierLowIntensityRSD (int): The RSD to use when the mean average intensities are lower than the replicate mean intensity cut off (repMeanOutCorCutOff).
        parameters (pandas DataFrame): A DataFrame produced by the import of the parameters.csv file
        peakAdjacentFrameMaxRT (float): The maximum time difference (mins) between a feature edge and an adjacent frame where the adjacent frame could be considered for inclusion in the same feature.
        peakConcatenateAllFrames (bool): This toggle indicates which peak concatenation should be performed (TRUE:Concatenate all peak frame intensities into the peak centre, FALSE:Concatenate only the most intense peak frame into the peak centre)
        peakMaxRTWidth (float): The maximum allowable retention time (mins) a single lipid peak can span
        peakMinFoldCutOff (float): The minimum fold difference greater than the adjacent candidate frame's intensity that the current frame intensity must be in order for the current frame to be considered part of a peak.
        positiveModeAdductList (list): The pairs of positive adducts relating to the index in the adducts.csv file
        processedDataFrame (pandas DataFrame): The current DataFrame, this will be processed next by any of the PeakFilter modules.
        QCHighRSD (int): Percentage upper relative standard deviation cut off used in QC calculation.
        QCLowRSD (int): Percentage lower relative standard deviation cut off used in QC calculation.
        removeAdducts (bool): This toggle allows the user to specify whether adduct removal should be executed  (TRUE: Remove adducts; FALSE: Do NOT Remove adducts)
        removeContaminants (bool): This toggle allows the user to specify whether contaminant removal should be executed  (TRUE: Remove contaminants; FALSE: Do NOT Remove contaminants)
        removeSolvent (bool): Solvent removal toggle (TRUE: Remove solvent intensity from sample replicates; FALSE: Do NOT Remove solvent from sample replicates)
        removeStacks (bool): This toggle allows the user to specify whether stack removal should be executed  (TRUE: Remove stacks; FALSE: Do NOT Remove stacks)
        retentionTimeHighCutOff (float): Frames with a retention time greater than this value will be removed from the dataset
        retentionTimeLowCutOff (float): Frames with a retention time less than this value will be removed from the dataset
        rtCorrectMeans (bool): This toggle allows the user to specify whether they want to RT correct the means of the sample replcates (for use with biological replcates ) (TRUE: RT correct means; FALSE: Do NOT RT correct means)
        rtCorrectStDev (float): The maximum deviation from the mean intensity of the target frame that the donor intensity may be in RT correction.
        rtTolMultipler (float): A multiplier for peakAdjacentFrameMaxRT to allow a smaller tolerance in certain circumstances (e.g. when looking for stacks)
        solventFoldCutOff (float): The minimum fold difference greater than the solvent intensity a sample replicate intensity must be in order to be considered significant enough to process further during solvent removal
        sourcePath (str): The location of the source files
        timeStamp (str): A time date date suffix
    """

    def __init__(self):
        self.timeStamp = time.strftime("_%Y%m%d-%H%M%S")
        self.__getParameters__()
        self.__setFileSuffix__()
        self.currentPath = os.getcwd()
        self.sourcePath = self.currentPath + os.sep + "sourceFiles"
        self.destFolder = "outputFiles"
        self.destPath = self.currentPath + os.sep + self.destFolder
        self.finalFileName = "finalOutputFile"
        self.numberOfFiles = 1

    def __getParameters__(self):
        """Set all parameters to the values imported from parameters.csv
        """
        # Read in default self.parameters csv file
        self.parameters = pd.read_table('parameters.csv', sep=',')

        self.firstRepOffset = int(self.parameters.ix[0][3])

        self.numberOfSamples = int(self.parameters.ix[1][3])
        self.numberOfTechReps = int(self.parameters.ix[2][3])

        self.setLastRepOffset()

        self.numberOfQCReps = int(self.parameters.ix[3][3])
        self.numberOfSolventReps = int(self.parameters.ix[4][3])

        self.filePolarityMode = str(self.parameters.ix[5][3])
        self.columnType = str(self.parameters.ix[6][3])

        self.QCLowRSD = float(self.parameters.ix[7][3])
        self.QCHighRSD = float(self.parameters.ix[8][3])

        self.removeSolvent = True if self.parameters.ix[
            9][3].upper() == "TRUE" else False
        self.solventFoldCutOff = float(self.parameters.ix[10][3])

        self.intensitySignificanceCutOff = int(self.parameters.ix[11][3])

        self.mzFixedError = float(self.parameters.ix[12][3])
        self.mzSizeErrorPPM = float(self.parameters.ix[13][3])

        self.peakMaxRTWidth = float(self.parameters.ix[14][3])
        self.peakMinFoldCutOff = float(self.parameters.ix[15][3])
        self.peakAdjacentFrameMaxRT = float(self.parameters.ix[16][3])
        self.peakConcatenateAllFrames = True if self.parameters.ix[
            17][3].upper() == "TRUE" else False

        self.removeContaminants = True if self.parameters.ix[
            18][3].upper() == "TRUE" else False
        self.removeAdducts = True if self.parameters.ix[
            19][3].upper() == "TRUE" else False
        self.adductAddition = True if self.parameters.ix[
            20][3].upper() == "TRUE" else False
        self.removeStacks = True if self.parameters.ix[
            21][3].upper() == "TRUE" else False
        self.maxStackGap = int(self.parameters.ix[22][3])
        self.lipidStackAddition = True if self.parameters.ix[
            23][3].upper() == "TRUE" else False

        # Potential change name
        self.rtTolMultipler = float(self.parameters.ix[24][3])

        self.outlierHighIntensityValue = int(self.parameters.ix[25][3])
        self.outlierLowIntensityRSD = int(self.parameters.ix[26][3])
        self.outlierHighIntensityRSD = int(self.parameters.ix[27][3])

        # Option to RT correct
        # Option to Outlier correct
        self.featureLevelMassAssignment = True if self.parameters.ix[
            28][3].upper() == "TRUE" else False

        # Adduct lists
        self.negativeModeAdductList = list(
            ast.literal_eval(self.parameters.ix[29][3]))
        self.positiveModeAdductList = list(
            ast.literal_eval(self.parameters.ix[30][3]))

        # Broad contaminant parameters
        self.broadContsdMult = float(self.parameters.ix[31][3])
        self.broadContMinPoints = int(self.parameters.ix[32][3])
        self.broadContRSDCutOff = int(self.parameters.ix[33][3])
        self.broadContrtSDCutOff = int(self.parameters.ix[34][3])

        # Retention time cut off allows user to consider frames within a subset
        # of the whole data, this reduces any error sieve errors (as wider
        # range can be calculated in sieve) and may help with sieve alignment
        # before LipidFinder. Frames before and after the cut offs will be used
        # to help identify frames within cut off
        self.retentionTimeLowCutOff = float(self.parameters.ix[35][3])
        self.retentionTimeHighCutOff = float(self.parameters.ix[36][3])

        self.rtCorrectStDev = float(self.parameters.ix[37][3])
        self.rtCorrectMeans = True if self.parameters.ix[
            38][3].upper() == "TRUE" else False

    # Suffix string for outputfiles to help user identify runs
    def __setFileSuffix__(self):
        self.fileSuffix = "_" + self.columnType + \
            "_" + self.filePolarityMode + self.timeStamp

    # function to determine whether or not the input file is from Sieve or
    # XCMS. Incorrect answer will prompt the function to call itself and ask
    # the user for a correct input.
    def checkSourceFileOrigin(self):
        """Prompts the user for the type of preprocessing performed on the input file(s)

        Returns:
            str: 'S' for SIEVE or 'X' for XCMS
        """
        fileOrigin = input(
            'Please specify Sieve or XCMS (enter "S" for Sieve or "X" for XCMS):')
        if fileOrigin.upper() == 'X' or fileOrigin.upper() == 'S':
            return fileOrigin.upper()
        else:
            print("\n")
            print("Please enter an 'X' or 'S'.")
            self.checkSourceFileOrigin()

    def checkDataParameters(self):
        """Run several prompts one after giving the user opportunity to verify key
        parameters
        """
        while True:
            overideSamplesFlag = input("Default data parameters are:-\n\n   Number of samples = %d\n   Number of replicates = %d\n   Number of QC replicates = %d\n   Number of solvent replicates = %d\n   File polarity mode = %s\n   Column type = %s\n\nChange these sample parameters (Y or N)? [N]: " %
                                       (self.numberOfSamples, self.numberOfTechReps, self.numberOfQCReps, self.numberOfSolventReps, self.filePolarityMode, self.columnType)) or 'N'
            if overideSamplesFlag.upper() == 'N':
                break
            elif overideSamplesFlag.upper() == 'Y':
                self.setNumSamples()
                self.__setNumberOfTechReps__()
                self.__setNumberOfQCReps__()
                self.__setNumberOfSolventReps__()
                self.__setFilePolarityMode__()
                self.__setColumnType__()
                break
            else:
                print("\nY or N required!")

    def setNumSamples(self):
        """Prompts the user for the number of samples

        Returns:
            int: The number of samples
        """
        while True:
            try:
                self.numberOfSamples = int(input(
                    "   How many samples? [%s]: " % self.numberOfSamples) or self.numberOfSamples)
                if self.numberOfSamples >= 1:
                    break
                else:
                    print ("   Integer >= 1 needed!")
            except ValueError:
                print ("   Integer >= 1 needed!")

    # Function to set number of technical replicates in the run and handle
    # input errors
    def __setNumberOfTechReps__(self):
        while True:
            try:
                self.numberOfTechReps = int(input(
                    "   How many technical replicates per sample? [%s]: " % self.numberOfTechReps) or self.numberOfTechReps)
                if self.numberOfTechReps >= 1:
                    break
                else:
                    print ("   Integer >= 1 needed!")
            except ValueError:
                print ("   Integer >= 1 needed!")

    # Function to set number of QC replicates in the run and handle input
    # errors
    def __setNumberOfQCReps__(self):
        while True:
            try:
                self.numberOfQCReps = int(input(
                    "   How many technical replicates of the QC sample? (use 0 for no QC samples) [%s]: " % self.numberOfQCReps) or self.numberOfQCReps)
                if self.numberOfQCReps >= 0:
                    break
                else:
                    print ("   Integer >= 0 needed!")
            except ValueError:
                print ("   Integer >= 0 needed!")

    # Function to set number of solvent replicates in the run and handle input
    # errors
    def __setNumberOfSolventReps__(self):
        while True:
            try:
                self.numberOfSolventReps = int(input(
                    "   How many technical replicates of the solvent? (use 0 for no solvent samples) [%s]: " % self.numberOfSolventReps) or self.numberOfSolventReps)
                if self.numberOfSolventReps >= 0:
                    break
                else:
                    print ("   Integer >= 0 needed!")
            except ValueError:
                print ("   Integer >= 0 needed!")

    # Function to set polarity of the file (POS or NEG)  and to handle input
    # errors
    def __setFilePolarityMode__(self):
        while True:
            tempFilePolarityMode = input(
                "   What is the polarity of the file? (use 'P' for Positive and 'N' for Negative) [%s]: " % self.filePolarityMode) or self.filePolarityMode
            upperTempFilePolarityMode = tempFilePolarityMode.upper()
            if upperTempFilePolarityMode in ['P', 'N']:
                self.filePolarityMode = upperTempFilePolarityMode
                self.__setFileSuffix__()
                break
            else:
                print("\n   'P' or 'N' required!")

    # Function to set column type (Polar or Non-Polar) and to handle input
    # errors
    def __setColumnType__(self):

        while True:
            tempColumnType = input(
                "   What is the column type? (use 'PO' for Polar and 'NP' for Non-polar) [%s]: " % self.columnType) or self.columnType
            upperTempColumnType = tempColumnType.upper()
            if upperTempColumnType in ['PO', 'NP']:
                self.columnType = upperTempColumnType
                self.__setFileSuffix__()
                break
            else:
                print("\n   'PO' or 'NP' required!")

    def setSourcePath(self, offset=0):
        """Gives user opprotunity to modify the path for the source files

        Args:
            offset (int, optional): The number of character offset from the left the prompt text will be printed

        """
        while True:
            tempSourcePath = input(
                offset * " " + "Specify any change to the default source path [%s]: " % self.sourcePath) or self.sourcePath
            if os.path.exists(tempSourcePath):
                self.sourcePath = tempSourcePath
                break
            else:
                print("Path does not exist!")

    def setDestFolder(self, offset=0):
        """Gives user opprotunity to modify the path for file output destination, if the
        folder doesen't exist it is created

        Args:
            offset (int, optional): The number of character offset from the left the prompt text will be printed
        """
        while True:
            tempDest = input(
                offset * " " + "Specify a sub-folder name to save the output files [%s]: " % self.destFolder) or self.destFolder

            # If the folder does not exist, try to create it
            if not os.path.exists(self.currentPath + os.sep + tempDest):
                try:
                    os.mkdir(tempDest)
                    self.destFolder = tempDest
                    self.destPath = self.currentPath + os.sep + self.destFolder
                    break
                except OSError:
                    print("Invalid folder name!")

            # If it does exist set the destPath to it
            else:
                self.destFolder = tempDest
                self.destPath = self.currentPath + os.sep + self.destFolder
                break

    def setFinalFileName(self, offset=0):
        """Gives user opportunity to modify the final file name suffix

        Args:
            offset (int, optional): The number of character offset from the left the prompt text will be printed
        """
        while True:
            tempFinalFileName = input(
                offset * " " + "Specify a final file name [%s]: " % self.finalFileName) or self.finalFileName
            # Check 'tempFinalFileName' for bad formatting
            if len(re.findall(r'[^A-Za-z0-9_-]', tempFinalFileName)) == 0:
                # If no illegal characters then change finalFileName
                self.finalFileName = tempFinalFileName
                break
            else:
                print("Invalid file name!")

    def setNumberOfFiles(self):
        """Prompts the user for the number of source files
        """
        while True:
            try:
                tempNumberOfFiles = int(input(
                    "How many files make up the run? [%d]: " % self.numberOfFiles) or self.numberOfFiles)
                if tempNumberOfFiles >= 1:
                    self.numberOfFiles = tempNumberOfFiles
                    break
                else:
                    print ("Integer >= 1 needed!")
            except ValueError:
                print ("Integer >= 1 needed!")

    # create a loop to prompt user to input file names according to number of
    # files to append together
    def importDataFrame(self, file_origin):
        """Prompts user for each input file building up the raw DataFrame as it loops
        through according to number of source files

        Args:
            file_origin (str): 'S' for SIEVE or 'X' for XCMS

        """
        # A list to hold each file's DataFrame
        self.listRawFileDataFrame = []

        # Use while loop to add DataFrames of input files to a list
        fileNo = 1
        while fileNo in range(1, self.numberOfFiles + 1):
            # Prompt user for file name
            fileName = self.sourcePath + os.sep + \
                input("File " + str(fileNo) + ":-\n\n    ")
            try:

                # Read file into DataFrame
                if file_origin == 'X':
                    rawFileDataFrame = xcmsReformat.reformat(fileName)

                else:
                    rawFileDataFrame = pd.read_table(fileName, sep=',')

                # Add additional series to frame listing file number before m/z
                # column, to facillitate record removal based on overlap.
                rawFileDataFrame.insert(
                    self.firstRepOffset - 2, 'File', "File " + str(fileNo))

                # Put time bucket (Retention time as an integer)
                rawFileDataFrame.insert(
                    self.firstRepOffset + 1, 'Time Bucket', rawFileDataFrame['Time'].astype(int))

                # For all files except the first and last remove the first and
                # last minute in the file
                if fileNo > 1 and fileNo < self.numberOfFiles:
                    rawFileDataFrame = rawFileDataFrame[(rawFileDataFrame['Time Bucket'] != rawFileDataFrame[
                                                         'Time Bucket'].max()) & (rawFileDataFrame['Time Bucket'] != rawFileDataFrame['Time Bucket'].min())]
                # For the first file where more than 1 file remove the last
                # minute in the file for overlap
                elif fileNo == 1 and self.numberOfFiles > 1:
                    rawFileDataFrame = rawFileDataFrame[rawFileDataFrame[
                        'Time Bucket'] != rawFileDataFrame['Time Bucket'].max()]
                # For last file where more than 1 file remove the first minute
                # in the file for overlap
                elif fileNo == self.numberOfFiles and self.numberOfFiles > 1:
                    rawFileDataFrame = rawFileDataFrame[rawFileDataFrame[
                        'Time Bucket'] != rawFileDataFrame['Time Bucket'].min()]
                # Where only 1 file do not remove any file frames
                # Append the new DataFrame to the list
                self.listRawFileDataFrame.append(rawFileDataFrame)
                # Move to next slot in DataFrame list
                fileNo += 1
            except IOError as e:
                print("\n    " + str(e) + "!")
                # Prompt user again for same file number

        # Make a copy of the DataFrame at index 0 in the list of input file
        # DataFrames, this will always exist
        self.fullRawData = self.listRawFileDataFrame[0].copy()

        # Loop through the rest of the list of input file DataFrames and
        # concatenate to the DataFrame from above
        for rawFileDataFrameIndex in range(1, len(self.listRawFileDataFrame)):
            self.fullRawData = pd.concat([self.fullRawData, self.listRawFileDataFrame[
                                         rawFileDataFrameIndex]], ignore_index=True)

        # Remove overlap
        # Count the number of each time bucket for each file
        countFullRawData = pd.DataFrame({'count': self.fullRawData.groupby(
            ['File', 'Time Bucket']).size()}).reset_index()
        # Get a groupby object to index the highest count for each time bucket
        # in countFullRawData
        timeBucketSelectFile = countFullRawData.groupby(
            'Time Bucket')['count'].agg(lambda col: col.idxmax())
        # A DataFrame listing which file to use for each time bucket
        fileTimeBucketToUse = countFullRawData.ix[
            timeBucketSelectFile, ['File', 'Time Bucket']]
        # Merge fileTimeBucketToUse with fullRawData leaving most abundant over
        # lapping file frames
        self.fullRawData = pd.merge(self.fullRawData, fileTimeBucketToUse, on=[
                                    'File', 'Time Bucket'])

        # 'Time Bucket' serves no further purpose so drop for clarity
        self.fullRawData.drop('Time Bucket', inplace=True, axis=1)
        # This DataFrame is the curent highest processed DataFrame
        self.processedDataFrame = self.fullRawData.copy()

        # Set the firstRepOffset to match the processedDataFrame
        self.incFirstRepOffset()

    def incFirstRepOffset(self, inc=1):
        """Increments the position of the sample, first replicate

        Args:
            inc (int, optional): The number of index positions to increment by

        """
        self.firstRepOffset += inc
        self.setLastRepOffset()

    def setLastRepOffset(self):
        """Set the index postion of the last sample, last increment

        """
        self.lastRepOffset = self.firstRepOffset + \
            (self.numberOfSamples * self.numberOfTechReps)

    def getFirstQCRepOffset(self):
        """Get the index of the first QC sample

        Returns:
            int: The index of the first QC sample
        """
        self.firstQCRepOffset = self.firstRepOffset + \
            (self.numberOfSamples * self.numberOfTechReps)
        return self.firstQCRepOffset

    def getLastQCRepOffset(self):
        """Get the index of the last QC sample

        Returns:
            int: The index of the last QC sample
        """
        self.lastQCRepOffset = self.firstRepOffset + \
            (self.numberOfSamples * self.numberOfTechReps) + self.numberOfQCReps
        return self.lastQCRepOffset

    def getFirstSolRepOffset(self):
        """Get the index of the first solvent sample

        Returns:
            int: The index of the first solvent sample
        """
        self.firstSolRepOffset = self.firstRepOffset + \
            (self.numberOfSamples * self.numberOfTechReps) + self.numberOfQCReps
        return self.firstSolRepOffset

    def getLastSolRepOffset(self):
        """Get the index of the last solvent sample

        Returns:
            int: The index of the last solvent sample
        """
        self.lastSolRepOffset = self.firstRepOffset + \
            (self.numberOfSamples * self.numberOfTechReps) + \
            self.numberOfQCReps + self.numberOfSolventReps
        return self.lastSolRepOffset

    def lowerMZLimit(self, mz, mzSizeErrorPPM=None):
        """Find the lower tolerance limit for a m/z

        Args:
            mz (float): A m/z
            mzSizeErrorPPM (None, optional): The option to use a different mass size dependant PPM error to the one imported from parameters.csv

        Returns:
            float: The lower tolerance limit for the queried m/z
        """
        mzSizeErrorPPM = mzSizeErrorPPM or self.mzSizeErrorPPM
        return round((mz - (self.mzFixedError + (mz * mzSizeErrorPPM * 0.000001))), 5)

    def upperMZLimit(self, mz, mzSizeErrorPPM=None):
        """Find the upper tolerance limit for a m/z

        Args:
            mz (TYPE): A m/z
            mzSizeErrorPPM (None, optional): The option to use a different mass size dependant PPM error to the one imported from parameters.csv

        Returns:
            float: The upper tolerance limit for the queried m/z
        """
        mzSizeErrorPPM = mzSizeErrorPPM or self.mzSizeErrorPPM
        return round((mz + (self.mzFixedError + (mz * mzSizeErrorPPM * 0.000001))), 5)

    def mzRange(self, mz, mzSizeErrorPPM=None):
        """The mass difference (da) from lower tolerance limit to the upper tolerance limit

        Args:
            mz (TYPE): A m/z
            mzSizeErrorPPM (None, optional): The option to use a different mass size dependant PPM error to the one imported from parameters.csv

        Returns:
            float: The size (da) of the range from the lower tolerance limit to the upper tolerance limit
        """
        mzSizeErrorPPM = mzSizeErrorPPM or self.mzSizeErrorPPM
        return self.upperMZLimit(mz, mzSizeErrorPPM) - self.lowerMZLimit(mz, mzSizeErrorPPM)

    def lowerRTLimit(self, rt, rtTolMultipler=1, peakAdjacentFrameMaxRT=None):
        """Find the lower tolerance limit for a retention time

        Args:
            rt (float): A retention time
            rtTolMultipler (int, optional): A multiplier for to facillitate adjustment RT tolerance
            peakAdjacentFrameMaxRT (None, optional): The option to use a different maximum time difference (mins) between a feature edge and an adjacent frame other than the one imported from parameters.csv

        Returns:
            float: The lower tolerance limit for the queried retention time
        """
        peakAdjacentFrameMaxRT = peakAdjacentFrameMaxRT or self.peakAdjacentFrameMaxRT
        return round((rt - (peakAdjacentFrameMaxRT * rtTolMultipler)), 5)

    def upperRTLimit(self, rt, rtTolMultipler=1, peakAdjacentFrameMaxRT=None):
        """Find the upper tolerance limit for a retention time

        Args:
            rt (float): A retention time
            rtTolMultipler (int, optional): A multiplier for to facillitate adjustment RT tolerance
            peakAdjacentFrameMaxRT (None, optional): The option to use a different maximum time difference (mins) between a feature edge and an adjacent frame other than the one imported from parameters.csv

        Returns:
            float: The upper tolerance limit for the queried retention time
        """
        peakAdjacentFrameMaxRT = peakAdjacentFrameMaxRT or self.peakAdjacentFrameMaxRT
        return round((rt + (peakAdjacentFrameMaxRT * rtTolMultipler)), 5)
