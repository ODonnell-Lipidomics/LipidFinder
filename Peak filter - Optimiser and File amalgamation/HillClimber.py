from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import random
import pandas as pd
import numpy as np
import collections as col
from math import e

import solventCalcs as sol
import clustering as clus
import peakFinder as pf


class HillClimber(object):
    """Summary

    Attributes:
        all_parameters (list): A list of 4 lists of possible values for each parameter
        current_parameters (list): The current highest scoring list of parameters in the current hill climbing iteration
        current_parameters_llf_count (int): The number of lipid-like features found by the 'current_parameters'
        current_parameters_score (float): The score found when 'current_parameters' were tested
        cycle_1opt_best_parameters (list): The highest scoring parameter set found after starting the current 1-opt cycle.
        cycle_1opt_best_parameters_llf_count (int): The number of lipid-like features found after starting the current 1-opt cycle by using the using 'cycle_1opt_best_parameters'
        cycle_1opt_best_parameters_score (float): The highest score found after starting the current 1-opt cycle using 'cycle_1opt_best_parameters'
        cycle_1opt_count (int): Number of 1opt cycles completed in each the hill climber iteration
        debug (bool): If True then progress will be printed to the screen during run time
        export_file_name (str): The name given to the exported csv file of hill climbing results
        hc_parameters_dimensions (str): The name of the hill climbing parameters dimensions file
        hc_run_count (int): A count of the number of random restart hill climbing iterations already completed
        hc_run_limit (int): The number of random restart hill climbing iterations to be performed
        lfd (LipidFinderData): A LipidFinderData object
        max_cycle_1opt_count (int): The maximum number of 1-opt parameter cycles to in each hill climber iteration
        overall_best_parameters (list): The current highest scoring list of parameters in ALL hill climbing iterations
        overall_best_parameters_llf_count (int): The number of lipid-like features found by the 'overall_best_parameters'
        overall_best_parameters_score (float): The score found when 'overall_best_parameters' were tested
        pDF (pandas DataFrame): The single replicate data set processed using the current 1-opt parameters
        ps_change_count (int): A tally of the number of times the best parameter set changes in each hill climbing run
        results (pandas DataFrame): A DataFrame to record the results of the hill climbing process as it proceeds
        target_len (int): The number of target features in the target file
        tDF (pandas DataFrame): The target features csv file imported into a DataFrame
        timer (Timer): A Timer object to allow the user monitor progress at run time
        tInt (numpy array): The array of target intensities
        tMZ (numpy array): The array of target m/z
        tRep (str): The target replicate name
        tRT (numpy array): The array of target retention time values
        visited (collections OrderedDict()): A dictionary of visited parameters sets along with their scores and lipid-like features found
    """

    def __init__(self, lfd, tDF, tRep, timer, hc_run_limit=1, max_cycle_1opt_count=5, debug=True):

        self.lfd = lfd
        self.tDF = tDF
        self.tRep = tRep
        self.timer = timer
        self.hc_run_limit = hc_run_limit
        self.max_cycle_1opt_count = max_cycle_1opt_count
        self.debug = debug

        # Get number of target features in target file
        self.target_len = len(self.tDF)

        # Get dataframe to process
        self.pDF = self.__get_pDF__()

        # Set the default hc_parameters_dimensions file name
        self.hc_parameters_dimensions = 'hc_parameters_dimensions'

        # Set the default export file name
        self.export_file_name = 'hc_exported_file'

        # Update sample parameters in lfd as current processedDataFrame will be
        # overwritten with new processedDataFrame consisting of the single
        # sample replicate that the target list has been derived from.
        self.lfd.numberOfSamples = 1
        self.lfd.numberOfTechReps = 1
        self.lfd.numberOfQCReps = 0
        self.lfd.numberOfSolventReps = 1

        # Update lfd with pDF as processedDataFrame
        self.__reset_lfd_pDF__()

        # Get relevant arrays from target file data frame
        self.tMZ = self.tDF['MZ'].values
        self.tRT = self.tDF['Time'].values
        self.tInt = self.tDF['Intensity'].values

        # The list of parameter lists
        self.all_parameters = self.__get_all_parameters__()

        # Select a valid random starting state, done here so starting parameter
        # can be optionally changed by user prior to processing if necessary
        self.current_parameters = self.__set_random_parameters__()

    # This returns a subset of the lfd.processedDataframe. Only 1 replicate
    # and solvent mean are retained and renamed to give an easier to
    # manipulate file
    def __get_pDF__(self):
        # Get mean solv
        mSolv = sol.getSolvMeanName(self.lfd)

        # Get target replicate name
        nTRep = self.tRep

        df = self.lfd.processedDataFrame[
            ['id', 'MZ', 'Time', nTRep, mSolv]].copy()
        return df.rename(columns={nTRep: 'Intensity', mSolv: 'SolventMean'})

    # Make dataframe to process the processedDataFrame and change the lfd
    # settings to reflect this
    def __reset_lfd_pDF__(self):
        self.lfd.processedDataFrame = self.pDF.copy()
        self.lfd.firstRepOffset = 3
        self.lfd.setLastRepOffset()  # check for lastRepOffset in all funcs

    # Creates 2d list of all possible values for each tunable parameters from
    # specified or default hc_parameters_dimensions file.
    def __get_all_parameters__(self, hc_parameters_dimensions=None):

        if hc_parameters_dimensions is None:
            hc_parameters_dimensions = self.hc_parameters_dimensions

        # Read in the parameters dimensions (pd) file to a pandas dataframe
        HCAllParameters = pd.read_table(
            hc_parameters_dimensions + '.csv', sep=',')

        # For each tunable parameter get the Lower Limit, Higher Limit and
        # Increments
        pd_adjacentFrameMaxTimeDist = HCAllParameters.ix[HCAllParameters[
            'Parameter'] == 'adjacentFrameMaxTimeDist', 2:5].values[0]
        pd_peakTimeWidthCutOff = HCAllParameters.ix[HCAllParameters[
            'Parameter'] == 'peakTimeWidthCutOff', 2:5].values[0]
        pd_minPeakFoldCutOff = HCAllParameters.ix[HCAllParameters[
            'Parameter'] == 'minPeakFoldCutOff', 2:5].values[0]
        pd_mzSizeErrorPPM = HCAllParameters.ix[HCAllParameters[
            'Parameter'] == 'mzSizeErrorPPM', 2:5].values[0]

        # Use comprehension to formulate all possible lists parameter values
        # (pv)
        pv_adjacentFrameMaxTimeDist = [round(x * pd_adjacentFrameMaxTimeDist[2] + pd_adjacentFrameMaxTimeDist[0], 3) for x in range(
            int(round(old_div((pd_adjacentFrameMaxTimeDist[1] - pd_adjacentFrameMaxTimeDist[0]), pd_adjacentFrameMaxTimeDist[2]) + 1)))]
        pv_peakTimeWidthCutOff = [round(x * pd_peakTimeWidthCutOff[2] + pd_peakTimeWidthCutOff[0], 3) for x in range(
            int(round(old_div((pd_peakTimeWidthCutOff[1] - pd_peakTimeWidthCutOff[0]), pd_peakTimeWidthCutOff[2]) + 1)))]
        pv_minPeakFoldCutOff = [round(x * pd_minPeakFoldCutOff[2] + pd_minPeakFoldCutOff[0], 3) for x in range(
            int(round(old_div((pd_minPeakFoldCutOff[1] - pd_minPeakFoldCutOff[0]), pd_minPeakFoldCutOff[2]) + 1)))]
        pv_mzSizeErrorPPM = [round(x * pd_mzSizeErrorPPM[2] + pd_mzSizeErrorPPM[0], 3) for x in range(
            int(round(old_div((pd_mzSizeErrorPPM[1] - pd_mzSizeErrorPPM[0]), pd_mzSizeErrorPPM[2]) + 1)))]

        return [pv_adjacentFrameMaxTimeDist, pv_peakTimeWidthCutOff, pv_minPeakFoldCutOff, pv_mzSizeErrorPPM]

    # Create a random, valid parameter state
    def __set_random_parameters__(self,):
        ps = []
        ps.append(random.choice(self.all_parameters[0]))

        # Ensure peakTimeWidthCutOff is >= 3* adjacentFrameMaxTimeDist
        ps.append(random.choice(
            [value for value in self.all_parameters[1] if value >= 3 * ps[0]]))
        ps.append(random.choice(self.all_parameters[2]))
        ps.append(random.choice(self.all_parameters[3]))
        return ps

    # Hill climbing with a multiple runs
    def hc_mr(self):
        """Runs the hill climbing process

        """
        # Store results of each scoring
        self.results = pd.DataFrame()

        # Number of hill climb (hc) runs
        self.hc_run_count = 1

        # Number of 1opt cycles completed in each hc run
        self.cycle_1opt_count = 0

        # Number of times the best parameter set changes in each hc run
        self.ps_change_count = 0

        # Score the starting parameters (self.current_parameters)
        self.current_parameters_score, self.current_parameters_llf_count = self.__test_fitness__(
            self.current_parameters)

        # Set overall best parameters and scores to starting parameter set and
        # score
        self.overall_best_parameters = list(self.current_parameters)
        self.overall_best_parameters_score = self.current_parameters_score
        self.overall_best_parameters_llf_count = self.current_parameters_llf_count

        # Run hc repeated times
        while self.hc_run_count <= self.hc_run_limit:

            # Create a list of parameter sets already quantified, this will
            # help us get past shoulders where traditional hill climbing would
            # get stuck owing to not being able to find a better solution and
            # allowing us to move to an equally scoring parameter set in the
            # state space may cause an infinite loop back and forth on the
            # shoulder. This is stopped by using a visited list so we cannot
            # revisit a previously visited state and so will search for a
            # higher or equally scoring parameter set that has not been
            # visited. If there are none ie (a plateau) then programe will
            # exit, also storage space for results
            self.visited = col.OrderedDict()

            # Set best from previous 1opt cycle parameters to current hc run
            # starting parameter set and score
            self.cycle_1opt_best_parameters = list(self.current_parameters)
            self.cycle_1opt_best_parameters_score = self.current_parameters_score
            self.cycle_1opt_best_parameters_llf_count = self.current_parameters_llf_count

            # Append starting data to visited list
            self.__append_visited__(self.cycle_1opt_best_parameters,
                                    self.cycle_1opt_best_parameters_score, self.current_parameters_llf_count)

            # Print to screen
            if self.debug:
                print(
                    "\n\nStart of hill climbing run %d\n******************************\n" % self.hc_run_count)
                print("  Initial parameters = %s" %
                      str(self.cycle_1opt_best_parameters))
                print("  Initial score = %s" %
                      str(self.cycle_1opt_best_parameters_score))
                print("  Initial lipid like features = %d" %
                      self.cycle_1opt_best_parameters_llf_count)

            # Run through hill climbing once
            self.__perform_hc__()

            # Append the full results of the hill climb to a pandas DataFrame
            # before it's reset for the next hill climb
            self.__append_results__()

            # If best parameters in last hill climb iteration score equal to or
            # better than current overall best parameters reset current overall
            # best parameters and score
            if self.cycle_1opt_best_parameters_score >= self.overall_best_parameters_score:
                self.overall_best_parameters = list(
                    self.cycle_1opt_best_parameters)
                self.overall_best_parameters_score = self.cycle_1opt_best_parameters_score
                self.overall_best_parameters_llf_count = self.cycle_1opt_best_parameters_llf_count

            if self.debug:
                print(
                    "\nEnd of hill climbing run %d\n****************************" % self.hc_run_count)
                print("\n    Best parameters found in this hill climbing run %s" % str(
                    self.cycle_1opt_best_parameters))
                print("    Best parameters score in this hill climbing run  = %s" % str(
                    self.cycle_1opt_best_parameters_score))
                print("    Best parameters lipid like features in this hill climbing run = %d\n" %
                      self.cycle_1opt_best_parameters_llf_count)
                print("**************************************************************************\n**************************************************************************")
                print("\nCurrent overall best parameters found %s" %
                      str(self.overall_best_parameters))
                print("Current overall best parameters score %s" %
                      str(self.overall_best_parameters_score))
                print("Current overall best parameters lipid like features = %d" %
                      self.overall_best_parameters_llf_count)

            # Reset variables
            self.hc_run_count += 1
            self.cycle_1opt_count = 0
            self.ps_change_count = 0
            self.current_parameters = self.__set_random_parameters__()
            self.current_parameters_score, self.current_parameters_llf_count = self.__test_fitness__(
                self.current_parameters)

        self.__export_results__()

    # Perform 1 full hill climbing run
    def __perform_hc__(self):

        # Perform 1opt scoring cycles up to the maximum limit
        while self.cycle_1opt_count < self.max_cycle_1opt_count:

            self.cycle_1opt_count += 1

            if self.debug:
                print("\n  Start of 1opt cycle number %d\n  ============================" %
                      self.cycle_1opt_count)

            # Run through a single 1opt cycle
            self.__score_cycle_1opt__()

            # Check if a full 1opt has been performed without a parameter set
            # change, either stuck at local max or visited every 1opt set
            # previously (no move possible)
            if self.current_parameters != self.cycle_1opt_best_parameters:
                self.cycle_1opt_best_parameters = self.current_parameters
                self.cycle_1opt_best_parameters_score = self.current_parameters_score
                self.cycle_1opt_best_parameters_llf_count = self.current_parameters_llf_count

                if self.debug:
                    print("\n  End of 1opt cycle number %d\n  ==============================" %
                          self.cycle_1opt_count)
                    print("\n    Current best parameters %s" %
                          str(self.cycle_1opt_best_parameters))
                    print("    Current best parameters score %s" %
                          str(self.cycle_1opt_best_parameters_score))
                    print("    Current best parameters lipid like features = %d" %
                          self.cycle_1opt_best_parameters_llf_count)

            else:
                if self.debug:
                    print(
                        "\n  Hill climb run ended owing to no move for a full cycle_1opt\n")
                break

        else:
            if self.debug:
                print(
                    "\n  Hill climb run ended owing to reaching max_cycle_1opt_count\n")

    # Run through a single 1opt cycle
    def __score_cycle_1opt__(self):
        self.timer.mark()
        for p_ind in range(len(self.current_parameters)):
            if self.debug:
                print("\n    Parameter %d\n    ------------" % p_ind)

            # Get 1opt parameter set for current parameter index
            ps_1opt = self.__get_ps_1opt__(p_ind)
            self.__score_ps_1opt__(ps_1opt)

    # Process 1 individual 1opt parameter set
    def __score_ps_1opt__(self, ps_1opt):

        # Score each ps in ps_1opt, if empty because all in visited list, next
        # set of ps_1opt will be scored, if at end of full __score_cycle_1opt__ no
        # change to current_parameters has been made then end, no improvement
        # possible ever from here
        for ps in ps_1opt:

            # Score current parameter set and get lipid like features returned
            score, llf_count = self.__test_fitness__(ps)

            # Is current score higher than any previous scores for current hill
            # climb run?
            if score > self.current_parameters_score:

                # Allows exploration of plateaus, will only stop when no
                # parameters to explore because all seen before, as using then
                # curSSB2b then if same eventually we'll loop back to beginning
                # visiting all possible
                self.current_parameters = ps
                self.current_parameters_score = score
                self.current_parameters_llf_count = llf_count
                self.ps_change_count += 1

            # End of scoring
            self.timer.mark()

            # Append current data to visited list
            self.__append_visited__(ps, score, llf_count)

            if self.debug:
                print("\n      Parameter set %s" % str(ps))
                print("      --Score = %s" % str(score))
                print(
                    "      --Lipid like features = %d\n      -------------------------------------" % llf_count)

    def __append_visited__(self, parameters, score, llf_count):
        self.visited.update([(tuple(parameters), [self.hc_run_count, self.cycle_1opt_count,
                                                  self.ps_change_count, score, llf_count, self.timer.total()])])

    # Get score for a parameter set
    def __test_fitness__(self, parameters):

        # Run parameters through clustering and peak analysis in LipidFinder
        self.__proc_LF__(parameters)

        # Get an index array of lipid like features (intensity indices which
        # have a non zero value)
        llf_index = np.nonzero(self.lfd.processedDataFrame[
                               'Intensity'].values)[0]

        # Get copies of MZ, RT and intensity arrays from
        # self.lfd.processedDataFrame for llf (where intensities > 0)
        pMZ = np.copy(self.lfd.processedDataFrame['MZ'].values[llf_index])
        pRT = np.copy(self.lfd.processedDataFrame['Time'].values[llf_index])
        pInt = np.copy(self.lfd.processedDataFrame[
                       'Intensity'].values[llf_index])

        # Reset self.lfd.processedDataFrame to state at the beginning of
        # current function
        self.__reset_lfd_pDF__()

        # Used to build up score of parameters
        score = 0

        # Loop through each target in target list
        for t in range(self.target_len):

            # Find matches within retention time and m/z tolerance, more than 1
            # candidate is possible sometimes
            candidate_array = np.where(((pMZ >= self.lfd.lowerMZLimit(self.tMZ[t])) & (pMZ <= self.lfd.upperMZLimit(
                self.tMZ[t]))) & ((pRT >= self.lfd.lowerRTLimit(self.tRT[t])) & (pRT <= self.lfd.upperRTLimit(self.tRT[t]))))[0]

            # Variable to record best scoring candidate
            candidate_max_score = 0

            # Score each candidate found and retain best score
            for c in candidate_array:

                # Find offset Score 0.0 - 1.0 for each dimension, best is 0
                MZ_offset = old_div(abs(
                    self.tMZ[t] - pMZ[c]), (self.lfd.upperMZLimit(self.tMZ[t]) - self.tMZ[t]))
                RT_offset = old_div(abs(
                    self.tRT[t] - pRT[c]), (self.lfd.upperRTLimit(self.tRT[t]) - self.tRT[t]))
                Int_offset = old_div(min(
                    abs(self.tInt[t] - pInt[c]), self.tInt[t]), self.tInt[t])

                # Get gaussian curve height (0.0 - 1.0) for each dimension (less deviation will score closer to 1.0 and more will score closer 0.0), https://en.wikipedia.org/wiki/Gaussian_function
                # a = height of the curve's peak = 1
                # b = position of centre = 0
                # c = standard deviation = 0.5, so 95% contained within +1 and -1
                # x = deviation
                gaussMZScore = e**-(old_div(MZ_offset**2, 0.5))
                gaussRTScore = e**-(old_div(RT_offset**2, 0.5))
                gaussIntScore = e**-(old_div(Int_offset**2, 0.5))

                # Put gaussian scores into an array
                array_scores = np.array(
                    [gaussIntScore, gaussMZScore, gaussRTScore])

                # Calculate score by finding the magnitude of vector between 0
                # (centre of 3d cube) and 3**0.5 (distance from centre to
                # corner of 3d cube), the higher the number, the closer to the
                # target
                candidate_score = old_div(
                    np.linalg.norm(array_scores), (3**0.5))

                # Keep best scoring candidate hit
                if candidate_max_score < candidate_score:
                    candidate_max_score = candidate_score

            # Add best candidate score to running total of score
            score += candidate_max_score

        # Return final score and number of lipid like features found
        return old_div(score, float(self.target_len)), len(llf_index)

    def __proc_LF__(self, parameters):

        # Change parameters ready for new run
        self.lfd.peakAdjacentFrameMaxRT = parameters[0]
        self.lfd.peakMaxRTWidth = parameters[1]
        self.lfd.peakMinFoldCutOff = parameters[2]
        self.lfd.mzSizeErrorPPM = parameters[3]

        # Perform LipidFinder clustering (saveFrame set to False to increase
        # speed)
        clus.mzCluster(self.lfd, False)
        clus.featureCluster(self.lfd, False)

        # Perform peak finding (saveFrame set to False to increase speed)
        pf.processAllFeatures(self.lfd, False)

    # Method to generate the full 1opt set of parameters (ps_1opt) using the
    # parameter index
    def __get_ps_1opt__(self, index):

        # List to store 1opt parameter sets produced for parameter at index
        ps_1opt = []

        # Use of list to get a local copy of self.current_parameters
        current_parameters = list(self.current_parameters)

        for parameter in self.all_parameters[index]:

            # Ensure ps will be legal according to rule regarding parameter 0
            # and 1, if not do not append to ps_1opt
            if index == 0 and current_parameters[1] < 3 * parameter:
                continue
            elif index == 1 and parameter < 3 * current_parameters[0]:
                continue

            # Change parameter at index position to current parameter
            current_parameters[index] = parameter

            # If parameter set has not been scored before append to list, if it
            # has do not add
            if tuple(current_parameters) not in self.visited:
                ps_1opt.append(list(current_parameters))
        return ps_1opt

    # Append the results to a pandas DataFrame for ultimate export to csv
    def __append_results__(self):

        for i in self.visited:
            temp = pd.DataFrame({'Parameters': [i], 'Hill climb iteration': self.visited[i][0], 'Parameter set cycle': self.visited[i][
                                1], 'Parameter change': self.visited[i][2], 'Score': self.visited[i][3], 'Lipid like features found': self.visited[i][4], 'Time': self.visited[i][5]})
            self.results = pd.concat([self.results, temp])

    # Export the results as a csv file
    def __export_results__(self, export_file_name=None):
        if export_file_name:
            self.export_file_name = export_file_name
        self.results.to_csv(self.export_file_name +
                            '.csv', sep=',', index=False)
