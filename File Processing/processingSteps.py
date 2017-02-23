"""This module is used to get the number of input files, filenames, original websearch 
filename, read all files into separate dataframes and then combine the results dataframes
into one overall dataframe that will be processed further
"""
from __future__ import print_function
from builtins import str
from builtins import input
from builtins import range
import pandas as pd
import time
import os


def numFiles():
    """This is to get the number of database results input files that need to be processed

    Returns:
        int: number of input files
    """
    while True:
        try:
            countFiles = int(input("Number of input files:-\n\n    "))
            return countFiles  # need to make sure have more than one file!!
        except ValueError:
            print ("Integer needed, try again!")


def inputFiles(numFiles):
    """Asks used the input the name of each of the database results files
    separately, reads them into separate dataframes and then concatenates them
    into one new dataframe

    Args:
        number of input files (int): the number of database results files 
        that will be read into dataframes

    Returns:
        dataframe: dataframe of oncatenated database results dataframe
    """
    d = {}
    inputFile = {}

    for i in range(1, numFiles + 1):
        d[i] = "inputFile{0}".format(i)
        inputFile[i] = d[i]

    for i in range(1, numFiles + 1):
        inputF = input("Name of file " + str(i) + ":-\n\n    ")
        inputFile[i] = pd.read_table(inputF, sep=',')

# default is to concatenate along axis=0, i.e. rows - this is what I want
# ignore index, as want to generate a new one based on concatenated data
    dfCombined = pd.concat(inputFile, ignore_index=True).fillna("NA")

    return (dfCombined)


def inputOriginalfile():
    """Asks user to input the name of the original websearch input file 
    and reads it into a dataframe

    Returns:
        dataframe, str: original websearch file as dataframe, original websearch 
        filename as string
    """
    originalFile = input(
        "Name of original masses/retention times file:-\n\n    ")
    originalDF = pd.read_table(originalFile, sep=',')
    # rename MZ and Time from the lipidFinder output
    originalDF.rename(columns={'MZ': 'ORIGINAL_MASS'}, inplace=True)
    originalDF.rename(columns={'Time': 'RETENTION_TIME'}, inplace=True)

    return originalDF, originalFile


def readParameters():
    """Read the parameters file ad extract parameter(s) used during processing

    Returns:
        float: delta ppm value used as a cutoff during processing
    """
    parmsFile = pd.read_table('websearch_parameters.csv', sep=',')

    deltaPPM = int((parmsFile.ix[3][1]))

    return (deltaPPM)
