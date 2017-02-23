"""This module gets the input data and websearch parameters
"""
from builtins import input
import pandas as pd


def readFile():
    """This function takes a file (.csv format) input by the user and reads it into a dataframe

    Returns:
        dataframe: dataframe containing the data from the input file
    """
    inputFile = pd.read_table(input('Name of file:-\n\n    '), sep=',')

    return (inputFile)


def readParameters():
    """This function loads the websearch_parameters.csv file into a dataframe and extracts the 
    parameter values and assigns them to variables

    Returns:
        float: m/z tolerance for LIPIDMAPS database search
        float: m/z tolerance for LipidHome database search
        float: m/z tolerance for HMDB database search
        str: database search option ("COM", "CUR" or "ALL")
    """
    parmsFile = pd.read_table('websearch_parameters.csv', sep=',')

    lipidmapsTol = float((parmsFile.ix[0][1]))
    lipidhomeTol = float((parmsFile.ix[1][1]))
    hmdbTol = float((parmsFile.ix[2][1]))
    searchCriteria = parmsFile.ix[4][1]

    return (lipidmapsTol, lipidhomeTol, hmdbTol, searchCriteria)
