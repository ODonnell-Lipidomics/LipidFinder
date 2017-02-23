"""This module filters the websearch results based on delta ppm value
adducts; assigns lipid classes; removes duplicate records and 
generates output results files
"""
from builtins import zip
import pandas as pd
import time
import os
import numpy as np

# only selects all records with delta_ppm values between value specified
# in parameters file


def deltaElim(mergeDF, deltaPPM):
    """Filters results based on the delta ppm value specified by the user
    in the websearch_parameters.csv file. Records with values above and below
    this delta ppm value are removed from the dataset

    Args:
        mergeDF (dataframe): input dataframe
        deltaPPM (float): delta ppm cutoff value specified by user 
        in websearch_parameters.csv

    Returns:
        dataframe: filtered dataframe
    """
    mergeDF = mergeDF[(mergeDF['DELTA_PPM'] <= deltaPPM) &
                      (mergeDF['DELTA_PPM'] >= -deltaPPM)]

    return mergeDF


# rename all categories listed in search results file to the standard
# categories from the category map file
def adductRename(mergeDF):
    """Adduct names are renamed to the same standard names (in old_new_map below)

    Args:
        mergeDF (dataframe): input dataframe

    Returns:
        dataframe: output dataframe
    """
    old_new_map = {'M-H': '[M-H]-', '[M-H]-': '[M-H]-', 'M+CH3COO': '[M+CH3COO]-', '[M+CH3COO]-': '[M+CH3COO]-', 'M+H': '[M+H]+',
                   '[M+H]+': '[M+H]+', 'M+Na': '[M+Na]+', '[M+Na]+': '[M+Na]+', 'M+NH4': '[M+NH4]+', '[M+NH4]+': '[M+NH4]+', 'NA': 'NA'}
    mergeDF['ADDUCT_ION'] = mergeDF['ADDUCT_ION'].map(old_new_map)

    return mergeDF


# only keeps the adducts listed in adductList and deletes all others
def adductElim(mergeDF):
    """Dataframe is filtered so that we only keep the adducts listed in adductList below

    Args:
        mergeDF (dataframe): input dataframe

    Returns:
        dataframe: output dataframe
    """
    adductList = ['[M-H]-', '[M+H]+', '[M+Na]+',
                  '[M+NH4]+', '[M+CH3COO]-', 'NA']
    mergeDFNew = mergeDF[(mergeDF['ADDUCT_ION'].isin(adductList))]

    return mergeDFNew


# rename all categories listed in search results file to the standard
# categories from the category map file
def categoryRename(mergeDF):
    """Lipid categories are renamed to the standard lipid category names
    as per LIPIDMAPS. The categories_map.csv is used to 'map' old category
    names to new category names

    Args:
        mergeDF (dataframe): input dataframe

    Returns:
        dataframe: output dataframe
    """
    categoryFileDF = pd.read_table(
        "categories_map.csv", sep=',', keep_default_na=False)

    # new way: make the categories df into a dictionary - much faster!!
    catMap = dict(list(zip(categoryFileDF.old_category,
                           categoryFileDF.new_category)))
    mergeDF['CATEGORY'] = mergeDF['CATEGORY'].map(catMap)

    return mergeDF


def removeDuplicateNames(mergeDF):
    """Records that match on mass, retention time and name are found:
    of each of the matches, only the one with the lowest delta ppm value is retained
    and the other are removed. Output file is generated of the overall merged results
    from all databases

    Args:
        mergeDF (dataframe): input dataframe

    Returns:
        dataframe: output dataframe
    """
    # compare on names reduced to lower case, use one with lower delta ppm
    groupedDF = mergeDF.groupby([mergeDF.ORIGINAL_MASS, mergeDF.RETENTION_TIME, mergeDF.NAME.str.lower()])[
        'DELTA_PPM'].agg(lambda col: col.idxmin())

    # Get the full record for each of the results in b
    sortedDF = mergeDF.ix[groupedDF]

    # Sort by mass and then name
    sortedDF.sort_values(
        by=['ORIGINAL_MASS', 'RETENTION_TIME', 'NAME'], inplace=True)

    dt = time.strftime("%Y%m%d-%H%M%S")
    sortedFile = "Merged_processed_file_" + dt + ".csv"
    sortedDF.to_csv(sortedFile, sep=',', float_format='%.6f', index=False)

    return sortedDF, sortedFile


def getMasscategorycounts(finalDF, originalDF):
    """works out the most commonly found lipid category per mass in original websearch 
    file and how many of that category were found

    Args:
        finalDF (dataframe): input dataframe
        originalDF (dataframe): original websearch file stored as a dataframe

    Returns:
        dataframe: output dataframe
    """
    subsetMerged = pd.DataFrame({'Count': finalDF.groupby(
        ['ORIGINAL_MASS', 'RETENTION_TIME', 'CATEGORY']).size()}).reset_index()

    groupedDF = subsetMerged.groupby(['ORIGINAL_MASS', 'RETENTION_TIME'])

    # set up a new DF for category counts results
    outFrame = pd.DataFrame()
    # group gives groupedDF as mass/rt groups (incl. category and count)
    for name, group in groupedDF:
        # if there is only one count in each group, then use this and write it
        # to the new DF
        if len(group) == 1:
            df = group
            outFrame = pd.concat([outFrame, df])
        else:
            # there is more than one count in a group, so pick the one with the maximum count value (for category),
            # but ignore 'other metabolites' cases - we don't want to record
            # this category unless its the only one
            df = group.ix[group['Count'][group['CATEGORY'] !=
                                         'other metabolites'].idxmax()].to_frame().transpose()
            outFrame = pd.concat([outFrame, df])

    # as only have some rows from original DF, we need to reset the index and
    # then re-sort by mass and rt
    outFrame.reset_index(inplace=True, drop=True)
    outFrame = outFrame.sort_values(
        by=['ORIGINAL_MASS', 'RETENTION_TIME'], ascending=[True, True])

    # concatenate the above results with the original input masses (this one has no categories)
    # we merge with the original dataframe that includes duplicates (multiple retention times per mass value) - one that was used for websearch
    # uses a left join, so returns all the rows from the original dataframe
    # and includes matches from outFrame dataframe
    finalFile = originalDF.merge(
        outFrame, on=['ORIGINAL_MASS', 'RETENTION_TIME'], how='left')
    finalFile.reset_index(inplace=True, drop=True)

    finalFile = finalFile.sort_values(
        by=['ORIGINAL_MASS', 'RETENTION_TIME'], ascending=[True, True])

    # replace 'NA' with 'unknown' in Category columns and 'NA' with '1' in
    # Count columns
    finalFile['CATEGORY'].fillna('unknown', inplace=True)
    finalFile['Count'].fillna(1, inplace=True)

    dt = time.strftime("%Y%m%d-%H%M%S")
    finalFile.to_csv("Category_counts_" + dt + ".csv",
                     sep=',', float_format='%.6f', index=False)

    return finalFile


def getDatabasecategorycounts(originalDF, finalDF):
    """Generates 2 summary files: one which shows the number of results found per mass
    per database and the other which shows a presence/absence of results found per mass 
    per database

    Args:
        originalDF (dataframe): original websearch file stored as a dataframe
        finalDF (dataframe): input dataframe

    Returns:
        dataframe: output dataframe
    """
    # get count per database per mass
    newDF = pd.DataFrame({'count': finalDF.groupby(
        ['ORIGINAL_MASS', 'DATABASE']).size()}).reset_index()

    # reshape the data so that that 1st column is mass and then each database
    # is a column
    newDF2 = newDF.pivot(index='ORIGINAL_MASS', columns='DATABASE')

    # reset the index after reshaping
    resultDF = newDF2.reset_index()

    # fill empty cells with 0s
    resultDF.fillna('0', inplace=True)
    resultDF2 = newDF2.reset_index()

    # check which databases we have in the file
    hmdbFlag = 0
    lipidhomeFlag = 0
    lipidmapsFlag = 0
    metlinFlag = 0

    # if we have the database in the file, set flag to 1
    newColumns = ['ORIGINAL_MASS']
    if any(finalDF.DATABASE == 'HMDB'):
        hmdbFlag = 1
        newColumns.append('HMDB')
    if any(finalDF.DATABASE == 'LIPIDHOME'):
        lipidhomeFlag = 1
        newColumns.append('LIPIDHOME')
    if any(finalDF.DATABASE == 'LIPIDMAPS'):
        lipidmapsFlag = 1
        newColumns.append('LIPIDMAPS')
    if any(finalDF.DATABASE == 'METLIN'):
        metlinFlag = 1
        newColumns.append('METLIN')

    # get rid of the multi-index and rename columns
    resultDF.columns = resultDF.columns.get_level_values(0)
    resultDF.columns = [''.join(col).strip()
                        for col in resultDF.columns.values]
    resultDF.columns = newColumns

    resultDF2.columns = resultDF2.columns.get_level_values(0)
    resultDF2.columns = [''.join(col).strip()
                         for col in resultDF2.columns.values]
    resultDF2.columns = newColumns

    finalFile = pd.merge(originalDF, resultDF, on=[
                         'ORIGINAL_MASS'], how='left')
    finalFile.fillna('0', inplace=True)
    finalFile.reset_index(inplace=True, drop=True)

    finalFile = finalFile.sort_values(
        by=['ORIGINAL_MASS', 'RETENTION_TIME'], ascending=[True, True])

    dt = time.strftime("%Y%m%d-%H%M%S")
    finalFile.to_csv("DB_mass_counts_" + dt + ".csv", sep=',',
                     float_format='%.6f', index=False)

    if hmdbFlag == 1:
        resultDF2.ix[(resultDF2.HMDB > 0), 'HMDB'] = 1
    if lipidhomeFlag == 1:
        resultDF2.ix[(resultDF2.LIPIDHOME > 0), 'LIPIDHOME'] = 1
    if lipidmapsFlag == 1:
        resultDF2.ix[(resultDF2.LIPIDMAPS > 0), 'LIPIDMAPS'] = 1
    if metlinFlag == 1:
        resultDF2.ix[(resultDF2.METLIN > 0), 'METLIN'] = 1

    # fill empty cells with 0s
    resultDF2.fillna('0', inplace=True)

    finalFile2 = pd.merge(originalDF, resultDF2, on=[
                          'ORIGINAL_MASS'], how='left')
    finalFile2.fillna('0', inplace=True)
    finalFile2.reset_index(inplace=True, drop=True)

    finalFile2 = finalFile2.sort_values(
        by=['ORIGINAL_MASS', 'RETENTION_TIME'], ascending=[True, True])

    dt = time.strftime("%Y%m%d-%H%M%S")
    finalFile2.to_csv("DB_mass_presence_absence_" + dt + ".csv",
                      sep=',', float_format='%.6f', index=False)

    return resultDF, resultDF2
