# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to summarise the putative profile:
    > create_summary():
        Create a summary XLSX file containing only one row per m/z and
        retention time with the most common lipid category.

Examples:
    >>> from Configuration import LFParameters
    >>> from LFDataFrame import LFDataFrame
    >>> from MSSearch import Summary
    >>> parameters = LFParameters('mssearch', 'parameters.json')
    >>> data = LFDataFrame('profile.csv', parameters)
    >>> Summary.create_summary(data, parameters)

    It can also work with a pandas.DataFrame as input:
    >>> import pandas
    >>> from Configuration import LFParameters
    >>> from MSSearch import Summary
    >>> parameters = LFParameters('mssearch', 'parameters.json')
    >>> data = pandas.read_csv('profile.csv')
    >>> Summary.create_summary(data, parameters)
"""

import os

import pandas


def create_summary(data, parameters, dst=''):
    # type: (object, LFParameters, str) -> None
    """Create a summary XLSX file containing only one row per m/z and
    retention time with the most common lipid category.

    'data' must have, at least, m/z, retention time (RT), "Main Class"
    and "Category" columns.
    If two or more categories have the same number of matches, the first
    one to appear in 'data' is chosen, unless it is "other metabolites":
    in that case the second one will be selected. The same rule applies
    for the main class.
    If 'dst' is not an absolute path, the current working directory will
    be used as starting point. If "mssearch_<db>_summary.xlsx"
    file already exists, it will be overwritten without warning.
    "<db>" stands for the selected LIPID MAPS database.

    Keyword Arguments:
        data       -- LFDataFrame or pandas.DataFrame instance
        parameters -- LipidFinder's MS Search parameters instance
        dst        -- destination directory where the XSLX file will be
                      created [default: current working directory]
    """
    mzCol = parameters['mzCol']
    rtCol = parameters['rtCol']
    # Get the list of unique m/z values from 'data'
    mzList = data[mzCol].unique().tolist()
    summary = pandas.DataFrame(columns=list(data))
    # Iterate over the groups of identifications for the same m/z and RT
    for name, group in data.groupby([mzCol, rtCol]):
        categories = group['Category'].value_counts()
        if (categories.empty):
            # No matches were returned for this group
            summary = summary.append(group, ignore_index=True)
        else:
            bestCategory = categories.index[0]
            # Keep only those rows for the most frequent category
            subgroup = group[group['Category'] == bestCategory]
            if (bestCategory == 'other metabolites'):
                if ((len(categories) == 1) or (categories[0] > categories[1])):
                    # "Other metabolites" has empty "Main Class", but it
                    # can have more than one row, so keep the first one
                    summary = summary.append(subgroup.head(1),
                                             ignore_index=True)
                    continue
                else:
                    # Change the best category for a specialized one
                    # with the same number of matches (next category)
                    bestCategory = categories.index[1]
                    subgroup = group[group['Category'] == bestCategory]
            # Get the most frequent main class and keep its first row
            mainClass = subgroup['Main Class'].value_counts().index[0]
            summary = summary.append(
                    subgroup[subgroup['Main Class'] == mainClass].head(1),
                    ignore_index=True)
    # Create XLSX file with the summary putative profiling in 'dst'
    fileName = 'mssearch_{0}_summary.xlsx'.format(
            parameters['database'].lower())
    summary.to_excel(os.path.join(dst, fileName), index=False,
                     engine='xlsxwriter')
