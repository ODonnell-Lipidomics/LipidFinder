# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Set of methods aimed to generate data plotting figures:
    > category_scatterplot():
        Generates a lipid category scatter plot from the input dataframe
        and exports it in the selected format (PDF by default).

Examples:
    >>> df = {'MZ': pd.Series([101.23, 403.78]),
    ...       'RT': pd.Series([12, 31]),
    ...       'Category': pd.Series(['PS', 'PC'])}
    >>> category_scatterplot(df, 'MZ', 'RT')

    The exported file with the scatter plot figure will be located at
    the current working directory under the name "category_scatterplot"
    with the extension of the chosen format. Each lipid category will
    always have the same color assigned. The color palette can be
    changed to color-blind friendly.
"""

import os
import string
import warnings

from cycler import cycler
import matplotlib.style
import matplotlib.pyplot as pyplot
import numpy
import pandas


# Ignore Future Warnings from pandas library
warnings.simplefilter(action='ignore', category=FutureWarning)
# List of lipid categories in reverse order (including "unknown")
CATEGORIES = ['unknown', 'sterol lipids', 'sphingolipids', 'saccharolipids',
              'prenol lipids', 'polyketides', 'other metabolites',
              'glycerophospholipids', 'glycerolipids', 'fatty acyls']


def category_scatterplot(data, parameters, dst):
    # type: (object, LFParameters, str) -> None
    """Create and export the lipid category scatter plot.

    Generates a scatter plot of the main lipid categories based on m/z
    vs retention time (RT) of the input dataframe, which must contain a
    column named "Category" (case sensitive). The plot is saved in the
    file format selected during the parameter configuration (PDF by
    default).

    Keyword arguments:
        data       -- LFDataFrame or pandas.DataFrame instance
        parameters -- LipidFinder's MS Search parameters instance
        dst        -- destination directory where the figure will be
                      exported
    """
    mzCol = parameters['mzCol']
    rtCol = parameters['rtCol']
    # Choose the color palette to assign to each lipid category
    if (parameters['figColors'] == 'standard'):
        colors = ['#111111', '#ffC0cb', '#1e90ff', '#00ffff', '#ffd700',
                  '#ff1493', '#ff8c00', '#32cd32', '#ff0000', '#9370db']
        markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
        sizes = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
        widths = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    elif (parameters['figColors'] == 'colorblind'):
        colors = ['#b6dbff', '#009292', '#490092', '#ff6db6', '#004bec',
                  '#111111', '#24ff24', '#db6d00', '#920000', '#ffff6d']
        markers = ['o', 'v', '^', '<', '>', 's', 'P', 'X', '*', 'd']
        sizes = [3, 4, 4, 4, 4, 3, 4, 4, 5, 4]
        widths = [0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0.1]
    # Replace NaN in "Category" column by "unknown"
    data['Category'].fillna('unknown', inplace=True)
    # Get the main category for each m/z and RT pair
    catData = _get_main_categories(data, mzCol, rtCol)
    # Configure the plot style, layout and parameters
    matplotlib.style.use('seaborn-paper')
    ax = pyplot.subplot(111)
    # Calculate the ceiling of the 102% of the maximum retention time
    maxRT = numpy.amax(catData[rtCol].values)
    maxX = numpy.ceil(1.02 * maxRT)
    # Calculate the ceiling of the 110% of the maximum m/z
    maxMZ = numpy.amax(catData[mzCol].values)
    maxY = numpy.ceil(1.1 * maxMZ)
    # Set range of X and Y axes
    pyplot.xlim([0, maxX])
    pyplot.ylim([0, maxY])
    # Set label text of X and Y axes
    pyplot.xlabel('Retention time (min)')
    pyplot.ylabel('m/z', fontstyle='italic')
    # Set a readable and color-blind friendly marker cycle
    ax.set_prop_cycle(
            cycler('marker', markers) + cycler('color', colors) +\
            cycler('markersize', sizes) + cycler('markeredgewidth', widths))
    # Load each category to the plot
    for i, category in enumerate(CATEGORIES):
        catMatches = catData.loc[catData['Category'].str.lower() == category]
        if (len(catMatches) == 0):
            # Skip to the next marker so each category will always have
            # the same one, allowing an ease comparison between plots
            next(ax._get_lines.prop_cycler)
        else:
            pyplot.plot(catMatches[rtCol], catMatches[mzCol], linestyle='None',
                        markeredgecolor='#666666',
                        label=string.capwords(category))
    # Get handles and labels for legend
    handles, labels = ax.get_legend_handles_labels()
    numCats = len(labels)
    # Change the axes position of the plot to leave some room under it
    # to display the legend based on its number of rows for an up to
    # 5-column layout
    box = ax.get_position()
    numRows = numpy.ceil(numCats / 5.0)
    yShift = numRows * 0.05
    ax.set_position([box.x0, box.y0 + box.height * yShift, box.width,
                     box.height * (1.0 - yShift)])
    # Reverse labels in legend to sort them alphabetically, locate the
    # legend below the plot and display the categories in 1 or 2 rows
    # and up to 5 columns
    numCols = int(numpy.ceil(numCats / numRows))
    ax.legend(handles[::-1], labels[::-1], loc='upper center',
              bbox_to_anchor=(0.5, -0.12), fancybox=False, shadow=False,
              ncol=numCols, numpoints=1)
    # Adjust plot parameters to fit it into the figure area
    pyplot.tight_layout()
    # Save the figure into the selected file format and close it
    figName = 'category_scatterplot_{0}.{1}'.format(
            parameters['database'].lower(), parameters['figFormat'])
    figPath = os.path.join(dst, figName)
    pyplot.savefig(figPath, dpi=600, bbox_inches='tight')
    pyplot.close()


def _get_main_categories(data, mzCol, rtCol, defaultCat='other metabolites'):
    # type: (object, str, str, str) -> pandas.DataFrame
    """Return a dataframe with the most frequent lipid category per m/z
    and retention time.

    The returned dataframe contains the most frequent lipid category per
    m/z and retention time (RT) in the input dataframe. If there is more
    than one category per m/z and RT tuple, the default category is
    excluded.

    Keyword arguments:
        data       -- LFDataFrame or pandas.DataFrame instance
        mzCol      -- name of m/z column in 'data'
        rtCol      -- name of retention time column in 'data'
        defaultCat -- default category [default: "other metabolites"]
    """
    # Check the default category is within the expected values
    if (defaultCat not in CATEGORIES):
        raise ValueError("'defaultCat' must be one of {0}".format(CATEGORIES))
    # Get count the number of matches per m/z, RT and category
    categoryCounts = pandas.DataFrame(
                {'Count': data.groupby([mzCol, rtCol, 'Category'],
                                       sort=True).size()})
    categoryCounts.reset_index(inplace=True)
    # Group the category counts by m/z and RT
    groupedData = categoryCounts.groupby([mzCol, rtCol])
    # Set up a new dataframe to save the most frequent categories
    catData = pandas.DataFrame()
    for name, group in groupedData:
        # Append the row with the most frequent category for each m/z
        # and RT (excluding "Count" column)
        if (len(group) == 1):
            catData = catData.append(group[[mzCol, rtCol, 'Category']],
                                     ignore_index=True)
        else:
            # When there are two or more categories, pick the
            # non-default most frequent one
            index = group.loc[group['Category'].str.lower() != defaultCat,
                              'Count'].idxmax()
            catData = catData.append(
                    group.loc[index, [mzCol, rtCol, 'Category']],
                    ignore_index=True)
    return catData
