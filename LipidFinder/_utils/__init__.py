# -*- coding: utf-8 -*-
# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Utility functions for LipidFinder."""

from __future__ import print_function

import os


def normalise_path(path):
    # type: (str) -> str
    """Return a valid path for the current OS.

    Keyword Arguments:
        path -- path to normalize
    """
    return os.path.normpath(os.path.realpath(path))


def traceless_warning(message, category, filename, lineno, file=None,
                      line=None):
    # type: (str, Warning, str, int, file, str) -> str
    """Return a warning message without the traceback information.

    Keyword Arguments:
        message  -- warning message
        category -- Warning instance
        filename -- name of the file where the warning was raised
        lineno   -- line number where the warning was raised
        file     -- file instance
        line     -- line object
    """
    return 'Warning{0}{1}{1}'.format(message, os.linesep)


def mz_delta(mz, fixederr, ppmerr, precision=5):
    # type: (float, float, float, int) -> float
    """Return the delta tolerance for the given m/z.

    Keyword Arguments:
        mz        -- m/z reference value
        fixederr  -- allowed fixed error
        ppmerr    -- mass-dependant PPM error to add to the fixed error
        precision -- number of decimal digits to use with floats (e.g. a
                     precision of 2 forces a difference of 0.01 between
                     two any consecutive float numbers) [default: 5]
    """
    return round(fixederr + (mz * ppmerr * 1e-6), precision)


def mz_tol_range(mz, fixederr, ppmerr, precision=5):
    # type: (float, float, float, int) -> float
    """Return lower and upper tolerance limits for the given m/z.

    Keyword Arguments:
        mz        -- m/z reference value
        fixederr  -- allowed fixed error
        ppmerr    -- mass-dependant PPM error to add to the fixed error
        precision -- number of decimal digits to use with floats (e.g. a
                     precision of 2 forces a difference of 0.01 between
                     two any consecutive float numbers) [default: 5]
    """
    delta = mz_delta(mz, fixederr, ppmerr, precision)
    return (round(mz - delta, precision), round(mz + delta, precision))


def rt_delta(maxdiff, precision=5):
    # type: (float, int) -> float
    """Return the delta tolerance for the given retention time.

    Keyword Arguments:
        maxdiff   -- maximum time difference between a feature edge and
                     an adjacent frame to be considered part of the same
                     feature
        precision -- number of decimal digits to use with floats (e.g. a
                     precision of 2 forces a difference of 0.01 between
                     any two consecutive float numbers) [default: 5]
    """
    return round(maxdiff, precision)


def rt_tol_range(rt, maxdiff, precision=5):
    # type: (float, float, int) -> float
    """Return lower and upper tolerance limits for the given retention
    time.

    Keyword Arguments:
        rt        -- retention time (RT) reference value
        maxdiff   -- maximum time difference between a feature edge and
                     an adjacent frame to be considered part of the same
                     feature
        precision -- number of decimal digits to use with floats (e.g. a
                     precision of 2 forces a difference of 0.01 between
                     any two consecutive float numbers) [default: 5]
    """
    delta = rt_delta(maxdiff, precision)
    return (round(rt - delta, precision), round(rt + delta, precision))


def print_progress_bar(iteration, total, prefix='', suffix='Completed',
                       length=34):
    # type: (int, int, str, str, int) -> None
    """Call in a loop to create terminal progress bar.

    Extracted from first answer at: https://stackoverflow.com/questions/
    3173320/text-progress-bar-in-the-console

    Keyword Arguments:
        iteration -- current iteration
        total     -- total iterations
        prefix    -- prefix string [default: no prefix]
        suffix    -- suffix string [default: "Completed"]
        length    -- character length of bar [default: 34]
    """
    percent = "{0:.1f}".format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = '#' * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if (iteration == total):
        print()
