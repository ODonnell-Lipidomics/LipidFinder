# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Python 3 compatibility tools.

The inclusion of this module provides full support for this tools in
Python 2.7 and Python 3.3 (or newer).

_py3k.__init__ file from MEvoLib [J. Alvarez-Jarreta et al. 2017] has
been used as template to create this file.
"""

import os
import sys


if (sys.version_info[0] >= 3):
    # Code for Python 3
    from builtins import open, zip, map, filter, range, input, round
    import codecs, io

    def _is_int_or_long(value):
        """Return True if 'value' is an integer, False otherwise."""
        # There are no longs on Python 3
        return isinstance(value, int)

    def viewkeys(dictionary):
        """Return a view of the keys of 'dictionary'."""
        return dictionary.keys()

    def viewvalues(dictionary):
        """Return a view of the values of 'dictionary'."""
        return dictionary.values()

    def viewitems(dictionary):
        """Return a view of the items of 'dictionary'."""
        return dictionary.items()

    # On Python 3 urllib, urllib2, and urlparse were merged
    from urllib.request import urlopen, Request, urlretrieve, urlparse
    from urllib.parse import urlencode, quote, quote_plus
    from urllib.error import HTTPError, URLError
    # On Python 3 subprocess.DEVNULL exists
    from subprocess import DEVNULL
    #On Python 3, this will be a unicode StringIO
    from io import StringIO
    from tempfile import TemporaryDirectory
else:
    # Code for Python 2
    from __builtin__ import open, basestring, unicode, round
    # Import Python 3 like iterator functions:
    from future_builtins import zip, map, filter
    from __builtin__ import xrange as range
    from __builtin__ import raw_input as input

    def _is_int_or_long(value):
        """Return True if 'value' is an integer or long, False
        otherwise.
        """
        return isinstance(value, (int, long))

    def viewkeys(dictionary):
        """Return a view of the keys of 'dictionary'."""
        return ( dictionary.viewkeys() )

    def viewvalues(dictionary):
        """Return a view of the values of 'dictionary'."""
        return dictionary.viewvalues()

    def viewitems(dictionary):
        """Return a view of the items of 'dictionary'."""
        return dictionary.viewitems()

    # Under urllib.request on Python 3:
    from urllib2 import urlopen, Request
    from urllib import urlretrieve
    from urlparse import urlparse
    # Under urllib.parse on Python 3:
    from urllib import urlencode, quote, quote_plus
    # Under urllib.error on Python 3:
    from urllib2 import HTTPError, URLError
    # On Python 2 subprocess.DEVNULL doesn't exist
    DEVNULL = open(os.path.devnull, 'w')
    # On Python 2 this will be a (bytes) string based handle.
    # Note this doesn't work as it is unicode based:
    #     from io import StringIO
    try:
        from cStringIO import StringIO
    except ImportError:
        from StringIO import StringIO
    try:
        input = raw_input
    except NameError:
        pass


if (sys.platform == "win32"):
    # Can't use commands.getoutput on Python 2, Unix only/broken:
    # http://bugs.python.org/issue15073
    # Can't use subprocess.getoutput on Python 3, Unix only/broken:
    # http://bugs.python.org/issue10197
    def getoutput(cmd):
        import subprocess
        child = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 universal_newlines=True, shell=False)
        stdout, stderr = child.communicate()
        # Remove trailing "\n" to match the Unix function
        return stdout.rstrip("\n")
elif (sys.version_info[0] >= 3):
    # Use subprocess.getoutput on Python 3
    from subprocess import getoutput
else:
    # Use commands.getoutput on Python 2
    from commands import getoutput
