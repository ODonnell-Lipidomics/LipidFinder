# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Collection of parameters used to set up each LipidFinder module."""

from collections import OrderedDict
import json
import os
import warnings

import pandas
import pkg_resources

from LipidFinder._py3k import viewitems
from LipidFinder._utils import normalise_path, traceless_warning


# List of LipidFinder modules
MODULES = ['peakfilter', 'amalgamator', 'mssearch']


class LFParameters:
    """A LFParameters object stores a set of parameters to be used by
    the given LipidFinder module.

    Attributes:
        _parameters  (Private[collections.OrderedDict])
            Dictionary where the parameters and their associated
            information are stored.
        _module  (Private[str])
            LipidFinder module: "peakfilter", "amalgamator" or
            "mssearch".

    Examples:
        LFParameters objects can be created in two different ways:
            >>> from Configuration import LFParameters
            >>> default = LFParameters('peakfilter')
            >>> params = LFParameters('peakfilter', 'parameters.json')

        The former will load the default value for each PeakFilter
        parameter whilst the latter will override each PeakFilter
        parameter's default value by the new value found in the JSON
        parameters file "parameters.json". If one of the values in the
        JSON file is not valid, a warning message will be displayed
        indicating what is wrong with that parameter and the default
        value will not be overwritten. The new parameters file can have
        all or only part of the whole set of parameters of the module.

        Once the parameter set is loaded, the user can access each
        parameter through dictionary-like features:
            >>> 'polarity' in default
            True
            >>> default['polarity']
            'Negative'

        Or modify each parameter's value with a new valid value:
            >>> default['polarity'] = 'Positive'
            >>> default['polarity']
            'Positive'

        If the new value is not valid, a warning message will we printed
        and the parameter's value will not be overwritten:
            >>> default['polarity'] = 'Neutral'
            Warning (polarity): 'Neutral' not in: Positive, Negative
            >>> default['polarity']
            'Positive'

        Finally, the edited parameter set can be saved to a JSON file:
            >>> default.write()
            >>> default.write('/home/user/new_paramaters.json')

        Only active parameters (those that will be used based on the
        specified module) will be saved in the file. The former call
        will save the parameters to "parameters.json" file in the
        current working directory. The latter will save them to
        "/home/user/new_parameters.json". If the destination file
        already exists, it will be overwritten without warning.
    """

    def __init__(self, module='peakfilter', src=''):
        # type: (str, str) -> LFParameters
        """Constructor of the class LFParameters.

        First, the template file is loaded and the module's parameters
        selected. Next, if a source JSON parameters file path is
        provided, the default values are overwritten by the
        corresponding new (valid) values.

        Keyword Arguments:
            module -- LipidFinder module [default: "peakfilter"]
            src    -- path of an existing JSON parameters file
                      [default: leave default values]
        """
        # Check if the selected module is included in LipidFinder
        if (module.lower() not in MODULES):
            raise ValueError("'module' must be one of {0}".format(MODULES))
        self._module = module.lower()
        # Show warnings always
        warnings.simplefilter('always')
        # Replace the warning message handler to remove traceback information
        warnings.formatwarning = traceless_warning
        # Load the parameters template that contains all the parameters
        # used by all modules, including their description, restrictions
        # and default values
        templatePath = pkg_resources.resource_filename(
                'LipidFinder', 'Data/parameters_template.json')
        with open(templatePath, 'r') as templateFile:
            self._parameters = json.load(templateFile,
                                         object_pairs_hook=OrderedDict)
        # Keep only the parameters of the selected module
        self._parameters = OrderedDict(
                [(k, v) for k, v in viewitems(self._parameters)
                        if module in v['modules']])
        # Set each parameter value to its default, or to the minimum/
        # first option if there is no default entry in the template
        for key, data in viewitems(self._parameters):
            if (data['type'] == 'bool'):
                data['value'] = data.get('default', False)
            elif (data['type'] == 'path'):
                if ('default' in data):
                    # Base the default path value on LipidFinder's root
                    # folder instead of the current working directory
                    data['value'] = normalise_path(
                            pkg_resources.resource_filename('LipidFinder',
                                                            data['default']))
                else:
                    data['value'] = None
            elif (data['type'] in ['str', 'int', 'float', 'selection']):
                data['value'] = data.get('default', None)
            else:
                data['value'] = data.get('default', '[]')
        if (src):
            # Load the parameter values of an existing JSON parameters
            # file
            with open(src, 'r') as srcFile:
                newParams = json.load(srcFile)
            # Check if the value is valid before overwriting the default
            # one to ensure the coherence between parameters
            for key in (x for x in self._parameters.keys() if x in newParams):
                self[key] = newParams[key]

    def __getitem__(self, key):
        # type: (str) -> object
        """Return the value of the parameter 'key'.

        Keyword Arguments:
            key -- name of the parameter
        """
        return self._parameters[key]['value']

    def __setitem__(self, key, value):
        # type: (str, object) -> None
        """Set parameter's new value only if it meets the restrictions.

        The restrictions are: expected type and/or data structure (in
        the case of lists), and value within range/list of options or
        valid path.

        Keyword Arguments:
            key   -- name of the parameter
            value -- new value to be assigned to the parameter
        """
        typeStr = self._parameters[key]['type']
        if (typeStr in ['int', 'float']):
            if (self._validate_number(key, value)):
                self._parameters[key]['value'] = value
        elif (typeStr == 'selection'):
            if (self._validate_selection(key, value)):
                self._parameters[key]['value'] = value
        elif (typeStr == 'path'):
            if (self._validate_path(key, value)):
                self._parameters[key]['value'] = value
        elif (typeStr in ['int range', 'float range']):
            if (self._validate_range(key, value)):
                self._parameters[key]['value'] = value
        elif (typeStr == 'multiselection'):
            if (self._validate_multiselection(key, value)):
                self._parameters[key]['value'] = value
        elif (typeStr == 'pairs'):
            if (self._validate_pairs(key, value)):
                self._parameters[key]['value'] = value
        else: # typeStr in ['str', 'bool']
            if (self._validate_literal(key, value)):
                self._parameters[key]['value'] = value

    def __contains__(self, key):
        # type: (str) -> bool
        """Return True if 'key' is a parameter, False otherwise.

        Keyword Arguments:
            key -- name of the parameter
        """
        return key in self._parameters

    def write(self, dst):
        # type: (str) -> None
        """Write the parameters values into a JSON file.

        The JSON file will contain a dictionary with active parameter
        names as key and their value. If 'dst' is not an absolute path,
        the current working directory will be used as starting point. If
        the destination file already exists, it will be overwritten.

        Keyword Arguments:
            dst -- destiny JSON parameters file path
        """
        # Create a new dictionary of pairs (<parameter_name>: <value>)
        paramsDict = OrderedDict()
        for key, data in self._parameters.items():
            if (self._is_active(key)):
                paramsDict[key] = data['value']
        # Write the dictionary into a JSON parameters file
        with open(normalise_path(dst), 'w') as paramsFile:
            json.dump(paramsDict, paramsFile, indent=4)

    def _validate_literal(self, key, value, verbose=True):
        # type: (str, object, bool) -> bool
        """Return True if 'value' is a valid literal for the chosen
        parameter, False otherwise.

        Keyword Arguments:
            key     -- name of the parameter
            value   -- new value to be assigned to the parameter
            verbose -- print name of the parameter in warning messages
                       [default: True]
        """
        header = ' ({0})'.format(key) if verbose else ''
        typeStr = self._parameters[key]['type']
        if (not isinstance(value, self._eval_type(typeStr))):
            # Condition for string types run in Python 2.7 and Windows
            if ((typeStr == 'str') and (not isinstance(value, basestring))):
                warnings.warn("{0}: '{1}' type expected, '{2}' found".format(
                        header, typeStr, value.__class__.__name__))
                return False
        if (value in [None, '']):
            warnings.warn("{0}: the literal can't be empty".format(header))
            return False
        return True

    def _validate_number(self, key, value, verbose=True):
        # type: (str, object, bool) -> bool
        """Return True if 'value' is a valid int/float for the given
        parameter, False otherwise.

        Keyword Arguments:
            key     -- name of the parameter
            value   -- new value to be assigned to the parameter
            verbose -- print name of the parameter in warning messages
                       [default: True]
        """
        header = ' ({0})'.format(key) if verbose else ''
        typeStr = self._parameters[key]['type']
        if (not isinstance(value, self._eval_type(typeStr))):
            warnings.warn("{0}: '{1}' type expected, '{2}' found".format(
                    header, typeStr, value.__class__.__name__))
            return False
        if (not (self._min(key) <= value <= self._max(key))):
            warnings.warn('{0}: value out of range [{1}, {2}]'.format(
                    header, self._min(key), self._max(key)))
            return False
        return True

    def _validate_selection(self, key, value, verbose=True):
        # type: (str, str, bool) -> bool
        """Return True if 'value' is a valid element for the given
        parameter, False otherwise.

        Keyword Arguments:
            key     -- name of the parameter
            value   -- new value to be assigned to the parameter
            verbose -- print name of the parameter in warning messages
                       [default: True]
        """
        header = ' ({0})'.format(key) if verbose else ''
        if ((not isinstance(value, str))
            and (not isinstance(value, basestring))):
            warnings.warn("{0}: 'str' type expected, '{1}' found".format(
                    header, value.__class__.__name__))
            return False
        options = self._parameters[key]['options']
        if (value not in options):
            warnings.warn("{0}: '{1}' not in: {2}".format(
                    header, value, ', '.join(str(x) for x in options)))
            return False
        return True

    @staticmethod
    def _validate_path(key, value, verbose=True):
        # type: (str, str, bool) -> bool
        """Return True if 'value' is a valid path, False otherwise.

        Keyword Arguments:
            key     -- name of the parameter
            value   -- new value to be assigned to the parameter
            verbose -- print name of the parameter in warning messages
                       [default: True]
        """
        header = ' ({0})'.format(key) if verbose else ''
        if ((not isinstance(value, str))
            and (not isinstance(value, basestring))):
            warnings.warn("{0}: 'str' type expected, '{1}' found".format(
                    header, value.__class__.__name__))
            return False
        path = normalise_path(value)
        if (not os.path.isfile(path)):
            warnings.warn("{0}: No such file: '{1}'".format(header, path))
            return False
        return True

    def _validate_range(self, key, value, verbose=True):
        # type: (str, list, bool) -> bool
        """Return True if 'value' is a valid numeric range for the given
        parameter, False otherwise.

        Keyword Arguments:
            key     -- name of the parameter
            value   -- new value to be assigned to the parameter
            verbose -- print name of the parameter in warning messages
                       [default: True]
        """
        header = ' ({0})'.format(key) if verbose else ''
        if (not isinstance(value, list)):
            warnings.warn("{0}: 'list' type expected, '{1}' found".format(
                    header, value.__class__.__name__))
            return False
        if (len(value) != 2):
            warnings.warn("{0}: list of 2 elements expected, {1} found".format(
                    header, len(value)))
            return False
        # Get the elements' expected type from the first part of the
        # parameter's "type" field
        typeStr = self._parameters[key]['type'].split(' ')[0]
        if (any(not isinstance(x, self._eval_type(typeStr)) for x in value)):
            warnings.warn(("{0}: ['{1}','{1}'] types expected, ['{2}','{3}']"
                           " found".format(header, typeStr,
                                           value[0].__class__.__name__,
                                           value[1].__class__.__name__)))
            return False
        if (not (self._min(key) <= value[0] < value[1] <= self._max(key))):
            warnings.warn(("{0}: unfulfilled condition: {1} <= {2} < {3} <= "
                           "{4}".format(header, self._min(key), value[0],
                                       value[1], self._max(key))))
            return False
        return True

    def _validate_multiselection(self, key, value, verbose=True):
        # type: (str, list, bool) -> bool
        """Return True if 'value' is a list of valid elements for the
        given parameter, False otherwise.

        Keyword Arguments:
            key     -- name of the parameter
            value   -- new value to be assigned to the parameter
            verbose -- print name of the parameter in warning messages
                       [default: True]
        """
        header = ' ({0})'.format(key) if verbose else ''
        if (not isinstance(value, list)):
            warnings.warn("{0}: 'list' type expected, '{1}' found".format(
                    header, value.__class__.__name__))
            return False
        options = self._parameters[key]['options']
        if (any((x not in options for x in value))):
            warnings.warn(("{0}: at least one of the elements is not in: "
                           "{1}".format(header,
                                        ', '.join(str(x) for x in options))))
            return False
        return True

    def _validate_pairs(self, key, value, verbose=True):
        # type: (str, list, bool) -> bool
        """Return True if 'value' is a 2-tuples list of valid elements
        for the given parameter, False otherwise.

        Keyword Arguments:
            key     -- name of the parameter
            value   -- new value to be assigned to the parameter
            verbose -- print name of the parameter in warning messages
                       [default: True]
        """
        header = ' ({0})'.format(key) if verbose else ''
        if (not isinstance(value, list)):
            warnings.warn("{0}: 'list' type expected, '{1}' found".format(
                    header, value.__class__.__name__))
            return False
        # Check if each element of 'value' is a 2-tuple (list type)
        if (any((not isinstance(x, list)) or (len(x) != 2) for x in value)):
            warnings.warn('{0}: list of 2-tuples expected'.format(header))
            return False
        # Load the list of valid values for each 2-tuple from the first
        # column of the file whose path is saved as value of another
        # parameter ("file" key in the current parameter's information)
        srcFilePath = self._parameters[self._parameters[key]['file']]['value']
        options = pandas.read_csv(srcFilePath).iloc[:, 0].tolist()
        if (any(not set(x).issubset(options) for x in value)):
            warnings.warn(("{0}: at least one of the elements is not in: "
                           "{1}".format(header,
                                        ', '.join(str(x) for x in options))))
            return False
        return True

    def _is_active(self, key):
        # type: (str) -> bool
        """Return True if every parameter's "triggers" restriction is
        True, False otherwise.

        If the current value is None, the parameter is considered
        inactive (return False).

        Keyword Arguments:
            key -- name of the parameter
        """
        if ('triggers' in self._parameters[key]):
            try:
                return all(map(self._eval_expr,
                               self._parameters[key]['triggers']))
            except TypeError:
                return False
        else:
            return True

    def _min(self, key):
        # type: (str) -> object
        """Return the largest value in the parameter's "min" list.

        Keyword Arguments:
            key -- name of the parameter
        """
        if ('min' in self._parameters[key]):
            # Return the largest item among the minimum thresholds
            return round(max(map(
                    self._eval_expr, self._parameters[key]['min'])), 6)
        else:
            typeObj = eval(self._parameters[key]['type'].split(' ')[0])
            return typeObj(-1e6)

    def _max(self, key):
        # type: (str) -> object
        """Return the smallest value in the parameter's "max" list.

        Keyword Arguments:
            key -- name of the parameter
        """
        if ('max' in self._parameters[key]):
            # Return the smallest item among the maximum thresholds
            return round(min(map(
                    self._eval_expr, self._parameters[key]['max'])), 6)
        else:
            typeObj = eval(self._parameters[key]['type'].split(' ')[0])
            return typeObj(1e6)

    def _eval_expr(self, expr):
        # type: (object) -> object
        """Return the result of evaluating the expression.

        The expression can be either a literal (int, float, str) or a
        list. For the former, the result will be either the parameter
        value (if 'expr' is a key of parameters dictionary) or the
        literal itself. If 'expr' is a list, it must be in the following
        format:
            [parameter_name, operator_string, value]
        and the method will return the result of evaluating the
        operation applied to the parameter's value and the other value.

        Keyword Arguments:
            expr -- expression to evaluate
        """
        if (isinstance(expr, list)):
            strExpr = repr(self._parameters[expr[0]]['value']) + expr[1] \
                      + repr(expr[2])
            return eval(strExpr)
        else:
            if (expr in self._parameters):
                return self._parameters[expr]['value']
            else:
                return expr

    def _eval_type(self, typeStr):
        # type: (str) -> object
        """Returns the type corresponding to evaluate typeStr.

        If typeStr is "float", returns the tuple (int, float). Returns
        the corresponding type otherwise.

        Keyword Arguments:
            typeStr -- string with the type to return
        """
        if typeStr == 'float':
            return (int, float)
        else:
            return eval(typeStr)
