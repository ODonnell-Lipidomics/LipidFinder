# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Command-line Interface (CLI) to manage the parameters' collection."""

import ast
import os
import re
import warnings

import pandas

from LipidFinder.Configuration import LFParameters
from LipidFinder._py3k import input
from LipidFinder._utils import normalise_path, traceless_warning


class LFParametersCLI(LFParameters):
    """A LFParametersCLI object stores a set of LipidFinder parameters
    to be used in the specified module.

    This subclass of LFParameters implements a question-answering system
    (through a command-line interface) that is executed during the
    object creation. It allows the user to check (and change) each
    active parameter's value interactively.

    Attributes:
        _parameters  (Private[collections.OrderedDict])
            Dictionary where the parameters and their associated
            information are stored.
        _showHelp  (Private[bool])
            Adds help information to each parameter's description (if
            available).

    Examples:
        LFParametersCLI objects can be created as follows:
            >>> from Configuration.LFParametersCLI import
            ...     LFParametersCLI
            >>> LFParametersCLI()
            The next questions will help you review the parameters for
            LipidFinder's peakfilter module. A short description of the
            parameter and its current value are displayed. Next, you can
            either introduce the new value(s) or press ENTER, leaving
            the default one(s).


            File polarity mode:
            Options: Positive, Negative
            ...

        Alternatively, a specific module can be introduced as argument:
            >>> from Configuration.LFParametersCLI import
            ...     LFParametersCLI
            >>> LFParametersCLI(module='mssearch')
            The next questions will help you review the parameters for
            LipidFinder's mssearch module. A short description of the
            parameter and its current value are displayed. Next, you can
            either introduce the new value(s) or press ENTER, leaving
            the default one(s).


            Name (header) of the column that contains the m/z values:
            ...

        It also accepts a parameters JSON file, which will override the
        default values before starting the question-answering system:
            >>> from Configuration.LFParametersCLI import
            ...     LFParametersCLI
            >>> LFParametersCLI(src='/home/user/my_parameters.json')
            The next questions will help you review the parameters for
            LipidFinder's peakfilter module. A short description of the
            parameter and its current value are displayed. Next, you can
            either introduce the new value(s) or press ENTER, leaving
            the default one(s).


            File polarity mode:
            Options: Positive, Negative
            ...

        If the input is wrong, a warning message will be displayed and
        the input will be requested again:
            >>> LFParametersCLI()
            ...
            File polarity mode:
            Options: Positive, Negative
            Pos
            Warning: 'Pos' not in: Positive, Negative


            File polarity mode:
            Options: Positive, Negative
            ...

        At the end, the user will be asked for a file path were to save
        the new parameters (in JSON format):
            >>> LFParametersCLI()
            ...
            Where do you want to save the new set of parameters (path)?
            /home/user/new_parameters.json

            >>>

        The previous example would save the new parameters' values for
        PeakFilter to "/home/user/new_parameters.json".
    """

    def __init__(self, support=True, **kwargs):
        # type: (bool, ...) -> LFParametersCLI
        """Constructor of the class LFParametersCLI.

        First, the module's parameters template file is loaded. Next, if
        a source JSON parameters file path is provided, the default
        values are overwritten by the corresponding new (valid) values.
        Finally, the question-answering system is executed.

        Keyword Arguments:
            support -- display help information together with the
                       description [default: True]
        """
        LFParameters.__init__(self, **kwargs)
        self._showHelp = support
        # Launch a question-answering system to review the parameters
        print(("The next questions will help you review the parameters for "
               "LipidFinder's{0}{1} module. A short description of the "
               "parameter and its current value{0}are displayed. Next, you can "
               "either introduce the new value(s) or press ENTER,{0}leaving the"
               " default one(s).{0}").format(os.linesep, self._module))
        for key in (x for x in self._parameters.keys() if self._is_active(x)):
            data = self._parameters[key]
            typeStr = data['type']
            # Compose the main request message to display
            text = '{0}{1}{2}{0}'.format(os.linesep, data['description'],
                                         ' [y/n]' if typeStr == 'bool' else '')
            if (self._showHelp):
                # Show help if available for the parameter
                if (('help_cli' in data) and data['help_cli']):
                    # Some parameters might have an option available
                    # only for GUI, e.g. "Leave empty to include all."
                    text += 'Help: {0}{1}'.format(data['help_cli'], os.linesep)
                elif ('help' in data):
                    text += 'Help: {0}{1}'.format(data['help'], os.linesep)
            # Request the parameter value with a question-answering
            # system
            if (typeStr == 'bool'):
                # Add the current value to the request message
                text += '[current: {0}]  '.format('YES' if self[key] else 'NO')
                self._request_bool(key, text)
            elif (typeStr in ['int', 'float']):
                # Add the current value to the request message
                if (self[key] is not None):
                    # Add the current value to the request message
                    text += '[current: {0}]  '.format(self[key])
                self._request_number(key, text)
            elif (typeStr == 'selection'):
                # Add the available options to the request message
                text += 'Options: {0}{1}'.format(
                        self._from_list_to_str(data['options'], typeStr),
                        os.linesep)
                if (self[key]):
                    # Add the current value to the request message
                    text += '[current: {0}]  '.format(self[key])
                self._request_str_input(key, text)
            elif (typeStr in ['int range', 'float range']):
                if (self[key]):
                    # Add the current value to the request message
                    text += '[current: {0}]  '.format(
                            self._from_list_to_str(self[key], typeStr))
                self._request_range(key, text)
            elif (typeStr == 'multiselection'):
                # Add the available options to the request message
                text += 'Options: {0}{1}'.format(
                        self._from_list_to_str(data['options'], typeStr),
                        os.linesep)
                if (self[key]):
                 # Add the current value to the request message
                    text += '[current: {0}]  '.format(
                            self._from_list_to_str(self[key], typeStr))
                self._request_list_input(key, text)
            elif (typeStr == 'pairs'):
                # Load the list of valid values for each 2-tuple from
                # the first column of the file whose path is saved as
                # value of another parameter ("file" key in the current
                # parameter's information)
                srcFilePath = self._parameters[data['file']]['value']
                options = pandas.read_csv(srcFilePath).iloc[:, 0].tolist()
                # Add the available options to the request message
                text += 'Options: ' + self._from_list_to_str(options, typeStr) \
                        + os.linesep
                if (self[key]):
                    # Add the current value to the request message
                    text += '[current: {0}]  '.format(
                            self._from_list_to_str(self[key], typeStr))
                self._request_list_input(key, text)
            else:
                if (self[key]):
                    # Add the current value to the request message
                    text += '[current: {0}]  '.format(self[key])
                self._request_str_input(key, text)
        # Question-answering system to save current parameters' values
        text = ("{0}{0}Where do you want to save the new set of parameters "
                "(path)?{0}  ".format(os.linesep))
        while True:
            answer = input(text)
            if (answer):
                self.write(normalise_path(answer))
                break

    def _request_bool(self, key, text):
        # type: (str, str) -> None
        """Request a boolean value as "yes" or "no" answer.

        Keyword Arguments:
            key  -- name of the parameter
            text -- request message to display
        """
        while True:
            answer = input(text)
            if (not answer):
                # Leave default value
                break
            else:
                if (answer.lower() in ['yes', 'y', 'no', 'n']):
                    answer = answer.lower()
                    # Translate "yes"/"y" to True and "no"/"n" to False
                    self._parameters[key]['value'] = answer.startswith('y')
                    break
                else:
                    warnings.warn((": 'yes', 'y', 'no' or 'n' expected, '{0}' "
                                   "inputted").format(answer))

    def _request_number(self, key, text):
        # type: (str, str) -> None
        """Request a number value, either int or float.

        Keyword Arguments:
            key  -- name of the parameter
            text -- request message to display
        """
        while True:
            answer = input(text)
            if (not answer):
                if (self[key] is None):
                    warnings.warn(': there is no default value')
                else:
                    # Leave default value
                    break
            else:
                try:
                    answer = ast.literal_eval(answer)
                except ValueError:
                    typeStr = self._parameters[key]['type']
                    warnings.warn(": '{0}' type expected".format(typeStr))
                else:
                    # The validation process will warn about what was
                    # wrong if 'answer' is not valid
                    if (self._validate_number(key, answer, False)):
                        self._parameters[key]['value'] = answer
                        break

    def _request_str_input(self, key, text):
        # type: (str, str) -> None
        """Request a string, either a path, one among several options or
        a plain string.

        Keyword Arguments:
            key  -- name of the parameter
            text -- request message to display
        """
        while True:
            answer = input(text)
            if (not answer):
                if (not self[key]):
                    warnings.warn(': there is no default value')
                else:
                    # Leave default value
                    break
            else:
                typeStr = self._parameters[key]['type']
                # The validation process will warn about what was wrong
                # if 'answer' is not valid
                if ((typeStr == 'str')
                    and self._validate_literal(key, answer, False)):
                    self._parameters[key]['value'] = answer
                    break
                elif ((typeStr == 'selection')
                      and self._validate_selection(key, answer, False)):
                    self._parameters[key]['value'] = answer
                    break
                elif ((typeStr == 'path')
                      and self._validate_path(key, answer, False)):
                        self._parameters[key]['value'] = normalise_path(answer)
                        break

    def _request_range(self, key, text):
        # type: (str, str) -> None
        """Request a numeric range of two elements (either integers or
        floats).

        Both numbers must be of the same type.

        Keyword Arguments:
            key  -- name of the parameter
            text -- request message to display
        """
        while True:
            answer = input(text)
            if (not answer):
                if (not self[key]):
                    warnings.warn(': there is no default value')
                else:
                    # Leave default value
                    break
            else:
                answer = self._from_str_to_list(answer)
                try:
                    # Convert each string in 'answer' to its
                    # corresponding numeric type (int or float)
                    answer = [ast.literal_eval(x) for x in answer]
                except SyntaxError:
                    warnings.warn(": each number must be separated by ','")
                except ValueError:
                    typeStr = self._parameters[key]['type'].split(' ')[0]
                    warnings.warn(": list of two '{0}' types expected".format(
                            typeStr))
                else:
                    # The validation process will warn about what was
                    # wrong if 'answer' is not valid
                    if (self._validate_range(key, answer, False)):
                        self._parameters[key]['value'] = answer
                        break

    def _request_list_input(self, key, text):
        # type: (str, str) -> None
        """Request a list of values, either strings or 2-tuples of
        strings.

        Keyword Arguments:
            key  -- name of the parameter
            text -- request message to display
        """
        while True:
            answer = input(text)
            if (not answer):
                if (not self[key]):
                    warnings.warn(': there is no default value')
                else:
                    # Leave default value
                    break
            else:
                typeStr = self._parameters[key]['type']
                answer = self._from_str_to_list(answer, typeStr)
                # The validation process will warn about what was wrong
                # if 'answer' is not valid
                if ((typeStr == 'pairs')
                    and self._validate_pairs(key, answer, False)):
                    self._parameters[key]['value'] = answer
                    break
                elif ((typeStr == 'multiselection')
                      and self._validate_multiselection(key, answer, False)):
                    self._parameters[key]['value'] = answer
                    break

    @staticmethod
    def _from_str_to_list(value, typeStr):
        # type: (str, str) -> list
        """Return a Python's list from a user-friendly list in string
        format.

        Keyword Arguments:
            value   -- user-friendly string list
            typeStr -- expected parameter type
        """
        if (typeStr == 'pairs'):
            procValue = value.strip()
            procValue = re.sub(r"\([ \t]*", '["', procValue)
            procValue = re.sub(r"[ \t]*\)", '"]', procValue)
            procValue = re.sub(r"(?P<a>[^\] \t])[ \t]*,[ \t]*(?P<b>[^\[ \t])",
                               '\g<a>","\g<b>', procValue)
            procValue = '[{0}]'.format(procValue)
        elif (typeStr == 'multiselection'):
            procValue = value.strip()
            procValue = re.sub(r"[ \t]*,[ \t]*", '","', procValue)
            procValue = '["{0}"]'.format(procValue)
        try:
            return list(ast.literal_eval(procValue))
        except SyntaxError:
            # The user-friendly list in string format is not formed
            # correctly
            return [value]

    @staticmethod
    def _from_list_to_str(value, typeStr):
        # type: (list, str) -> str
        """Return a user-friendly list in string format from a Python's
        list.

        Keyword Arguments:
            value   -- Python's list
            typeStr -- expected parameter type
        """
        valueStr = repr(value)
        valueStr = valueStr.replace("u'", "")
        valueStr = valueStr.replace("'", "")
        if (typeStr != 'multiselection'):
            valueStr = valueStr.replace("[", "(")
            valueStr = valueStr.replace("]", ")")
        return valueStr[1:-1]
