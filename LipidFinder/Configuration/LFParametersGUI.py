# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.
"""Graphical User Interface (GUI) to manage the parameters' collection.
"""

from collections import OrderedDict
import os

from IPython.display import display
from ipywidgets import widgets, Layout
import pandas

from LipidFinder.Configuration import LFParameters
from LipidFinder._utils import normalise_path


class _TaggedToggleButton(widgets.ToggleButton):
    """Add "tag" attribute to widgets.ToggleButton class."""

    def __init__(self, tag, **kwargs):
        widgets.ToggleButton.__init__(self, **kwargs)
        self.tag = tag


class _TaggedCheckbox(widgets.Checkbox):
    """Add "tag" attribute to widgets.Checkbox class."""

    def __init__(self, tag, **kwargs):
        widgets.Checkbox.__init__(self, **kwargs)
        self.tag = tag


class _TaggedButton(widgets.Button):
    """Add "tag" attribute to widgets.Button class."""

    def __init__(self, tag, **kwargs):
        widgets.Button.__init__(self, **kwargs)
        self.tag = tag


class LFParametersGUI(LFParameters):
    """A LFParametersGUI object stores a set of LipidFinder parameters
    to be used in the specified module.

    This subclass of LFParameters implements a graphical interface using
    jupyter notebook's widgets, executed during the object creation. It
    allows the user to check, change and save each active parameter's
    value interactively.

    Attributes:
        _parameters  (Private[collections.OrderedDict])
            Dictionary where the parameters and their associated
            information are stored.
        _floatPointPrecision  (Private[int])
            Number of digits after the radix point in floats.
        _floatStep  (Private[float])
            Minimum difference between two consecutive float numbers.
        _style  (Private[dict])
            Dictionary with the default style settings for widgets.
        _inputWidth  (Private[str])
            String representation of the default width of input widgets.
        _widgets  (Private[collections.OrderedDict])
            Dictionary where the widgets for each parameter are stored.

    Examples:
        LFParametersGUI objects can be created as follows:
            >>> from Configuration.LFParametersGUI import
            ...     LFParametersGUI
            >>> LFParametersGUI()
            >>> LFParametersGUI(src='/home/user/my_parameters.json')

        The former will load the default PeakFilter parameters and will
        load and display the interface afterwards. The latter will load
        the default PeakFilter parameters, override them with the values
        found in the JSON file provided, and finally it will load and
        display the interface.

        Alternatively, a specific module can be introduced as argument:
            >>> from Configuration.LFParametersGUI import
            ...     LFParametersGUI
            >>> LFParametersGUI(module='mssearch')
    """

    def __init__(self, precision=4, **kwargs):
        # type: (int, ...) -> LFParametersGUI
        """Constructor of the class LFParametersGUI.

        First, the module's parameters template file is loaded. Next, if
        a source JSON parameters file path is provided, the default
        values are overwritten by the corresponding new (valid) values.
        Finally, the graphical user interface is displayed.

        Keyword Arguments:
            precision -- number of decimal digits to use with floats
                         (e.g. a precision of 2 forces a difference of
                          0.01 between any two consecutive float numbers)
                         [default: 4]
        """
        # Minimum difference between two consecutive float numbers
        self._floatPointPrecision = precision
        self._floatStep = 10 ** -(precision)
        # Load the parameters dictionary using parent class' constructor
        LFParameters.__init__(self, **kwargs)
        # Default style
        self._style = {'description_width': '0px'}
        # Default width of input widgets
        self._inputWidth = '26%'
        # Generate an ordered dict to store each parameter's set of
        # widgets in the same order as in the parameters' dict
        self._widgets = OrderedDict()
        # Create every widget of the GUI
        for key, data in self._parameters.items():
            disabled = not self._is_active(key)
            # Load the information of each parameter
            self._widgets[key] = [self._create_label(key, disabled),
                                  self._create_help_icon(key, disabled)]
            # Create the input widget or container of input widgets for
            # each parameter type
            if (data['type'] == 'bool'):
                self._widgets[key].append(
                        self._create_bool_widget(key, disabled))
            elif (data['type'] == 'int'):
                self._widgets[key].append(
                        self._create_int_widget(key, disabled))
            elif (data['type'] == 'float'):
                self._widgets[key].append(
                        self._create_float_widget(key, disabled))
            elif (data['type'] == 'selection'):
                self._widgets[key].append(
                        self._create_selection_widget(key, disabled))
            elif (data['type'] == 'path'):
                self._widgets[key].append(
                        self._create_path_widget(key, disabled))
            elif (data['type'] == 'int range'):
                self._widgets[key].append(
                        self._create_int_range_widget(key, disabled))
            elif (data['type'] == 'float range'):
                self._widgets[key].append(
                        self._create_float_range_widget(key, disabled))
            elif (data['type'] == 'multiselection'):
                self._widgets[key].append(
                        self._create_multiselection_widget(key, disabled))
            elif (data['type'] == 'pairs'):
                self._widgets[key].append(
                        self._create_pairs_widget(key, disabled))
            else: # data['type'] == 'str'
                self._widgets[key].append(
                        self._create_str_widget(key, disabled))
        # Display the GUI
        hboxLayout = Layout(align_items='center')
        for key, widgetList in self._widgets.items():
            display(widgets.HBox(widgetList, layout=hboxLayout))
        # Finally, create the save interface to allow the user to save
        # the current parameters values in a JSON file
        display(widgets.HBox([], layout=Layout(height='15px')))
        display(widgets.HBox([], layout=Layout(height='0px',
                                               border='2px solid lightgray')))
        display(widgets.HBox([], layout=Layout(height='2px')))
        self._widgets['save'] = self._create_save_widget()
        hboxLayout = Layout(justify_content='space-between',
                            align_items='center')
        display(widgets.HBox(self._widgets['save'], layout=hboxLayout))

    def _create_label(self, key, disabled):
        # type: (str, bool) -> widgets.HTML
        """Return an HTML widget with the parameter's description.

        If 'disabled' is False, the text will be in black, otherwise it
        will be in gray.

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        text = self._parameters[key]['description']
        label = ("<p style=\"font-size:110%; line-height:19px; color:{0};\">{1}"
                 "</p>").format('Gray' if disabled else 'Black', text)
        return widgets.HTML(value=label, style=self._style,
                            layout=Layout(width='50%'))

    def _create_help_icon(self, key, disabled):
        # type: (str, bool) -> widgets.HTML
        """Return an HTML widget with the parameter's help as tooltip of
        a help icon.

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        if ('help' in self._parameters[key]):
            code = ("<link rel=\"stylesheet\" href=\"https://fonts.googleapis.c"
                    "om/icon?family=Material+Icons\"><i class=\"material-icons"
                    "\" style=\"color:{0}; font-size:18px; display:inline"
                    "-flex; vertical-align:middle;\" title=\"{1}\">help</i>"
                    "").format("SteelBlue", self._parameters[key]['help'])
        else:
            code = ''
        layout = Layout(width='2%',
                        visibility='hidden' if disabled else 'visible')
        return widgets.HTML(value=code, style=self._style, layout=layout)

    def _create_str_widget(self, key, disabled):
        # type: (str, bool) -> widgets.Text
        """Return a Text widget with the parameter's value.

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        if ('example' in self._parameters[key]):
            example = self._parameters[key]['example']
        else:
            example = ''
        inputWidget = widgets.Text(
                value=self[key], description=key, placeholder=example,
                style=self._style, layout=Layout(width=self._inputWidth),
                continuous_update=False, disabled=disabled)
        # Add handler for when the "value" trait changes
        inputWidget.observe(self._default_handler, names='value')
        return inputWidget

    def _create_bool_widget(self, key, disabled):
        # type: (str, bool) -> widgets.HBox
        """Return an HBox containing a ToggleButton widget to represent
        the parameter's value.

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        inputWidget = _TaggedToggleButton(
                value=self[key], description='Yes' if self[key] else 'No',
                tag=key, style=self._style, layout=Layout(width='50%'),
                button_style='primary', disabled=disabled)
        # Add handler for when the "value" trait changes
        inputWidget.observe(self._bool_handler, names='value')
        layout = Layout(width=self._inputWidth, justify_content='center')
        return widgets.HBox([inputWidget], layout=layout)

    def _create_int_widget(self, key, disabled):
        # type: (str, bool) -> widgets.BoundedIntText
        """Return a BoundedIntText widget with the parameter's value.

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        inputWidget = widgets.BoundedIntText(
                value=self[key], description=key, min=self._min(key),
                max=self._max(key), style=self._style,
                layout=Layout(width=self._inputWidth), continuous_update=False,
                disabled=disabled)
        # Save the widget's value in case its constructor automatically
        # replaces an empty one given as argument
        self._parameters[key]['value'] = inputWidget.value
        # Add handler for when the "value" trait changes
        inputWidget.observe(self._default_handler, names='value')
        return inputWidget

    def _create_float_widget(self, key, disabled):
        # type: (str, bool) -> widgets.BoundedFloatText
        """Return a BoundedFloatText widget with the parameter's value.

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        inputWidget = widgets.BoundedFloatText(
                value=self[key], description=key, min=self._min(key),
                max=self._max(key), step=self._floatStep, style=self._style,
                layout=Layout(width=self._inputWidth), continuous_update=False,
                disabled=disabled)
        # Save the widget's value in case its constructor automatically
        # replaces an empty one given as argument
        self._parameters[key]['value'] = inputWidget.value
        # Add handler for when the "value" trait changes
        inputWidget.observe(self._default_handler, names='value')
        return inputWidget

    def _create_selection_widget(self, key, disabled):
        # type: (str, bool) -> widgets.Dropdown
        """Return a Dropdown widget with the parameter's options and its
        current value selected.

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        inputWidget = widgets.Dropdown(
                options=self._parameters[key]['options'], value=self[key],
                description=key, style=self._style,
                layout=Layout(width=self._inputWidth), disabled=disabled)
        # Add handler for when the "value" trait changes
        inputWidget.observe(self._default_handler, names='value')
        return inputWidget

    def _create_path_widget(self, key, disabled):
        # type: (str, bool) -> widgets.HBox
        """Return an HBox containing a Text widget with the parameter's
        value.

        If the Text widget is enabled and the file does not exist, a
        warning icon will be displayed next to it to alert the user.

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        inputWidget = widgets.Text(
                value=self[key], description=key, style=self._style,
                layout=Layout(width='92%'), continuous_update=False,
                disabled=disabled)
        # Add handler for when the "value" trait changes
        inputWidget.observe(self._path_handler, names='value')
        # Create an HTML widget with a warning icon that will be
        # displayed if the Text widget is enabled and the file does not
        # exist
        code = ("<link rel=\"stylesheet\" href=\"https://fonts.googleapis.com/i"
                "con?family=Material+Icons\"><i class=\"material-icons\" style="
                "\"font-size:18px; color:Red; display:inline-flex; vertical-ali"
                "gn:middle;\" title=\"File not found!\">warning</i>")
        warn = not disabled and not os.path.isfile(self[key])
        layout = Layout(width='5%',
                        visibility='visible' if warn else 'hidden')
        warnWidget = widgets.HTML(value=code, style=self._style, layout=layout)
        layout = Layout(width='46%', justify_content='space-between')
        return widgets.HBox([inputWidget, warnWidget], layout=layout)

    def _create_int_range_widget(self, key, disabled):
        # type: (str, bool) -> widgets.HBox
        """Return an HBox containing two BoundedIntText widgets with the
        parameter's range values.

        The widgets are created to fulfill the "int range" type
        condition: lower_bound < upper_bound

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        lowerBound = widgets.BoundedIntText(
                value=self[key][0], description=key, min=self._min(key),
                max=self[key][1] - 1, style=self._style,
                layout=Layout(width='50%'), continuous_update=False,
                disabled=disabled)
        # Save the widget's value in case its constructor automatically
        # replaces an empty one given as argument
        self._parameters[key]['value'][0] = lowerBound.value
        # Add handler for when the "value" trait changes
        lowerBound.observe(self._range_handler, names='value')
        upperBound = widgets.BoundedIntText(
                value=self[key][1], description=key, min=self[key][0] + 1,
                max=self._max(key), style=self._style,
                layout=Layout(width='50%'), continuous_update=False,
                disabled=disabled)
        # Save the widget's value in case its constructor automatically
        # replaces an empty one given as argument
        self._parameters[key]['value'][1] = upperBound.value
        # Add handler for when the "value" trait changes
        upperBound.observe(self._range_handler, names='value')
        return widgets.HBox([lowerBound, upperBound],
                            layout=Layout(width=self._inputWidth))

    def _create_float_range_widget(self, key, disabled):
        # type: (str, bool) -> widgets.HBox
        """Return an HBox containing two BoundedFloatText widgets with
        the parameter's range values.

        The widgets are created to fulfill the "float range" type
        condition: lower_bound < upper_bound

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        lowerBound = widgets.BoundedFloatText(
                value=self[key][0], description=key, min=self._min(key),
                max=self[key][1] - self._floatStep, step=self._floatStep,
                style=self._style, layout=Layout(width='50%'),
                continuous_update=False, disabled=disabled)
        # Save the widget's value in case its constructor automatically
        # replaces an empty one given as argument
        self._parameters[key]['value'][0] = lowerBound.value
        # Add handler for when the "value" trait changes
        lowerBound.observe(self._range_handler, names='value')
        upperBound = widgets.BoundedFloatText(
                value=self[key][1], description=key,
                min=self[key][0] + self._floatStep, max=self._max(key),
                step=self._floatStep, style=self._style,
                layout=Layout(width='50%'), continuous_update=False,
                disabled=disabled)
        # Save the widget's value in case its constructor automatically
        # replaces an empty one given as argument
        self._parameters[key]['value'][1] = upperBound.value
        # Add handler for when the "value" trait changes
        upperBound.observe(self._range_handler, names='value')
        return widgets.HBox([lowerBound, upperBound],
                            layout=Layout(width=self._inputWidth))

    def _create_multiselection_widget(self, key, disabled):
        # type: (str, bool) -> widgets.Box
        """Return a Box containing as many Checkbox widgets as
        parameter's options, with those in its "value" field checked.

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        itemWidgets = []
        for item in self._parameters[key]['options']:
            layoutWidth = '23%' if (len(item) <= 10) else '48%'
            inputWidget = _TaggedCheckbox(
                    value=item in self[key], description=item, tag=key,
                    style=self._style, layout=Layout(width=layoutWidth),
                    disabled=disabled)
            # Add handler for when the "value" trait changes
            inputWidget.observe(self._multiselection_handler, names='value')
            itemWidgets.append(inputWidget)
        layout = Layout(width='46%', display='flex', flex_flow='row wrap',
                        justify_content='space-between')
        return widgets.Box(itemWidgets, layout=layout)

    def _create_pairs_widget(self, key, disabled):
        # type: (str, bool) -> widgets.HBox
        """Return an HBox containing the interface to add and remove
        pairs of available elements.

        The term "available elements" refers to those elements in the
        first column of the CSV file's path stored under the parameter's
        "file" key. Users will not be able to add existing pairs or
        pairs formed by the same element twice.

        Keyword Arguments:
            key      -- name of the parameter
            disabled -- is the parameter/widget disabled?
        """
        # Load the list of available elements from the first column of
        # the CSV file saved under the parameter's "file" key
        srcFilePath = self[self._parameters[key]['file']]
        options = pandas.read_csv(srcFilePath).iloc[:, 0].tolist()
        # Create two Select widgets with the list of available elements
        leftSelect = widgets.Select(
                options=options, rows=4, style=self._style,
                layout=Layout(width='20%'), disabled=disabled)
        rightSelect = widgets.Select(
                options=options, rows=4, style=self._style,
                layout=Layout(width='20%'), disabled=disabled)
        # Create the add and remove buttons with the handler to add and
        # remove pairs, respectively
        addButton = _TaggedButton(
                description='Pair >>', tooltip='Add new pair', tag=key,
                layout=Layout(width='95%'), disabled=disabled)
        # Add handlerfor when the button is clicked
        addButton.on_click(self._pairs_add_handler)
        delButton = _TaggedButton(
                description='<< Remove', tooltip='Remove selected pair',
                tag=key, layout=Layout(width='95%'), disabled=disabled)
        # Add handler for when the button is clicked
        delButton.on_click(self._pairs_del_handler)
        layout = Layout(width='21%', justify_content='space-around')
        # Hold the buttons in a VBox to get the desired layout
        buttonsBox = widgets.VBox([addButton, delButton], layout=layout)
        # Create a Select widget with the parameter's list of pairs
        pairs = [' , '.join(x) for x in self[key]]
        pairsSelect = widgets.Select(
                options=pairs, rows=4, style=self._style,
                layout=Layout(width='28%'), disabled=disabled)
        layout = Layout(width='46%', justify_content='space-around')
        return widgets.HBox([leftSelect, rightSelect, buttonsBox, pairsSelect],
                            layout=layout)

    def _create_save_widget(self):
        # type: () -> list
        """Return a list containing the interface to save the current
        parameters values as a JSON file in an introduced path.
        """
        text = ("<p style=\"font-size:110%; line-height:19px; color:Black;\">"
                 "Where do you want to save the new set of parameters?</p>")
        label = widgets.HTML(value=text, style=self._style,
                             layout=Layout(width='38%'))
        # Create the path input widget (Text) with a default path and
        # file name
        defaultPath = normalise_path("parameters.json")
        inputWidget = widgets.Text(
                value=defaultPath, placeholder=defaultPath, style=self._style,
                layout=Layout(width='40%'), continuous_update=False)
        # Add handler for when the "value" trait changes
        inputWidget.observe(self._save_path_handler, names='value')
        # Create an HTML widget with a warning icon that will be
        # displayed if the directory path does not exist
        code = ("<link rel=\"stylesheet\" href=\"https://fonts.googleapis.com/i"
                "con?family=Material+Icons\"><i class=\"material-icons\" style="
                "\"font-size:18px; color:Red; display:inline-flex; vertical-ali"
                "gn:middle;\" title=\"Path not found!\">warning</i>")
        dirPath = os.path.split(inputWidget.value)[0]
        visibility = 'visible' if not os.path.isdir(dirPath) else 'hidden'
        layout = Layout(width='2%', visibility=visibility)
        warnWidget = widgets.HTML(value=code, style=self._style, layout=layout)
        # Create a save button that will be active only if every active
        # parameter is valid and the destination path exists
        saveButton = widgets.Button(
                description='Save', button_style='danger',
                tooltip='Save parameters in a JSON file',
                layout=Layout(width='12%', height='35px'),
                disabled=not self._valid_parameters())
        # Add handler for when the button is clicked
        saveButton.on_click(self._save_button_handler)
        return [label, inputWidget, warnWidget, saveButton]

    def _update(self):
        # type: () -> None
        """Return an HBox containing the interface to add and remove
        pairs of available elements.

        The term "available elements" refers to those elements in the
        first column of the CSV file's path stored under the parameter's
        "file" key. Users will not be able to add existing pairs or
        pairs formed by the same element twice. If the CSV file path
        changes, the pairs list will be emptied and the set of available
        elements will be updated.
        """
        # Update the status and/or visibility of each parameter's widget
        for key in self._parameters.keys():
            interface = self._widgets[key]
            disabled = not self._is_active(key)
            if (disabled):
                interface[0].value = interface[0].value.replace('Black', 'Gray')
            else:
                interface[0].value = interface[0].value.replace('Gray', 'Black')
            interface[1].layout.visibility = 'hidden' if disabled else 'visible'
            typeStr = self._parameters[key]['type']
            if (typeStr == 'bool'):
                interface[2].children[0].disabled = disabled
            elif (typeStr in ['int', 'float']):
                # Update minimum and maximum bounds too
                interface[2].min = self._min(key)
                interface[2].max = self._max(key)
                interface[2].disabled = disabled
            elif (typeStr == 'path'):
                interface[2].children[0].disabled = disabled
                # Display the warning widget if the parameter is enabled
                # and the file does not exist
                if (not disabled and not os.path.isfile(self[key])):
                    interface[2].children[1].layout.visibility = 'visible'
                else:
                    interface[2].children[1].layout.visibility = 'hidden'
            elif (typeStr in ['int range', 'float range']):
                # Update minimum and maximum bounds of the range too
                interface[2].children[0].min = self._min(key)
                interface[2].children[0].disabled = disabled
                interface[2].children[1].max = self._max(key)
                interface[2].children[1].disabled = disabled
            elif (typeStr == 'multiselection'):
                for child in interface[2].children:
                    child.disabled = disabled
            elif (typeStr == 'pairs'):
                interface[2].children[0].disabled = disabled
                interface[2].children[1].disabled = disabled
                for grandchild in interface[2].children[2].children:
                    grandchild.disabled = disabled
                interface[2].children[3].disabled = disabled
            else:
                interface[2].disabled = disabled
        # Ensure the save button should be available and ready to save
        # the new set of parameters
        self._widgets['save'][3].description = 'Save'
        self._widgets['save'][3].icon = ''
        self._widgets['save'][3].disabled = not self._valid_parameters()

    def _default_handler(self, change):
        # type: (dict) -> None
        """Handle the "value" trait change assigning the new value to
        the corresponding parameter.

        The update() method is launched at the end to ensure every
        widget is updated according to the change in this parameter.

        Keyword Arguments:
            change -- dict holding the information about the change
        """
        key = change['owner'].description
        self._parameters[key]['value'] = change['new']
        self._update()

    def _bool_handler(self, change):
        # type: (dict) -> None
        """Handle the "value" trait change assigning the new value to
        the corresponding "bool" type parameter.

        The update() method is launched at the end to ensure every
        widget is updated according to the change in this parameter.

        Keyword Arguments:
            change -- dict holding the information about the change
        """
        key = change['owner'].tag
        self._parameters[key]['value'] = change['new']
        # Change ToggleButton's description to "Yes" or "No" depending
        # on whether its new value is True or False, respectively
        change['owner'].description = 'Yes' if change['new'] else 'No'
        self._update()

    def _path_handler(self, change):
        # type: (dict) -> None
        """Handle the "value" trait change assigning the new value to
        the corresponding "path" type parameter.

        The update() method is launched at the end to ensure every
        widget is updated according to the change in this parameter.

        Keyword Arguments:
            change -- dict holding the information about the change
        """
        key = change['owner'].description
        self._parameters[key]['value'] = normalise_path(change['new'])
        # Replace the introduced path by its normalised version to
        # provide the user with more information in case there is
        # something wrong with the path
        change['owner'].value = self[key]
        # Get the "pairs" type parameter that has this parameter in its
        # "field" key to update the contents of its widgets
        for param, data in self._parameters.items():
            if ((data['type'] == 'pairs') and (data['file'] == key)):
                pairsWidget = self._widgets[param][2]
                if (os.path.isfile(self[key])):
                    # Update the information of available elements
                    options = pandas.read_csv(self[key]).iloc[:, 0].tolist()
                    pairsWidget.children[0].options = options
                    pairsWidget.children[1].options = options
                else:
                    # Since the file does not exist, there are no
                    # available elements
                    pairsWidget.children[0].options = []
                    pairsWidget.children[1].options = []
                # Since the file has changed, empty the list of pairs
                self._parameters[param]['value'] = []
                pairsWidget.children[3].options = []
                break
        self._update()

    def _range_handler(self, change):
        # type: (dict) -> None
        """Handle the "value" trait change assigning the new value to
        the corresponding "int/float range" type parameter.

        The update() method is launched at the end to ensure every
        widget is updated according to the change in this parameter.

        Keyword Arguments:
            change -- dict holding the information about the change
        """
        key = change['owner'].description
        # Both children have the same step
        step = self._widgets[key][2].children[0].step
        if (change['owner'].min == self._min(key)):
            # Trait changed in the widget corresponding to the lower
            # bound of the range
            self._parameters[key]['value'][0] = change['new']
            self._widgets[key][2].children[1].min = change['new'] + step
        else:
            # Trait changed in the widget corresponding to the upper
            # bound of the range
            self._parameters[key]['value'][1] = change['new']
            self._widgets[key][2].children[0].max = change['new'] - step
        self._update()

    def _multiselection_handler(self, change):
        # type: (dict) -> None
        """Handle the "value" trait change updating the list of values
        of the corresponding "multiselection" type parameter.

        The update() method is launched at the end to ensure every
        widget is updated according to the change in this parameter.

        Keyword Arguments:
            change -- dict holding the information about the change
        """
        key = change['owner'].tag
        if (change['new']):
            self._parameters[key]['value'].append(change['owner'].description)
        else:
            self._parameters[key]['value'].remove(change['owner'].description)
        self._update()

    def _pairs_add_handler(self, button):
        # type: (_TaggedButton) -> None
        """Handle when the button is clicked to add a pair to the
        corresponding "pairs" type parameter.

        The update() method is launched at the end to ensure every
        widget is updated according to the change in this parameter.

        Keyword Arguments:
            button -- clicked button widget instance
        """
        key = button.tag
        # Add selected elements in both Selection widgets as a new pair
        leftSel = self._widgets[key][2].children[0].value
        rightSel = self._widgets[key][2].children[1].value
        newPair = [leftSel, rightSel]
        # The pairs are considered sets, that is, the order of the
        # elements is ignored
        if ((leftSel != rightSel) and (newPair not in self[key])
            and (newPair[::-1] not in self[key])):
            self._parameters[key]['value'].append(newPair)
            # Since the "options" field is a tuple, build a new list
            # with the new pair
            self._widgets[key][2].children[3].options = \
                    [' , '.join(x) for x in self[key]]
        self._update()

    def _pairs_del_handler(self, button):
        # type: (_TaggedButton) -> None
        """Handle when the button is clicked to remove a pair of the
        corresponding "pairs" type parameter.

        The update() method is launched at the end to ensure every
        widget is updated according to the change in this parameter.

        Keyword Arguments:
            button -- clicked button widget instance
        """
        key = button.tag
        pairsWidget = self._widgets[key][2].children[3]
        # Get the selected pair from the pairs widget
        pairSel = pairsWidget.value
        if (pairSel):
            pair = pairSel.split(' , ')
            self._parameters[key]['value'].remove(pair)
            # Since the "options" field is a tuple, build a new list
            # without the deleted pair
            pairsWidget.options = [' , '.join(x) for x in self[key]]
            # Select the first pair to ensure coherence with the change
            if (pairsWidget.options):
                pairsWidget.value = pairsWidget.options[0]
        self._update()

    def _save_path_handler(self, change):
        # type: (dict) -> None
        """Handle the "value" trait change checking if the path where to
        save the parameters values exists.

        A warning sign will be displayed if the given directory path
        does not exist. The update() method is launched at the end to
        ensure every widget is updated according to the change in this
        parameter.

        Keyword Arguments:
            change -- dict holding the information about the change
        """
        newPath = normalise_path(change['new'])
        dirPath = os.path.split(newPath)[0]
        if (not os.path.isdir(dirPath)):
            self._widgets['save'][2].layout.visibility = 'visible'
        else:
            self._widgets['save'][2].layout.visibility = 'hidden'
        # Replace the introduced path by its normalised version to
        # provide the user with more information in case there is
        # something wrong
        change['owner'].value = newPath
        self._update()

    def _save_button_handler(self, button):
        # type: (widgets.Button) -> None
        """Handle when the button is clicked to save the parameters
        values in a JSON file.

        Keyword Arguments:
            button -- clicked button widget instance
        """
        self.write(self._widgets['save'][1].value)
        # Change the button's text to tell the user the JSON parameters
        # file has been correctly created
        button.description = 'Saved'
        button.icon = 'check'

    def _min(self, key):
        # type: (str) -> object
        """Return the largest value in the parameter's "min" list.

        Applies round() method to the output of LFParameter._max() to
        get a more comparable result regarding floating point arithmetic
        issues.

        Keyword Arguments:
            key -- name of the parameter
        """
        return round(LFParameters._min(self, key), self._floatPointPrecision)

    def _max(self, key):
        # type: (str) -> object
        """Return the smallest value in the parameter's "max" list.

        Applies round() method to the output of LFParameter._max() to
        get a more comparable result regarding floating point arithmetic
        issues.

        Keyword Arguments:
            key -- name of the parameter
        """
        return round(LFParameters._max(self, key), self._floatPointPrecision)

    def _valid_parameters(self):
        # type: () -> bool
        """Return True if every active parameter has a valid value,
        False otherwise.

        The list of valid parameters also includes "save" destination
        path, where the JSON parameters file will be saved.
        """
        enabledKeys = (x for x in self._parameters.keys() if self._is_active(x))
        for key in enabledKeys:
            data = self._parameters[key]
            # Only "multiselection" type parameters can be empty ([])
            if ((data['type'] != 'multiselection')
                and (data['value'] in [None, '', []])):
                return False
            # "path" type parameters must be checked manually, whilst
            # the rest are already controlled by their widget
            if ((data['type'] == 'path') and not os.path.isfile(data['value'])):
                return False
        # This method is also called when the save interface is being
        # created, so the "save" key will not exist yet
        if ('save' in self._widgets):
            # Check if the directory path where to save the JSON
            # parameters file exists
            dirPath = os.path.split(self._widgets['save'][1].value)[0]
            if (not os.path.isdir(dirPath)):
                return False
        return True
