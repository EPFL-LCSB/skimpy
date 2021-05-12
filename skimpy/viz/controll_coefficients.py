# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2021 Laboratory of Computational Systems Biotechnology (LCSB),
Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import numpy as np

from bokeh.plotting import figure, output_file
from bokeh.models import ColumnDataSource
from bokeh.io import curdoc, show


def plot_control_coefficients(control_coeffs,
                              filename='map.html',
                              top=10,
                              color='#3399FF',
                              backend='webgl',
                              x_label='',
                              y_label='',
                              **kwargs):
    """

    :param control_coeffs: pandas Series
    :param top:
    :param sorted:
    :param kwargs:
    :return:
    """
    if not top is None:
        temp_data = control_coeffs.abs().sort_values()
        temp_data = temp_data.iloc[-top:]
        data = control_coeffs[temp_data.index]
    else:
        data = control_coeffs

    # output to static HTML file
    output_file(filename)

    source = ColumnDataSource(dict(y=list(data.index),
                                   right=data.values,))

    plot = figure(
        y_range=list(data.index),
        **kwargs,)

    plot.hbar(y="y",
              right="right",
              left=0,
              height=0.5,
              fill_color=color,
              source=source)

    # Make the axsis pretty
    plot.xaxis.axis_label = x_label
    plot.yaxis.axis_label = y_label

    plot.output_backend = backend
    show(plot)
