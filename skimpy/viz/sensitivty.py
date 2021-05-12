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
from bokeh.models import ColumnDataSource, FactorRange
from bokeh.io import curdoc, show
from bokeh.transform import factor_cmap

def plot_sobol_coefficients(Si,
                            St,
                            filename='map.html',
                            colors=['#3399FF','#205f99'],
                            backend='webgl',
                            width=0.9,
                            x_label='',
                            y_label='',
                            **kwargs):
    """
    Plot sensitivity coefficients
    :param Si: pd.Series
    :param St: pd.Series
    :param filename: String
    :param colors:
    :param backend:
    :param x_label:
    :param y_label:
    :param kwargs:
    :return:
    """


    # output to static HTML file
    output_file(filename)
    index = ['Si','St']
    x = [ (str(v),i) for v in St.index for i in index]
    counts = sum(zip( Si.values, St.values), () )
    source = ColumnDataSource(dict(x=x,
                                   counts=counts))

    plot = figure( x_range=FactorRange(*x),
                   **kwargs,)

    plot.vbar(x="x",
              top="counts",
              width=width,
              fill_color=factor_cmap('x', palette=colors, factors=index, start=1, end=2),
              source=source)

    # Make the axsis pretty
    plot.xaxis.axis_label = x_label
    plot.yaxis.axis_label = y_label

    plot.y_range.start = 0
    plot.x_range.range_padding = 0.1

    plot.output_backend = backend
    show(plot)
