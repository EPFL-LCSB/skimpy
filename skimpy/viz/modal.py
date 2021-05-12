# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2020 Laboratory of Computational Systems Biotechnology (LCSB),
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
from math import pi
from collections import OrderedDict
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import BasicTicker, ColorBar, LinearColorMapper, PrintfTickFormatter
from bokeh.layouts import column

from bokeh.io import export_svgs

from scipy.cluster import hierarchy
import numpy as np
import pandas as pd

COLORS = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce",
          "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]


def plot_modal_matrix(M, filename='map.html',
                      phase=False,
                      colors=COLORS,
                      clustered=False,
                      x_label='',
                      y_label='',
                      backend='webgl',
                      **kwargs):

    # output to static HTML file
    output_file(filename)

    if phase:
        # calculate phase
        data = M.phi()
    else:
        data = M.abs()

    if clustered:
        # Cluster metabolites
        Z = hierarchy.linkage(data.T, 'complete', optimal_ordering=True)
        d = hierarchy.dendrogram(Z, no_plot=True)
        ix_col = np.array(d['ivl'], dtype=np.int)
        # Cluster eigenvalues
        Z = hierarchy.linkage(data, 'complete', optimal_ordering=True)
        d = hierarchy.dendrogram(Z, no_plot=True)
        ix_index = np.array(d['ivl'], dtype=np.int)
        data = data.iloc[ix_index, ix_col]

    TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

    # Format
    _lambda = ['{0:.2e} {1} {2:.2e}i'.format(np.real(n), '+-'[np.imag(n) < 0], abs(np.abs(np.imag(n))))
               for n in data.index]
    data.index = _lambda
    metabolites = list(data.columns)

    df = pd.DataFrame(data.stack(), columns=['value']).reset_index()
    df.columns = ['lambda', 'metabolite', 'value']

    mapper = LinearColorMapper(palette=colors, low=0, high=1.0)

    # PLOT
    p = figure(x_axis_location="above",
               tools=TOOLS,
               toolbar_location='below',
               x_range=metabolites,
               y_range=_lambda,
               **kwargs)

    p.rect(x="metabolite", y='lambda', width=1, height=1,
           source=df,
           fill_color={'field': 'value', 'transform': mapper},
           line_color=None)

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="7px",
                     ticker=BasicTicker(desired_num_ticks=len(colors)),
                     label_standoff=6, border_line_color=None, location=(0, 0))

    p.add_layout(color_bar, 'right')
    p.xaxis.major_label_orientation = pi / 2

    # Make the axsis pretty
    p.xaxis.axis_label = x_label
    p.yaxis.axis_label = y_label

    p.output_backend = backend
    show(p)
