# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2017 Laboratory of Computational Systems Biotechnology (LCSB),
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
from collections import OrderedDict
from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.layouts import column
from bokeh.palettes import Spectral11, viridis

import numpy as np
import pandas as pd


def timetrace_plot(time,data,filename='out.html', legend = None):
    """
    Classic time vs. Y-value plot.

    :param time:
    :param data:
    :param filename:
    :param legend:
    :return:
    """
    # Make cool plot functions maype there is also a cooler way?
    if filename == '':
        # show the plot
        pass
    else:
        # print to file
        # output to static HTML file
        output_file(filename)

        # MAKE COLORZ
        num_species = data.shape[1]
        if num_species > 11:
            palette = viridis(num_species)
        else:
            palette = Spectral11[:num_species]

        # PLOTS
        p = figure(x_axis_type="datetime")
        
        for e in range(num_species):
            if legend is not None:
                legend_str = legend[e]
            else:
                legend_str = 'Compound {}'.format(e)
            p.line(time,data[:,e], line_color = palette[e], legend=legend_str)

        # show the results
        show(p)


def boxplot(df, filename):

        mean = df.mean()
        not_nan = [i for i,e in enumerate(mean) if e is not np.nan]
        cats = df.mean().dropna().index.values
        # find the quartiles and IQR for each category
        groups = df[cats]
        q1 = groups.quantile(q=0.25)
        q2 = groups.quantile(q=0.5)
        q3 = groups.quantile(q=0.75)
        iqr = q3 - q1
        upper = q3 + 1.5 * iqr
        lower = q1 - 1.5 * iqr

        qmin = groups.quantile(q=0.00)
        qmax = groups.quantile(q=1.00)
        upper = [min([x, y]) for (x, y) in zip(list(qmax.loc[:]), upper)]
        lower = [max([x, y]) for (x, y) in zip(list(qmin.loc[:]), lower)]

        p = figure(tools="",
                   plot_height=1000,
                   plot_width=20*len(cats),
                   y_axis_type="log",
                   background_fill_color="#efefef",
                   y_range=(qmin.min(), qmax.max()),
                   x_range=cats,
                   toolbar_location=None)



        # stems
        p.segment(cats, upper, cats, q3, line_color="black")
        p.segment(cats, lower, cats, q1, line_color="black")

        # boxes
        p.vbar(cats, 0.7, q2, q3, fill_color="#E08E79", line_color="black")
        p.vbar(cats, 0.7, q1, q2, fill_color="#3B8686", line_color="black")

        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = "white"
        p.grid.grid_line_width = 2
        p.xaxis.major_label_text_font_size = "12pt"
        p.xaxis.major_label_orientation = np.pi / 2.

        output_file(filename)
        show(p)



def plot_population_per_variable(data, filename, stride = 1):
    """

    :param data:
    :param filename:
    :param stride: How many points to skip for the plot. `stride=100` will only plot every 100
    points
    :return:
    """

    plots = OrderedDict()

    grouped = data.groupby('solution_id')
    # colors = viridis(len(grouped))

    for var in data.columns:
        if var in ['solution_id','time']:
            continue

        p = figure(x_axis_type = 'datetime')

        for group,this_data in grouped:
            this_data = this_data.iloc[::stride]
            p.line(this_data['time'], this_data[var],
                   # line_color = colors[group],
                   line_alpha = 0.2)

        p.title.text = var

        plots[var] = p

        output_file(filename.format(var))
        show(p)
        curdoc().clear()

    # c = column([p for p in plots.values()])
    # show(c)
