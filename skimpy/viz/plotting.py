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

def timetrace_plot(time,data,filename='out.html', legend = None):
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


def plot_population_per_variable(data, filename, stride = 1):

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
