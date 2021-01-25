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
from collections import OrderedDict
from bokeh.plotting import figure, output_file, show, curdoc, ColumnDataSource
from bokeh.layouts import column
from bokeh.palettes import Spectral11, viridis
from bokeh.io import export_svgs

import numpy as np
import pandas as pd


def plot_flux_profiles(data,
                       filename='out.html',
                       min_max=None,
                       colors=None, **kwargs):
    """

    :param data: pd.DataFrame
    :param filename: string
    :param min_max: tfa min/max output
    :return:
    """

    # output to static HTML file
    output_file(filename)

    # MAKE COLORS
    if colors is None:
        num_profiles = data.shape[0]
        if num_profiles > 11:
            colors = viridis(num_profiles)
        else:
            colors = Spectral11[:num_profiles]

    if not min_max is None:
        this_min_max = min_max.loc[data.columns]
        mid = np.arange(data.columns)
        left = mid - 0.2
        right = mid + 0.2

        p.quad(top=this_min_max['minimum'], bottom=this_min_max['maximum'], left=[1, 2, 3],
               right=[1.2, 2.5, 3.7], color="#B3DE69")


    source = ColumnDataSource(data)

    plot = figure(**kwargs)


