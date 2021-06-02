# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2019 Laboratory of Computational Systems Biotechnology (LCSB),
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

from skimpy.utils.tabdict import TabDict
from sympy import Symbol
from pandas import Series
import h5py

import numpy as np
import pandas as pd

class ParameterValues(object):
    """
    Parameters set for kinetic models wich can be indexed with symbols or
    """
    def __init__(self,parameter_values, kmodel=None):
        """

        :param kmodel: KineticModel class
        :param parameter_values: a dict contaning parameter names and values
        """
        if parameter_values.__class__ is Series:
            parameter_values = parameter_values.to_dict()

        self._parameter_values = TabDict([(str(p),v) for p,v in parameter_values.items()])

        # Check if this is a good solution

        if kmodel is not None:
            model_params =  kmodel.parameters
            self._str_to_sym = { p:model_params[p].symbol for p in model_params}
        else:
            self._str_to_sym = { p:Symbol(p) for p in parameter_values}

        self._sym_to_str = {self._str_to_sym[p]:p for p in self._parameter_values}

        #self._str_to_param = {p:model_params[p] for p in self._parameter_values}

    def __getitem__(self, item):
        if item.__class__ is  str:
            return self._parameter_values[item]
        if item.__class__ is Symbol:
            return self._parameter_values[self._sym_to_str[item]]

    def __setitem__(self, item, value):
        if item.__class__ is  str:
                self._parameter_values[item] = value
        if item.__class__ is Symbol:
                self._parameter_values[self._sym_to_str[item]] = value

    def items(self):
        return self._parameter_values.items()

    def keys(self):
        return  self._parameter_values.keys()

    def values(self):
        return  self._parameter_values.values()


class ParameterValuePopulation(object):
    def __init__(self,data, kmodel=None, index=None):

        self.kmodel = kmodel

        if type(data) == list:
            # Todo check for indexable
            self._data = [ ParameterValues(d,kmodel=kmodel) for d in data]
            if index is None:
                self._index = TabDict((str(i),i) for i,_ in enumerate(data))
            else:
                self._index = TabDict((k,i) for i,k in enumerate(index))

        elif type(data) == pd.DataFrame:
            raise NotImplemented("Type {} is not supported yet in progress".format(type(data)))

        elif type(data) == ParameterValues:
            # Todo check for indexable
            self._data = [ d for d in data]
            if index is None:
                self._index = TabDict((str(i),i) for i,_ in enumerate(data))
            else:
                self._index = TabDict((k,i) for i,k in enumerate(index))

        else:
            raise TypeError("Type {} is not supported".format(type(data)))

    def __getitem__(self,  index):
        if self._index is None:
            return self._data[index]
        else:
            return self._data[self._index[index]]

    def __len__(self,):
        return len(self._data)

    # Define the iterator
    def __iter__(self):
        self.n = 0
        return self

    def __next__(self):
        if self.n < len(self._data):
            result = self._data[self._index[self.n]]
            self.n += 1
            return result
        else:
            raise StopIteration


    def mean(self):
        """
        :return Computes the mean parameter values for the population:
        """
        if self._index is None:
            this_mean = pd.DataFrame(data= [dict(self._data[i]) for i in self._data.keys() ]).mean()
        else:
            this_mean = pd.DataFrame(data=[dict(self._data[self._index[i]]) for i in self._index ]).mean()

        return ParameterValues(this_mean, kmodel)


    def save(self,filename):
        """
        Saves the parameter population as hdf5 file
        :param filename: string XXX.h5 / XXX.hdf5
        :return:
        """
        f = h5py.File(filename, 'w') #TODO catch existing file?

        # TODO more central way ?
        param_names = np.array([k for k,v  in self._data[0]._parameter_values.items() if not (v is None)],
                               dtype=object)

        string_dt = h5py.special_dtype(vlen=str)

        f.create_dataset('parameter_names', data=param_names, dtype=string_dt)
        f.create_dataset('num_parameters_sets', data=len(self._data))
        f.create_dataset('index', data=np.array([k for k in self._index],dtype=object),
                         dtype=string_dt )


        for i,this_data in enumerate(self._data):
            this_data = np.array([this_data._parameter_values[p]
                                 for p in param_names], dtype=np.float64)
            f.create_dataset('parameter_set_{}'.format(i), data=this_data)

        f.close()


## TODO Lets see this should maybe
def load_parameter_population(filename, lower_index=None, upper_index=None):
    f = h5py.File(filename, 'r')
    data = []
    if lower_index is None:
        lower_index = 0
    if upper_index is None:
        upper_index = int(np.array(f.get('num_parameters_sets')))

    # deconde
    param_names = f.get('parameter_names')[:].astype(np.unicode_)

    try:
        index = f.get('index')[:].astype(np.unicode_)
        if index[0] is None:
            index = None
    except: # Put an error here
        index = None

    for i in range(lower_index,upper_index):
        this_param_set = 'parameter_set_{}'.format(i)
        param_values = f.get(this_param_set)[:].astype(np.float)
        this_data = {k:v for k,v in zip(param_names,param_values)}
        data.append(this_data)
    if index is None:
        param_population = ParameterValuePopulation(data)
    else:
        param_population = ParameterValuePopulation(data,
                                                    index=index[lower_index:upper_index])

    f.close()

    return param_population

