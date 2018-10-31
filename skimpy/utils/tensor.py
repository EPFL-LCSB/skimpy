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
from collections import deque, OrderedDict

import numpy as np
import pandas as pd

class Tensor(object):

    def __init__(self, data, indexes, *args, **kwargs):

        assert len(indexes) == 3
        for ix in indexes:
            assert isinstance(ix, pd.Index)

        self._data = data

        # Define the three axis indexes
        self.complementary_indexes = dict()
        self._i = indexes[0]
        self._j = indexes[1]
        self._k = indexes[2]

        self.build_index_dict([x.name for x in indexes])

    def build_index_dict(self, index_names):

        rotating_index_names = deque(index_names)

        self.indexes = OrderedDict(
            [(k,v) for k,v in zip(index_names, [self._i,
                                                self._j,
                                                self._k])])

        self.index_order = {k:v for k,v in zip(index_names, range(3))}

        # Use deque for circular permutation
        for k in index_names:
            index1 = self.indexes[rotating_index_names[-2]]
            index2 = self.indexes[rotating_index_names[-1]]
            self.complementary_indexes[k] = (index1, index2)
            rotating_index_names.rotate(-1)



    def slice_by(self, slicer, value):
        if isinstance(slicer,str):
            slicer = slicer.lower()
            slicer = self.indexes[slicer]
        else:
            slicer = self.indexes[slicer.name]

        index1, index2 = self.complementary_indexes[slicer.name]

        slicer_order = list(self.indexes.keys()).index(slicer.name)
        ix = slicer.get_loc(value)

        if slicer_order == 0:
            return pd.DataFrame(self._data[ix,:,:], index = index1, columns = index2)
        elif slicer_order == 1:
            return pd.DataFrame(self._data[:, ix,:].T, index = index1, columns = index2)
        elif slicer_order == 2:
            return pd.DataFrame(self._data[:,:,ix], index = index1, columns = index2)


if __name__ == '__main__':

    # Example: Control coefficient samples
    # 4 Samples * 3 parameters * 2 reactions

    sample_ix = pd.Index(range(4), name ='sample')
    param_ix = pd.Index(['Vmax1', 'Vmax2', 'km2'], name = 'param')
    rxn_ix = pd.Index(['r1','r2'], name = 'rxn')

    # Build the data
    data = np.ndarray(shape=(len(sample_ix), len(param_ix), len(rxn_ix)),
                      dtype = np.float32)

    sample0 = np.array( [[10, 20],
                        [30, 40],
                        [2 ,  7]] )

    for e,sample_no in enumerate(range(len(sample_ix))):
        data[sample_no,:,:] = sample0 + e/10.

    tdata = Tensor(data, indexes = [sample_ix, param_ix, rxn_ix])

    print('\n # Slice by sample:')
    print(tdata.slice_by('sample', 2))
    print('\n # Slice by param:')
    print(tdata.slice_by('param', 'Vmax2'))
    print('\n # Slice by reaction:')
    print(tdata.slice_by('rxn', 'r1'))





