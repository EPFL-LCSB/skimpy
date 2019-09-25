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

class TabDict(OrderedDict):
    """
    Really just an ordered dict with tab completion in interactive terminals
    """

    def __dir__(self):
        # Allow tab complete of items by their id
        attributes = dir(self.__class__)
        attributes.extend(self.keys())
        return attributes

    def __getattr__(self, attr):
        try:
            return TabDict.__getitem__(self, attr)
        except KeyError:
            raise AttributeError("TabDict has no attribute or entry %s" %
                                 attr)

    # TODO fix setting value by doing TabDict.key = value
    # def __setattr__(self, attr, value):
    #     try:
    #         if attr in self.keys():
    #             return TabDict.__setitem__(self, attr,value)
    #         else:
    #             raise KeyError()
    #     except KeyError:
    #         raise AttributeError("TabDict has no attribute or entry %s" %
    #                              attr)

    def iloc(self,ix):
        return list(self.items())[ix]



def iterable_to_tabdict(iterable, use_name = True):
    """
    Takes the items from an iterable and puts them in a TabDict, indexed by the
    elements' .name property

    :param iterable:
    :return:
    """
    if iterable is None:
        return TabDict()

    if use_name:
        return TabDict([(x.name, x) for x in iterable])
    else:
        return TabDict([(x.__str__(), x) for x in iterable])
