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

import logging
import os
import time

LOGFOLDERNAME = 'logs'

def get_bistream_logger(name):
    """
    Sets up a logger that outputs INFO+ messages on stdout and DEBUG+ messages
    in the log file

    :param name: a class __name__ attribute
    :return:
    """
    logger = logging.getLogger(name)

    try:
        os.makedirs(LOGFOLDERNAME)
    except FileExistsError:
        pass

    if not logger.handlers:
        logger.setLevel(logging.DEBUG)

        # create a file handler
        filename = name+get_timestr()+'.log'
        path = os.path.join(LOGFOLDERNAME, filename)
        file_handler = logging.FileHandler(path)
        file_handler.setLevel(logging.DEBUG)

        # create a stream handler
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)

        # create a logging format
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        stream_handler.setFormatter(formatter)

        # add the handlers to the logger
        logger.addHandler(file_handler)
        logger.addHandler(stream_handler)

    return logger

def get_timestr():
    timestr = time.strftime("_%Y%m%d_%H%M%S")
    return timestr