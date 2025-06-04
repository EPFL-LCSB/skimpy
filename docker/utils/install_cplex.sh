#!/bin/sh
# Thanks to cdiener: https://hub.docker.com/r/cdiener/cobra-docker/~/dockerfile/
# For the solution of simply getting the bins and python hooks

#
echo "Installing and Moving CPLEX files"

## Default Py3.9 install
if [ -d /solvers/ibm ]; then cd /solvers/ibm/ILOG/CPLEX_Studio221/cplex/python/3.9/x86-64_linux/ && \
	python3 setup.py install && \
	cp /solvers/ibm/ILOG/CPLEX_Studio221/cplex/bin/x86-64_linux/cplex /usr/bin/; fi

