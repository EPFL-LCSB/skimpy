# Install CPLEX 

Download the binaries CPLEX Studio version compatible with your operating system and the python version you used for the SKiMpy installation from the [IBM](https://www.ibm.com/academic/home) website (Free academic version avialable).
We recommend to use CPLEX Studio 221 as the package was tested in combindation with this version of CPLEX.

## OSx 

Unzip the archive an launch the installer. Following the installation instructons. 
After installing CPLEX Studio you need to install the python API for the respective virtual environment.

```bash
conda activate skimpy-env
python /Applications/CPLEX_Studio221/python/setup.py install
```

## Linux/WSL 

Execute the binary file you downloaded from the IBM website and follow the installtion instructions. 
Then navigate to the installation directory of the CPLEX Solver Suite, `path/to/yourCplexHome` and install the python API for the respective virtual environment.

```bash
cd path/to/yourCplexhome/python/PYTHON-VERSION/PLATFORM
conda activate skimpy-env
python setup.py install
```
