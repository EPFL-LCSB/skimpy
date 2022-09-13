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

Run the 

