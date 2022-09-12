# SKiMpy Docker

This Docker offers a suitable environment to run SKiMpy.

## Requirements

Make sure [docker](https://www.docker.com/) is installed.
You might want to install a commercial solver. Gurobi and CPLEX are supported. See [solver/instructions.txt](https://github.com/EPFL-LCSB/pytfa/blob/master/docker/solvers/instructions.txt) for more informations on how to install them.

Note that SKiMpy requires a python version >= 3.7 the default docker version is currently 3.9. 
It is important that your solver supports your installed python version. We recommend installing CPLEX Studio221 which 
supports python versions: 3.7, 3.8, 3.9 and 3.10


## Running the Docker

First, build the container with `build.bat` or `. build`.
Then start the container with `run.bat` or `. run`.
```bash
. build
. run
```

You can run the examples in /skimpy/tutorials:
```bash
cd /skimpy/tutorials/metabolic/
python tutorial_metabolic.py

```

You can also run them inside IPython to experiment and play with the objects:

```bash
ipython
run tutorial_metabolic.py
kmodel.reactions
kmodel.reactants

```

## Additional information

If you are running your Docker container in a Unix-based environment, you might get permission errors on the `.sh` scripts.
This is because permissions are inherited from the host environment. You can fix this by running in the `docker` folder:
```bash
chmod +x utils/*.sh
```
