# Rippled shock simulation

### Quick start

The simulation itself is in three parts: rippled_shock.py, simulation_functions.py and simulation_parameters.py. The jupyter notebook has the plotting codes and .hdf5 file is an example file with simulation results.

Change the parameters in the simulation_parameters.py and run the rippled_shock.py for the simulation. If you have a strict anti-malware protecion on your computer it might stop the simulation since the parallelization of the simulation spawns new processes as needed. The simulation is poorly optimized and with large particle counts (> 2 million, on a fairly modern laptop) RAM will become the limiting factor.

To plot the results read the jupyter notebook.

## User manual

### Prerequisites

The code is in Python 3 and to run the simulation an up-to-date python 3 package is recommended. Required external libraries are: [numpy](https://numpy.org/), [numba](http://numba.pydata.org/numba-doc/latest/index.html#) and [h5py](https://www.h5py.org/). The base python 3 should include the rest of the needed libraries: os, time, random and warnings. 

To use the plotting code you will need the [jupyter notebook](https://jupyter.org/).

For beginners the easiest way to start is to install [Anaconda](https://www.anaconda.com/). This package manager has everything you need and is easy to keep up-to-date. Just remember to download the Python 3 version.

### Download the code

From the top of this Github-page you can download the rippled shock simulation code. If you just want the simulation you need the files: simulation_parameters.py, simulation_functions.py and rippled_shock.py. 

The jupyter notebook (the .ipynb file) and the example result file (the .hdf5 file) are only used for plotting.

### Test the code

To run simulations you only need to change the parameters in the simulation_parameters.py file. Do not change the rippled_shock.py or simulation_functions.py unless you know what you are doing.

To run the simulation just run the rippled_shock.py file.

Before running any real simulations I recommend you test the code with a small simulation of 10 particles and 10.0 simulation time. Even with the most ancient computer this should be done in under 5 seconds.

If this is not the case and the code never finishes then the anti-malware on your computer might be supressing the code. The parallelization on the simulation is done using the multiprocessing library and some anti-malware programs understandably don't like processes spawning more processes. If this is the case just google "Add Exception to ..." (Windows Defender/Mac firewall/etc.) and follow the instructions.

If the test run works correcly then you can continue to larger simulations. Just be careful to not run too many particles in one simulation. The writing to file is done poorly and running too many particles in one simulation will crash the simulation once the RAM is filled up. To my knowledge there is no upper limit to the simulation time.

### Tips

- Don't change the file name each time. Learn to use the groups in hdf5-format and life (or atleast plotting these simulation results) will be much easier.

- To know how many cores you computer has just open python and write:
```python
import os
print(os.cpu_count())
```

- If you accidentally run too many particles and the simulation crashes the next time you try to open the hdf5-file you will get an EOF (end of file) error and the whole file is corrupted. Write a new file when testing how many particles you can run.


### Plotting the results

Read the .ipynb to learn how to plot the results.
