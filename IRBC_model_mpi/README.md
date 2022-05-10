# The IRBC Benchmark Model

This script provides code, parallelized with [MPI4PY](https://mpi4py.readthedocs.io/en/stable/), to solve the benchmark model in section 3 of the review article paper by [Brumm, Krause, Schaab, & Scheidegger (2022)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3979412). The simple parallelization approach we follow here adopts ideas from [Brumm & Scheidegger (2017), Section 4](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA12216), and can, in the case of large sparse grids, speed up the time-to-solution by orders of magnitude.

### Prerequisites / Installation

To run the benchmark IRBC model in parallel, it is necessary that Tasmanian as well as mpi4py is installed. Tasmanian is included in the
Python Pip index: https://pypi.org/project/Tasmanian/

```shell
$ python3 -m pip install scikit-build packaging numpy --user (required dependencies)
$ python3 -m pip install Tasmanian --user                    (user installation)
$ python3 -m pip install Tasmanian                           (virtual env installation) 
$ python3 -m pip install pip install mpi4py                  
```
Further information on alternative installation procedures can be found here: https://tasmanian.ornl.gov/documentation/md_Doxygen_Installation.html

## Usage
There are two modes to run this code:

   1. start the time iteration process from scratch (or restart from a predetermined iteration),
   2. inspect the final results without computing the model solution.

We have simplified the code such that the only user input is the desired running mode. To run, follow
these instructions.

Access this folder in the terminal:

```shell
    $ cd <PATH to the repository>/SparseGrids_in_econ_handbook/IRBC_model
```

### Mode 1: Start from scratch (or restart from a predetermined iteration)
If you want to start the time iteration from scratch or restart from a predetermined iteration, use the
following command in the terminal:

```shell
    $ python mainASG.py
```
In parameters.py, there are two parameters that are relevant here:

a) `numstart` determines at which iteration the main loop starts

b) `savefreq` determines the frequency at which the grid structure is stored

The default setting for savefreq is 10 such that the grid structure from every 10th iteration is stored. 
We provide these in `./data`. A list of .py-files that are involved in the iterative process can be found below.


#### MPI

The solutions during time-iteration can be embarassingly parallelized via MPI. To use this,
make sure mpi4py is set up, `export OMP_NUM_THREADS=1`, and then run (e.g. with 4 phyiscal cores):

```
mpirun -np 4 python mainASG.py
```

Depending on the particular installation of MPI, one sometimes has to slightly alter the run command as follows to avoid performance bottlenecks:
```
mpirun --bind-to core -np 4 python mainASG.py
```
The flag [--bind-to core](https://www.open-mpi.org/doc/v3.0/man1/mpirun.1.php) binds processes to cores.


### Mode 2: Postprocessing only
If you want to run the postprocessing without computing the model solution, you can use the following
command in the terminal:

```shell
    $ python mainASG.py --final_grid
```
The grid structure after the last time iteration step can be found in `./data` and will be loaded directly.


## Model Versions

This code solves the IRBC model both with and without irreversible investment. In parameters.py, there is the 
variable `typeIRBC` which can take the values `non-smooth` and `smooth`, indicating the respective model version.

There are two folders in `./data` containing plots and .txt-files of policy functions as well as stored grid
structures in .txt-files for both model versions.


## Code Organisation

The code on the IRBC model consists of the following files:

- mainASG.py        -- the main iterative loop; executes the postprocessing

- time_iteration.py -- contains time iteration step as well as the grid refinement process

- postprocessing.py -- contains postprocessing functions such as error measures and storage of results

- parameters.py     -- contains economic and technical parameters

- setup_asg.py      -- contains parameters and construction of the adaptive sparse grid structure

- setup_num_int.py  -- constructs nodes and weights for numerical integration

- ExpectFOC.py      -- constructs the integrands for the system of equations

- SOE.py            -- contains the system of equations

- aux_fcts.py       -- contains auxilliary functions


## Postprocessing

After convergence of the time iteration, the following postprocessing steps are executed:
    
- Computation of average and maximum errors along the simulation path
   
- Computation of average and maximum errors in the state space
    
- Plotting the policy functions in 1D
    
- Storing the policy functions


