"""

This Python code accompanies the review article by Brumm, Krause, Schaab, & Scheidegger (2021)
and corresponds to the International Real Business Cycle (IRBC) model with irreversible invest-
ment, as presented in section 3.2.

We demonstrate the capabilities of adaptive sparse grids as these occasionally binding constraints
induce non-differentiabilities in policy functions, a challenge for approximation that adaptivity is
able to address.

Then, we generate a sparse grid as well as adaptive sparse grid approximations to the function
and measure their accuracy as follows: We randomly generate 10,000 uniformly distributed test
points, and compute the maximum error and the L1-error.

To do so, we use TASMANIAN library to showcase the general workflow of setting up a sparse grid
and approximate a multivariate function.

--------------------------------------------------------------------------------------------------

The parameters that govern the adaptivity of grid points are located in setup_asg.py. In particular,
the following parameters can be set:

- surplThreshold is the threshold that determines whether new grid points are added.

- maxRef is the number of maximum refinements of the sparse grid. A value of 0 corresponds to an
    ordinary grid without any refinement.

- dimRef determines the model outputs that are considered in the refinement process. A value of -1
    implies that all outputs are used.

- typeRefinement corresponds to the refinement strategy. The following options are available:

    ['classic', 'parents', 'direction', 'fds', 'stable'] (see TASMANIAN documentation for details)

- scaleCorr is a scale correction applied before the grid refinement. This implies that model outputs
    can be weighted in their contribution to the grid refinement. Note that this has to be applied to
    all existing points in the grid such that scaleCorr has to be an array of the same size as the grid.
    A value of 0 then implies that a given model output is ignored in the refinement process.

--------------------------------------------------------------------------------------------------

The code is organized as follows:

    1) Solution of the IRBC model with irreversible investment, using either an ordinary or an
        adaptive sparse grid (given the parametrization in setup_asg).

    2) Computation of max- and L1-errors along the simulation path.

    3) Computation of max- and L1-errors for 10,000 points in the state space.


"""

################################################################################
#                        Load necessary libraries                              #
################################################################################

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import copy

from parameters import *
from setup_num_int import *
from setup_asg import *
import time_iteration as time_iter
import postprocessing as post

import Tasmanian

import argparse

# TeX Support for plots
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Parser settings for terminal calls
parser = argparse.ArgumentParser(description='Run the IRBC model.')
parser.add_argument('--final_grid', dest='load_flag', action='store_true',\
                        help='Postprocessing only')

args = parser.parse_args()
load_flag = args.load_flag


################################################################################
#                                  Main Loop                                   #
################################################################################

if args.load_flag:

    print('Postprocessing only.')
    
    # Construct the grid structure
    gridFinal = Tasmanian.TasmanianSparseGrid()
    # Read the properties from grid_final.txt
    gridFinal.read(data_location + "grid_final.txt")
    # Set as interpolation grid
    grid1 = Tasmanian.copyGrid(gridFinal)

else:

    print('Start time iteration from iter = ', numstart)

    # If numstart==0 => start from scratch; numstart>0 => Restart from given iteration

    if (numstart>0):

        # Construct a new grid structure
        gridRestart = Tasmanian.TasmanianSparseGrid()
        # Read properties from file for restart
        gridRestart.read(data_location + "grid_iter_" + str(numstart) + ".txt")
        # Set as interpolation grid
        grid0 = Tasmanian.copyGrid(gridRestart)

    for iter0 in range(numstart,maxiter+1):

        polGuess1 = copy.copy(polGuess)

        grid1 = time_iter.fresh_grid()

        # Index of current grid level to control the number of refinements
        ilev = gridDepth

        while ((grid1.getNumNeeded() > 0) and (ilev<=maxRefLevel)):

            grid1 = time_iter.ti_step(grid1,polGuess1,grid0)

            # We start the refinement process after a given number of iterations
            if (iter0>iterRefStart):
                grid1,polGuess1 = time_iter.refine(grid1)

            # Track the grid level
            ilev += 1

        ## Calculate (approximate) errors on tomorrow's policy grid
        metric, polGuess, grid0 = time_iter.policy_update(grid0,grid1)

        print("Iteration: %2d, Grid pts: %2d, Level: %2d, Metric: %.4E" % (iter0, grid0.getNumPoints(),ilev,metric))

        if (np.mod(iter0+1,savefreq)==0):
            time_iter.save_grid(grid0,iter0)

        if (metric<tol_ti):
            break

    grid1.write(data_location + "grid_final.txt")

    error_sim = post.errors_sim(grid1)
    errors_ss = post.errors_ss(grid1)


post.plot_policies(1,grid1)
