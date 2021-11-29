"""
This file sets the relevant parameters for the numerical integration and constructs
nodes and weights that are used to compute expextations.

"""

import numpy as np
import Tasmanian
from parameters import *


# Number of dimensions (One shock per country and the aggregate shock)
int_dim = nCountries+1
# Number of outputs
int_out = 1
# Level of integration grid (Number of points in each dimension; starts at 0)
int_depth = 2
# Tensor selection strategy
int_type = "tensor"
# Integration rule
int_rule = "gauss-hermite"

# Generate the grid structure
int_grid = Tasmanian.makeGlobalGrid(int_dim, int_out, int_depth, int_type, int_rule)
# Transform the nodes accordingly
GHnodes = np.sqrt(2)*sigE*int_grid.getNeededPoints()
# Get the corresponding quadrature weights
GHweights = int_grid.getQuadratureWeights()
# Get the number of nodes
numNodes = len(GHnodes)
