"""
This file sets the relevant parameters for the numerical integration and constructs
nodes and weights that are used to compute expextations.

"""

import numpy as np
import Tasmanian
from parameters import *


if typeInt=='GH-quadrature':

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
    IntNodes = np.sqrt(2)*sigE*int_grid.getNeededPoints()
    # Get the corresponding quadrature weights
    IntWeights = int_grid.getQuadratureWeights()
    # Number of integration nodes
    numNodes = len(IntNodes)


elif typeInt=='monomials_2d':

    # Number of integration nodes
    numNodes = 2*nShocks
    # Pre-allocate nodes array
    z1 = np.zeros((numNodes,nShocks))
    # Fill nodes array with [1.0;-1.0]
    for i1 in range(nShocks):
        z1[i1*2,i1] = 1.0
        z1[i1*2+1,i1] = -1.0
    # Compute integration nodes
    IntNodes = z1*np.sqrt(nShocks)*sigE
    # Compute integration weights
    IntWeights = np.ones(numNodes)*1.0/numNodes


elif typeInt=='monomials_power':

    # Number of integration nodes
    numNodes = 2*nShocks**2 + 1

    z0 = np.zeros((numNodes,nShocks))

    # Deviations in one dimension (note that the origin is row zero)
    for i1 in range(nShocks):
        z0[i1*2+1,i1] = 1.0
        z0[i1*2+2,i1] = -1.0

    i0 = 0
    # Deviations in two dimensions
    for i1 in range(nShocks):
        for i2 in range(i1+1,nShocks):
            z0[2*nShocks+1+i0*4,i1] = 1.0
            z0[2*nShocks+2+i0*4,i1] = 1.0
            z0[2*nShocks+3+i0*4,i1] = -1.0
            z0[2*nShocks+4+i0*4,i1] = -1.0
            z0[2*nShocks+1+i0*4,i2] = 1.0
            z0[2*nShocks+2+i0*4,i2] = -1.0
            z0[2*nShocks+3+i0*4,i2] = 1.0
            z0[2*nShocks+4+i0*4,i2] = -1.0
            i0 += 1

    # Nodes
    IntNodes = np.zeros((numNodes,nShocks))
    IntNodes[1:nShocks*2+1,:] = z0[1:nShocks*2+1,:]*np.sqrt(2.0+nShocks)*sigE
    IntNodes[nShocks*2+1:] = z0[nShocks*2+1:]*np.sqrt((2.0+nShocks)/2.0)*sigE

    # Weights
    IntWeights = np.zeros(numNodes)

    IntWeights[0] = 2.0/(2.0+nShocks)
    IntWeights[1:nShocks*2+1] = (4-nShocks)/(2*(2+nShocks)**2)
    IntWeights[nShocks*2+1:] = 1.0/(nShocks+2)**2
    

else:

    print('Error in determining the type of numerical integration')
    exit()
