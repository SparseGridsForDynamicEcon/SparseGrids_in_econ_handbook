import numpy as np
import Tasmanian
from parameters import *

## Gauss-Hermite quadrature

GHnum = 3
GHdim = nCountries+1

# Parameters for numerical integration

int_dim = GHdim
int_out = 1
int_depth = GHnum-1
int_type = "tensor"
int_rule = "gauss-hermite"
int_grid = Tasmanian.makeGlobalGrid(int_dim, int_out, int_depth, int_type, int_rule)
GHnodes = np.sqrt(2)*sigE*int_grid.getNeededPoints()
GHweights = int_grid.getQuadratureWeights()

numNodes = len(GHnodes)
