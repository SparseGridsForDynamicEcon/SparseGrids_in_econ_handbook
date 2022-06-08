"""
This file sets the relevant parameters for the adaptive sparse grid and initializes
the grid structure that is used in the main file. It also includes parameters that
govern the adaptivity process.

"""

import numpy as np
import Tasmanian
from parameters import *
from scipy import optimize


################################################################################
#                               Grid construction                              #
################################################################################

# Number of dimensions (capital stock and tfp for each country)
gridDim = nCountries*2
# Number of outputs (capital policy & multiplier for each country + ARC multiplier)
gridOut = nPols
# Grid level (we start with a sparse grid of level 4)
gridDepth = 3
# 1=linear, 2=quadratic, 3=cubic
gridOrder = 1
# Type of base functions
gridRule = "localp"

# Set the grid domain to [kmin,kmax]^n x [amin,amax]^n
gridDomain = np.zeros((gridDim,2))

gridDomain[0:nCountries,0] = kMin
gridDomain[0:nCountries,1] = kMax
gridDomain[nCountries:,0] = aMin
gridDomain[nCountries:,1] = aMax

# Generate the grid structure
grid0 = Tasmanian.makeLocalPolynomialGrid(gridDim,gridOut,gridDepth,gridOrder,gridRule)
# Transform the domain
grid0.setDomainTransform(gridDomain)
# Get the points that require function values
aPoints = grid0.getPoints()
# Get the number of points that require function values
aNum = grid0.getNumPoints()


################################################################################
#                        Adaptivity parameters                                 #
################################################################################

# Surplus threshold
surplThreshold = 1e-3
# Number of maximum refinements
maxRef = 1
# Maximum Level of ASG
maxRefLevel = gridDepth + maxRef
# Outputs that are considered in the refinement process (-1 implies that all outputs are considered)
dimRef = -1
# Type of grid refinements
typeRefinement = 'classic'
# Scale correction in the refinement process:
# We only let the capital policies of each country determine the addition of grid points
scaleCorr = np.zeros(nPols)
scaleCorr[0:nCountries] = 1


################################################################################
#                             Initialization                                   #
################################################################################

polGuess = np.zeros((aNum,nPols))

# Guesses for the capital policies and investment constraint multipliers
for i0 in range(nCountries):
    polGuess[:,i0] = aPoints[:,i0]*(1.0-delta)

    if typeIRBC=='non-smooth':
        polGuess[:,nCountries+1+i0] = -delta*aPoints[:,i0]

# Guess for the aggregate ressource constraint multiplier
# We use a nonlinear equation solver to find the ARC multipliers that are consistent
# with the guesses for the capital policies and investment constraint multipliers

def ARC_zero(lam_gues):
    
    res = 0.0
    
    for i1 in range(nCountries):
        res += np.exp(aPoints[i0,nCountries+i1])*A_tfp*aPoints[i0,i1]**zeta - (-delta*kappa/2.0)**2 - (lam_gues/pareto[i1])**(-gamma[i1])
    
    return res

for i0 in range(aNum):
    root = optimize.root(ARC_zero, 0.1, method='lm')
    polGuess[i0,nCountries] = root.x
    # Print the status in case the solver did not find the root
    if root.success!= 1:
        print(root.message)


# Load function values into grid structure
grid0.loadNeededPoints(polGuess)
