"""
This file computes the residuals of the system of equilibrium conditions.

Inputs:
1) x -- array that contains the (guessed) values for all policies
2) state -- array that contains the values of the state variables
3) grid -- grid structure that contains the function values for interpolation

Outputs:
1) res -- array of size [nCountries*2+1] that contains the residuals of the
    system of equilibrium conditions.

"""

import numpy as np
from Expect_FOC import *
from aux_fcts import *


def sysOfEqs(x,state,grid):

    # State variables
    capStates = state[0:nCountries]
    tfpStates = state[nCountries:]

    # Policy values
    capPolicies = x[0:nCountries]
    lamb = x[nCountries]

    if typeIRBC=='non-smooth':

        gzAlphas = x[nCountries+1:]

        # Garcia-Zengwill transformation of the occasionally binding constraints
        gzAlphaPlus = np.maximum(0.0,gzAlphas)
        gzAlphaMinus = np.maximum(0.0,-gzAlphas)

    # Computation of integrands
    Integrands = ExpectFOC(capPolicies, state, grid)

    IntResult = np.empty(nCountries)

    for iint in range(nCountries):
        IntResult[iint] = np.dot(IntWeights,Integrands[:,iint])

    res = np.zeros(nPols)

    # Computation of residuals of the equilibrium system of equations

    if typeIRBC=='non-smooth':

        # Euler equations & GZ alphas
        for ires in range(nCountries):
            res[ires] = (betta*IntResult[ires] + gzAlphaPlus[ires])\
                            /(1.0 + AdjCost_ktom(capStates[ires],capPolicies[ires])) - lamb
            res[nCountries+1+ires] = capPolicies[ires] - capStates[ires]*(1.0-delta) - gzAlphaMinus[ires]

    else:

        # Euler equations
        for ires in range(nCountries):
            res[ires] = betta*IntResult[ires]/(1.0 + AdjCost_ktom(capStates[ires],capPolicies[ires])) - lamb


    # Aggregate resource constraint
    for ires2 in range(nCountries):
        res[nCountries] += F(capStates[ires2],tfpStates[ires2]) + (1.0-delta)*capStates[ires2] - capPolicies[ires2]\
                            - AdjCost(capStates[ires2],capPolicies[ires2]) - (lamb/pareto)**(-1.0/gamma)



    return res
