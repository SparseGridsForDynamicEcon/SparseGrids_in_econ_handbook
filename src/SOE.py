import numpy as np
from Expect_FOC import *
from aux_fcts import *


def sysOfEqs(x,state,grid):

    capStates = state[0:nCountries]
    tfpStates = state[nCountries:]

    capPolicies = x[0:nCountries]
    lamb = x[nCountries]
    gzAlphas = x[nCountries+1:]

    gzAlphaPlus = np.maximum(0.0,gzAlphas)
    gzAlphaMinus = np.maximum(0.0,-gzAlphas)

    Integrands = ExpectFOC(capPolicies, state, grid)

    IntResult = np.empty(nCountries)

    for iint in range(nCountries):
        IntResult[iint] = np.dot(GHweights,Integrands[:,iint])

    res = np.zeros(nCountries*2+1)

    # Euler equations and irreversibility constraints

    for ires in range(nCountries):
        res[ires] = (betta*IntResult[ires] + gzAlphaPlus[ires])/(1.0 + AdjCost_ktom(capStates[ires],capPolicies[ires])) - lamb
        res[nCountries+1+ires] = capPolicies[ires] - capStates[ires]*(1.0-delta) - gzAlphaMinus[ires]

    # Aggregate resource constraint

    for ires2 in range(nCountries):
        res[nCountries] += F(capStates[ires2],tfpStates[ires2]) + (1.0-delta)*capStates[ires2] - capPolicies[ires2]\
                            - AdjCost(capStates[ires2],capPolicies[ires2]) - (lamb/pareto)**(-1.0/gamma)



    return res
