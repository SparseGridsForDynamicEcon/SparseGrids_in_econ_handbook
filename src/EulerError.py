from parameters import *
import numpy as np
import Tasmanian
from aux_fcts import *
from setup_num_int import *
from Expect_FOC import *

def Euler_error(x,state,pols):

    capStates = state[0:nCountries]
    tfpStates = state[nCountries:]

    capPolicies = x[0:nCountries]
    lamb = x[nCountries]
    gzAlphas = x[nCountries+1:]

    gzAlphaPlus = np.maximum(0.0,gzAlphas)
    gzAlphaMinus = np.maximum(0.0,-gzAlphas)

    Integrands = np.empty((numNodes,nCountries))
    newstate = np.empty(nCountries)
    captomtom = np.empty(nCountries)
    MPKtom = np.empty(nCountries)
    gzAlpTomPl = np.empty(nCountries)

    #Compute Density
    density = np.pi**(-(nCountries+1) * 0.5)

    for i_int in range(numNodes):

        for icou in range(nCountries):
            newstate[icou] = rhoZ*tfpStates[icou] + (GHnodes[i_int,icou] + GHnodes[i_int,nCountries])

        state_Sim = np.empty(nCountries*2)
        state_Sim[0:nCountries] = capPolicies
        state_Sim[nCountries:] = newstate

        lambtom = pols.evaluate(state_Sim)[nCountries]

        for icou in range(nCountries):
            captomtom[icou] = pols.evaluate(state_Sim)[icou]
            gzAlpTomPl[icou] = np.maximum(0.0,pols.evaluate(state_Sim)[nCountries+1+icou])
            MPKtom[icou] = 1.0 - delta + Fk(capPolicies[icou],newstate[icou]) - AdjCost_k(capPolicies[icou],captomtom[icou])
            Integrands[i_int,icou] = (lambtom*MPKtom[icou] - (1.0-delta)*gzAlpTomPl[icou]) * density


    IntResult = np.empty(nCountries)

    for iint in range(nCountries):
        IntResult[iint] = np.dot(GHweights,Integrands[:,iint])

    res = np.zeros(nCountries)

    for ires in range(nCountries):
        res[ires] = betta*IntResult[ires]/(lamb*(1.0 + AdjCost_ktom(capStates[ires],capPolicies[ires]))) - 1.0


    return res
