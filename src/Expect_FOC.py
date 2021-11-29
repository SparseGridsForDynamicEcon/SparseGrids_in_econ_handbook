import numpy as np
from parameters import *
from setup_num_int import *
from aux_fcts import *

def ExpectFOC(ktemp, state, grid):

    # 1) Determine next period's tfp states

    newstate = np.empty((numNodes,nCountries))

    for itfp in range(nCountries):
        newstate[:,itfp] = rhoZ*state[nCountries+itfp] + (GHnodes[:,itfp] + GHnodes[:,nShocks-1])
        newstate[:,itfp] = np.where(newstate[:,itfp] > aMin, newstate[:,itfp], aMin)
        newstate[:,itfp] = np.where(newstate[:,itfp] < aMax, newstate[:,itfp], aMax)

    # 2) Determine next period's state variables

    evalPt = np.empty((numNodes,nCountries*2))

    evalPt[:,0:nCountries] = ktemp
    evalPt[:,nCountries:] = newstate


    # 3) Determine relevant variables within expectations operator

    capPrPr = grid.evaluateBatch(evalPt)[:,0:nCountries]
    lambPr = grid.evaluateBatch(evalPt)[:,nCountries]
    gzAlphaPr = grid.evaluateBatch(evalPt)[:,nCountries+1:]

    gzAplusPr = np.maximum(0.0,gzAlphaPr)

    #Compute MPKtom
    MPKtom = np.empty((numNodes,nCountries))
    for impk in range(nCountries):
        MPKtom[:,impk] = 1.0 - delta + Fk(ktemp[impk],newstate[:,impk]) - AdjCost_k(ktemp[impk],capPrPr[:,impk])

    #Compute Density
    density = np.pi**(-(nCountries+1) * 0.5)

    #Specify Integrand
    ExpectFOC = np.empty((numNodes,nCountries))

    for iexp in range(nCountries):
        ExpectFOC[:,iexp] = (MPKtom[:,iexp]*lambPr - (1.0-delta)*gzAplusPr[:,iexp]) * density


    return ExpectFOC
