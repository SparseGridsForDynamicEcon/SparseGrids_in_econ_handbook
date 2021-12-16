"""
This file computes the expectation terms in the Euler equations of each country.

Inputs:
1) ktemp -- array that contains the (guessed) values for the capital policies
2) state -- array that contains the values of the state variables
3) grid -- grid structure that contains the function values for interpolation

Outputs:
1) ExpectFOC -- array of size [int_depth^(nCountries+1)]x[nCountries] that contains
    the expectation terms for each country in each possible state tomorrow

"""


import numpy as np
from parameters import *
from setup_num_int import *
from aux_fcts import *


def ExpectFOC(ktemp, state, grid):

    # 1) Determine next period's tfp states

    newstate = np.zeros((numNodes,nCountries))

    for itfp in range(nCountries):
        newstate[:,itfp] = rhoZ*state[nCountries+itfp] + (IntNodes[:,itfp] + IntNodes[:,nShocks-1])
        newstate[:,itfp] = np.where(newstate[:,itfp] > aMin, newstate[:,itfp], aMin)
        newstate[:,itfp] = np.where(newstate[:,itfp] < aMax, newstate[:,itfp], aMax)

    # 2) Determine next period's state variables

    evalPt = np.zeros((numNodes,nCountries*2))

    evalPt[:,0:nCountries] = ktemp
    evalPt[:,nCountries:] = newstate


    # 3) Determine relevant variables within the expectations operator

    capPrPr = grid.evaluateBatch(evalPt)[:,0:nCountries]
    lambPr = grid.evaluateBatch(evalPt)[:,nCountries]

    if typeIRBC=='non-smooth':
        gzAlphaPr = grid.evaluateBatch(evalPt)[:,nCountries+1:]
        gzAplusPr = np.maximum(0.0,gzAlphaPr)

    # Compute tomorrow's marginal productivity of capital
    MPKtom = np.zeros((numNodes,nCountries))
    for impk in range(nCountries):
        MPKtom[:,impk] = 1.0 - delta + Fk(ktemp[impk],newstate[:,impk]) - AdjCost_k(ktemp[impk],capPrPr[:,impk])


    #Compute Density
    if typeInt=='GH-quadrature':
        density = np.pi**(-(nCountries+1) * 0.5)
    else:
        density = 1.0

    #Specify Integrand
    ExpectFOC = np.zeros((numNodes,nCountries))

    if typeIRBC=='non-smooth':

        for iexp in range(nCountries):
            ExpectFOC[:,iexp] = (MPKtom[:,iexp]*lambPr - (1.0-delta)*gzAplusPr[:,iexp]) * density

    else:

        for iexp in range(nCountries):
            ExpectFOC[:,iexp] = MPKtom[:,iexp]*lambPr * density


    return ExpectFOC
