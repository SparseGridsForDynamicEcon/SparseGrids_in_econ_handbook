import numpy as np
import copy
from scipy import optimize
import matplotlib.pyplot as plt

from parameters import *
from setup_num_int import *
from setup_asg import *
from Expect_FOC import *
from SOE import *
from EulerError import *

import sys
sys.path.append("/usr/local/share/Tasmanian/python")

import Tasmanian
import time





# Initial Guess

polGuess = np.zeros((aNum,2*nCountries+1))

for i0 in range(nCountries):
    polGuess[:,i0] = aPoints[:,i0]
    polGuess[:,nCountries+1+i0] = -delta*aPoints[:,i0]

for i1 in range(nCountries):
    polGuess[:,nCountries] += np.exp(aPoints[:,nCountries+i1])*A_tfp*aPoints[:,i1]**zeta - delta*aPoints[:,i1]

polGuess[:,nCountries] = (polGuess[:,nCountries]/(nCountries*pareto**(1.0/gamma)))**(-gamma)



# Load function values into grid structure
grid0.loadNeededPoints(polGuess)

# Main Loop

time_start = time.time()

for iter0 in range(maxiter):

    aPoints = grid0.getPoints()
    polGuess1 = copy.copy(polGuess)

    # Set up grid1 -> grid structure that is allowed to change within the iteration
    grid1 = Tasmanian.makeLocalPolynomialGrid(gridDim,gridOut,gridDepth,gridOrder,gridRule)
    grid1.setDomainTransform(gridDomain)

    # Loop over refinement steps

    for iref in range(maxRef):

        if (grid1.getNumNeeded()>0):

            # Get the grid points that need function values
            aPoints1 = grid1.getNeededPoints()
            # Get the number of grid points that need function values
            aNumAdd = grid1.getNumNeeded()

            # Preallocation of policy update
            polUpdate = np.empty((aNumAdd,2*nCountries+1))

            # Array for intermediate update step
            polInt = np.empty((aNumAdd,2*nCountries+1))

            # Time Iteration Step
            for ii1 in range(aNumAdd):

                state = aPoints1[ii1]
                pol = polGuess1[ii1,:]
                root = optimize.fsolve(sysOfEqs, pol, args=(state,grid0))
                polInt[ii1,:] = root

            # Add the new function values to grid1
            grid1.loadNeededPoints(polInt)

            # We start the refinement process after a given number of iterations
            if (iter0>iterRefStart):
                aNumLoad = grid1.getNumLoaded()
                # Scaling to only allow for those policies that are supposed to
                # determine the refinement process
                scaleCorrMat = np.zeros((aNumLoad,2*nCountries+1))
                scaleCorrMat[:,0:2*nCountries+2] = scaleCorr

                grid1.setSurplusRefinement(tol, dimRef, "classic", [], scaleCorrMat)

                if (grid1.getNumNeeded()>0):

                    # Get vector with coordinates
                    nwpts = grid1.getNeededPoints()
                    aNumNew = grid1.getNumNeeded()

                    # We assign (for now) function values through interpolation
                    polGuess1 = np.zeros((aNumNew,2*nCountries+1))
                    polGuess1 = grid1.evaluateBatch(nwpts)

            else:
                break

        else:
            break


    # Evaluate errors on tomorrow's policy grid
    aPointsEval = grid1.getPoints()
    aNumTot = grid1.getNumPoints()

    polGuessTr1 = grid1.evaluateBatch(aPointsEval)
    polGuessTr0 = grid0.evaluateBatch(aPointsEval)

    # 1) Sup-Norm

    metricAux = np.empty(2*nCountries+1)

    for imet in range(2*nCountries+1):
        metricAux[imet] = np.amax(np.abs(polGuessTr0[:,imet]-polGuessTr1[:,imet]))


    metricSup = np.amax(metricAux)

    # 2) L2-Norm

    metricL2 = 0.0

    for imetL2 in range(2*nCountries+1):
        metricL2 += np.sum((np.abs(polGuessTr0[:,imetL2]-polGuessTr1[:,imetL2]))**2)

    metricL2 = (metricL2/(aNumTot*(2*nCountries+1)))**0.5

    # Now update polGuess and grid0

    polGuess = np.zeros((aNumTot,2*nCountries+1))

    for iupd in range(2*nCountries+1):
        polGuess[:,iupd] = 0.5*polGuessTr0[:,iupd] + 0.5*polGuessTr1[:,iupd]

    grid0 = Tasmanian.copyGrid(grid1)

    print(iter0,aNumTot,metricL2,metricSup)

    if (metricL2<tol_ti):
        break



time_end = time.time()

print(time_end-time_start)


##### EULER ERRORS (SIMULATION) ####

np.random.seed(seed=24022021)

TT = 11000

shocks = np.random.normal(0,1,(TT,nCountries+1))

tfpSim = np.empty((TT,nCountries))
capSim = np.empty((TT,nCountries))
lMultSim = np.empty(TT)
gzAlphaSim = np.empty((TT,nCountries))

errorEESim = np.empty((TT,nCountries))
errorICSim = np.empty((TT,nCountries))

errorTotSim = np.empty((TT,nCountries))

tfpSim[0,:] = 0.0
capSim[0,:] = 1.0


for ia in range(TT-1):

    for icou in range(nCountries):
        tfpSim[ia+1,icou] = rhoZ*tfpSim[ia,icou] + sigE*(shocks[ia+1,icou]+shocks[ia+1,nCountries])

    stateSim = np.empty(nCountries*2)
    stateSim[0:nCountries] = capSim[ia]
    stateSim[nCountries:] = tfpSim[ia+1]

    for icou in range(nCountries):
        capSim[ia+1,icou] = grid1.evaluate(stateSim)[icou]
        gzAlphaSim[ia+1,icou] = grid1.evaluate(stateSim)[nCountries+1+icou]

    lMultSim[ia+1] = grid1.evaluate(stateSim)[nCountries]

    policiesSim = np.empty(nCountries*2+1)
    policiesSim[0:nCountries] = capSim[ia+1,:]
    policiesSim[nCountries] = lMultSim[ia+1]
    policiesSim[nCountries+1:] = gzAlphaSim[ia+1,:]

    for icou in range(nCountries):
        errorICSim[ia+1,icou] = 1.0 - capSim[ia+1,icou]/(capSim[ia,icou]*(1.0-delta))

    errorEESim[ia+1,:] = Euler_error(policiesSim,stateSim,grid1)

    for icou in range(nCountries):
        errorTotSim[ia+1,icou] = max(errorEESim[ia+1,icou],errorICSim[ia+1,icou],\
                                        np.minimum(-errorEESim[ia+1,icou],-errorICSim[ia+1,icou]))


errorTotSim = errorTotSim[1000:,:]
errorTotSimDist = np.sort(errorTotSim,axis=None)

print('')
print('Errors in Simulation:')
print('Max error: ', np.log10(np.amax(np.abs(errorTotSimDist[:-100]))))
print('Avg error: ', np.log10(np.mean(np.abs(errorTotSimDist[:]))))


##### EULER ERRORS (STATE SPACE) ####

np.random.seed(seed=21042021)

TT = 10000

tfpSS = np.zeros((TT,nCountries))
capSS = np.zeros((TT,nCountries))

kpSS = np.zeros((TT,nCountries))
lMultSS = np.zeros(TT)
gzAlphaSS = np.zeros((TT,nCountries))

errorEESS = np.zeros((TT,nCountries))
errorICSS = np.zeros((TT,nCountries))
errorTotSS = np.zeros((TT,nCountries))

for i1 in range(TT):

    for icou in range(nCountries):
        tfpSS[i1,icou] = aMin + (aMax-aMin)*np.random.rand(1)

    for icou in range(nCountries):
        capSS[i1,icou] = kMin + (kMax-kMin)*np.random.rand(1)

    stateSS = np.empty(nCountries*2)
    stateSS[0:nCountries] = capSS[i1,:]
    stateSS[nCountries:] = tfpSS[i1,:]

    for icou in range(nCountries):
        kpSS[i1,icou] = grid1.evaluate(stateSS)[icou]
        gzAlphaSS[i1,icou] = grid1.evaluate(stateSS)[nCountries+1+icou]

    lMultSS[i1] = grid1.evaluate(stateSS)[nCountries]

    policiesSS = np.empty(nCountries*2+1)
    policiesSS[0:nCountries] = kpSS[i1,:]
    policiesSS[nCountries] = lMultSS[i1]
    policiesSS[nCountries+1:] = gzAlphaSS[i1,:]

    errorEESS[i1,:] = Euler_error(policiesSS,stateSS,grid1)

    for icou in range(nCountries):
        errorICSS[i1,icou] = 1.0 - kpSS[i1,icou]/(capSS[i1,icou]*(1.0-delta))

    for icou in range(nCountries):
        errorTotSS[i1,icou] = max(errorEESS[i1,icou],errorICSS[i1,icou],\
                                        np.minimum(-errorEESS[i1,icou],-errorICSS[i1,icou]))


errorTotSSDist = np.sort(errorTotSS,axis=None)

print('')
print('Errors in State Space:')
print('Max error: ', np.log10(np.amax(np.abs(errorTotSSDist[:-100]))))
print('Avg error: ', np.log10(np.mean(np.abs(errorTotSSDist[:]))))
