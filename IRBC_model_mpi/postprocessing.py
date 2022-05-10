"""
This library contains routines to postprocess the model solutions.
In particular, it consists of:

	- computing the error measures along the simulation path
	- computing the error measures in the state space
	- computing Euler equation errors
	- plotting and storing policy functions

"""

import numpy as np
from parameters import *
from setup_asg import *
from aux_fcts import *
from setup_num_int import *
from Expect_FOC import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
from mpi import rank

# TeX Support for plots
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

import Tasmanian



################################################################################
#                   Errors along the simulation path                           #
################################################################################

def errors_sim(grid):
    if rank !=0:
        return
    # Set random seed and draw shocks
    np.random.seed(seed=24022021)
    shocks = np.random.normal(0,1,(TT+burnin,nCountries+1))

    # Preallocation of simulation arrays
    tfpSim = np.zeros((TT+burnin,nCountries))
    capSim = np.zeros((TT+burnin,nCountries))
    lMultSim = np.zeros(TT+burnin)

    if typeIRBC=='non-smooth':
        gzAlphaSim = np.zeros((TT+burnin,nCountries))

    errorEESim = np.zeros((TT+burnin,nCountries))
    if typeIRBC=='non-smooth':
        errorICSim = np.zeros((TT+burnin,nCountries))
    errorTotSim = np.zeros((TT+burnin,nCountries))

    # Initialization: Steady state values
    tfpSim[0,:] = 0.0
    capSim[0,:] = 1.0

    # Simulation loop
    for ia in range(TT+burnin-1):

        for icou in range(nCountries):
            tfpSim[ia+1,icou] = rhoZ*tfpSim[ia,icou] + sigE*(shocks[ia+1,icou]+shocks[ia+1,nCountries])

        stateSim = np.zeros(nCountries*2)
        stateSim[0:nCountries] = capSim[ia]
        stateSim[nCountries:] = tfpSim[ia+1]

        for icou in range(nCountries):
            capSim[ia+1,icou] = grid.evaluate(stateSim)[icou]
            if typeIRBC=='non-smooth':
                gzAlphaSim[ia+1,icou] = grid.evaluate(stateSim)[nCountries+1+icou]

        lMultSim[ia+1] = grid.evaluate(stateSim)[nCountries]

        policiesSim = np.zeros(nCountries*2+1)
        policiesSim[0:nCountries] = capSim[ia+1,:]
        policiesSim[nCountries] = lMultSim[ia+1]
        if typeIRBC=='non-smooth':
            policiesSim[nCountries+1:] = gzAlphaSim[ia+1,:]

        if typeIRBC=='non-smooth':
            for icou in range(nCountries):
                errorICSim[ia+1,icou] = 1.0 - capSim[ia+1,icou]/(capSim[ia,icou]*(1.0-delta))

        errorEESim[ia+1,:] = Euler_error(policiesSim,stateSim,grid)

        if typeIRBC=='non-smooth':
            for icou in range(nCountries):
                errorTotSim[ia+1,icou] = max(errorEESim[ia+1,icou],errorICSim[ia+1,icou],\
                                            np.minimum(-errorEESim[ia+1,icou],-errorICSim[ia+1,icou]))
        else:
            errorTotSim[ia+1,:] = errorEESim[ia+1,:]


    errorTotSim = errorTotSim[1000:,:]
    errorTotSimDist = np.sort(errorTotSim,axis=None)

    print('')
    print('Errors in Simulation:')
    print('Max error: ', np.log10(np.amax(np.abs(errorTotSimDist[:-100]))))
    print('Avg error: ', np.log10(np.mean(np.abs(errorTotSimDist[:]))))


    return


################################################################################
#                       Errors in the state space                              #
################################################################################

def errors_ss(grid):
    if rank != 0:
        return
    # Set random seed
    np.random.seed(seed=21042021)

    # Preallocation of simulation arrays
    tfpSS = np.zeros((TT,nCountries))
    capSS = np.zeros((TT,nCountries))

    kpSS = np.zeros((TT,nCountries))
    lMultSS = np.zeros(TT)
    if typeIRBC=='non-smooth':
        gzAlphaSS = np.zeros((TT,nCountries))

    errorEESS = np.zeros((TT,nCountries))
    if typeIRBC=='non-smooth':
	    errorICSS = np.zeros((TT,nCountries))
    errorTotSS = np.zeros((TT,nCountries))

    for i1 in range(TT):

        for icou in range(nCountries):
            tfpSS[i1,icou] = aMin + (aMax-aMin)*np.random.rand(1)

        for icou in range(nCountries):
            capSS[i1,icou] = kMin + (kMax-kMin)*np.random.rand(1)

        stateSS = np.zeros(nCountries*2)
        stateSS[0:nCountries] = capSS[i1,:]
        stateSS[nCountries:] = tfpSS[i1,:]

        for icou in range(nCountries):
            kpSS[i1,icou] = grid.evaluate(stateSS)[icou]
            if typeIRBC=='non-smooth':
                gzAlphaSS[i1,icou] = grid.evaluate(stateSS)[nCountries+1+icou]

        lMultSS[i1] = grid.evaluate(stateSS)[nCountries]

        policiesSS = np.zeros(nPols)
        policiesSS[0:nCountries] = kpSS[i1,:]
        policiesSS[nCountries] = lMultSS[i1]
        if typeIRBC=='non-smooth':
            policiesSS[nCountries+1:] = gzAlphaSS[i1,:]
        errorEESS[i1,:] = Euler_error(policiesSS,stateSS,grid)

        if typeIRBC=='non-smooth':
            for icou in range(nCountries):
                errorICSS[i1,icou] = 1.0 - kpSS[i1,icou]/(capSS[i1,icou]*(1.0-delta))
            for icou in range(nCountries):
                errorTotSS[i1,icou] = max(errorEESS[i1,icou],errorICSS[i1,icou],\
                                    np.minimum(-errorEESS[i1,icou],-errorICSS[i1,icou]))
        else:
            errorTotSS[i1,:] = errorEESS[i1,:]


    errorTotSSDist = np.sort(errorTotSS,axis=None)

    print('')
    print('Errors in State Space:')
    print('Max error: ', np.log10(np.amax(np.abs(errorTotSSDist[:-100]))))
    print('Avg error: ', np.log10(np.mean(np.abs(errorTotSSDist[:]))))

    return



################################################################################
#                       Euler equation errors                                  #
################################################################################

def Euler_error(x,state,pols):

    capStates = state[0:nCountries]
    tfpStates = state[nCountries:]

    capPolicies = x[0:nCountries]
    lamb = x[nCountries]

    if typeIRBC=='non-smooth':

        gzAlphas = x[nCountries+1:]
        gzAlphaPlus = np.maximum(0.0,gzAlphas)
        gzAlphaMinus = np.maximum(0.0,-gzAlphas)

    Integrands = np.empty((numNodes,nCountries))
    newstate = np.empty(nCountries)
    captomtom = np.empty(nCountries)
    MPKtom = np.empty(nCountries)

    if typeIRBC=='non-smooth':

        gzAlpTomPl = np.empty(nCountries)

    #Compute Density
    if typeInt=='GH-quadrature':
        density = np.pi**(-(nCountries+1) * 0.5)
    else:
        density = 1.0

    for i_int in range(numNodes):

        for icou in range(nCountries):
            newstate[icou] = rhoZ*tfpStates[icou] + (IntNodes[i_int,icou] + IntNodes[i_int,nCountries])

        state_Sim = np.empty(nCountries*2)
        state_Sim[0:nCountries] = capPolicies
        state_Sim[nCountries:] = newstate

        lambtom = pols.evaluate(state_Sim)[nCountries]

        for icou in range(nCountries):
            captomtom[icou] = pols.evaluate(state_Sim)[icou]
            if typeIRBC=='non-smooth':
                gzAlpTomPl[icou] = np.maximum(0.0,pols.evaluate(state_Sim)[nCountries+1+icou])
            MPKtom[icou] = 1.0 - delta + Fk(capPolicies[icou],newstate[icou]) - AdjCost_k(capPolicies[icou],captomtom[icou])
            if typeIRBC=='non-smooth':
                Integrands[i_int,icou] = (lambtom*MPKtom[icou] - (1.0-delta)*gzAlpTomPl[icou]) * density
            else:
                Integrands[i_int,icou] = lambtom*MPKtom[icou] * density


    IntResult = np.empty(nCountries)

    for iint in range(nCountries):
        IntResult[iint] = np.dot(IntWeights,Integrands[:,iint])

    res = np.zeros(nCountries)

    for ires in range(nCountries):
        res[ires] = betta*IntResult[ires]/(lamb*(1.0 + AdjCost_ktom(capStates[ires],capPolicies[ires]))) - 1.0


    return res



################################################################################
#                   Plotting and storing policy functions                      #
################################################################################


def plot_policies(plotDim,grid):
    if rank != 0:
        return
    # Number of points in plot
    plotPoints = 50

    plotArray = np.linspace(gridDomain[plotDim,0],gridDomain[plotDim,1],plotPoints)


    # We fix every dimension at its steady state except plotDim
    steadyState = np.zeros(nCountries*2)
    steadyState[0:nCountries] = k_ss
    evalArray = steadyState*np.ones((plotPoints,nCountries*2))
    evalArray[:,plotDim] = plotArray

    policyArray = np.zeros((plotPoints,nCountries*2+1))
    # The first nCountries columns are grid points
    policyArray[:,:-1] = evalArray

    for iPlot in range(nPols):

        # The last column is the respective policy
        policyArray[:,-1] = grid.evaluateBatch(evalArray)[:,iPlot]

        if typeIRBC=='non-smooth':
            # Capital Policies
            if (iPlot<nCountries):
                np.savetxt(data_location + "Cap_policy_" + str(iPlot) + ".txt", policyArray)

                fig,ax = plt.subplots(figsize=(16,8))
                ax.plot(plotArray,policyArray[:,-1])
                plt.savefig(data_location + "Cap_policy_" + str(iPlot) + ".png", bbox_inches='tight')
                plt.close()

            # ARC multiplier
            elif (iPlot==nCountries):
                np.savetxt(data_location + "ARC_policy.txt", policyArray)

                fig,ax = plt.subplots(figsize=(16,8))
                ax.plot(plotArray,policyArray[:,-1])
                plt.savefig(data_location + "ARC_policy.png", bbox_inches='tight')
                plt.close()

            # Investement constraint multipliers
            else:
                np.savetxt(data_location + "IC_policy_" + str(iPlot-nCountries-1) + ".txt", policyArray)

                fig,ax = plt.subplots(figsize=(16,8))
                ax.plot(plotArray,policyArray[:,-1])
                plt.savefig(data_location + "IC_policy_" + str(iPlot-nCountries-1) + ".png", bbox_inches='tight')
                plt.close()

        else:
            # Capital Policies
            if (iPlot<nCountries):
                np.savetxt(data_location + "Cap_policy_" + str(iPlot) + ".txt", policyArray)

                fig,ax = plt.subplots(figsize=(16,8))
                ax.plot(plotArray,policyArray[:,-1])
                plt.savefig(data_location + "Cap_policy_" + str(iPlot) + ".png", bbox_inches='tight')
                plt.close()

            # ARC multiplier
            else:
                np.savetxt(data_location + "ARC_policy.txt", policyArray)

                fig,ax = plt.subplots(figsize=(16,8))
                ax.plot(plotArray,policyArray[:,-1])
                plt.savefig(data_location + "ARC_policy.png", bbox_inches='tight')
                plt.close()
