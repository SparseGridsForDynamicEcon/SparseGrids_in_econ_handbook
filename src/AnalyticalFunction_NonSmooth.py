import numpy as np
from scipy import interpolate,optimize
import matplotlib.pyplot as plt
import copy

import sys
sys.path.append("/usr/local/share/Tasmanian/python")

import Tasmanian

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


## LOCAL SPARSE Grid

gridDim = 2
gridOut = 1
gridDepth = [5,6,7,8,9,10,11,12,13,14,15]
gridOrder = 1
gridRule = "localp"

gridDomain = np.zeros((gridDim,2))
gridDomain[:,1] = 1.0

np.random.seed(seed=28092021)
TT = 10000

evalPoints = np.random.rand(TT,gridDim)

noDepth = len(gridDepth)

## LOCAL SPARSE Grid - Non-Smooth

l1ErrorNSG = np.zeros(noDepth)
l2ErrorNSG = np.zeros(noDepth)
lInfErrorNSG = np.zeros(noDepth)
NSGPoints = np.zeros(noDepth,dtype=int)

# Index for grid_depth
i0 = 0

for ig in gridDepth:

    gridNSG = Tasmanian.makeLocalPolynomialGrid(gridDim,gridOut,ig,gridOrder,gridRule)
    gridNSG.setDomainTransform(gridDomain)
    aPointsNSG = gridNSG.getNeededPoints()
    aNumNSG = gridNSG.getNumPoints()
    NSGPoints[i0] = aNumNSG

    fValuesNSG = np.empty((aNumNSG,gridOut))
    fTrueNSG = np.empty((TT,gridOut))
    fImpNSG = np.empty((TT,gridOut))

    for i1 in range(aNumNSG):

        fValuesNSG[i1] = np.maximum(0.0,1.0 - np.exp(0.5-np.prod(aPointsNSG[i1,:]+0.2)**(1.0/gridDim)))

    gridNSG.loadNeededPoints(fValuesNSG)

    for i2 in range(TT):

        fTrueNSG[i2] = np.maximum(0.0,1.0 - np.exp(0.5-np.prod(evalPoints[i2,:]+0.2)**(1.0/gridDim)))
        fImpNSG[i2] = gridNSG.evaluate(evalPoints[i2,:])[0]

    l1ErrorNSG[i0] = 1.0/TT*np.sum(np.abs(fTrueNSG-fImpNSG))
    l2ErrorNSG[i0] = (1.0/TT*np.sum(np.abs(fTrueNSG-fImpNSG)**2))**0.5
    lInfErrorOrder = np.sort(np.abs(fTrueNSG-fImpNSG),axis=None)
    lInfErrorNSG[i0] = np.amax(lInfErrorOrder[:-10])

    i0 += 1




################################################################################
#
#                       ADAPTIVE SPARSE GRID
#
################################################################################

# Adaptive Sparse Gridpoints

gridDim = 2
gridOut = 1
gridDepthASG = 2
gridOrder = 1
gridRule = "localp"

surplThreshold = 1e-6

l1ErrorNASG = np.zeros(noDepth)
l2ErrorNASG = np.zeros(noDepth)
lInfErrorNASG = np.zeros(noDepth)
NASGPoints = np.zeros(noDepth,dtype=int)


for i0 in range(noDepth):

    gridNASG = Tasmanian.makeLocalPolynomialGrid(gridDim,gridOut,gridDepthASG,gridOrder,gridRule)
    gridNASG.setDomainTransform(gridDomain)
    aPointsNASG = gridNASG.getNeededPoints()
    aNumNASG = gridNASG.getNumPoints()

    aNumNASG1 = copy.copy(aNumNASG)

    ilev = gridDepthASG-1

    while ((gridNASG.getNumNeeded() > 0) and (ilev<gridDepth[i0])):

        aPointsNASG1 = gridNASG.getNeededPoints()
        aNumAdd = gridNASG.getNumNeeded()

        fValuesNASG1 = np.zeros((aNumAdd,gridOut))

        for i1 in range(aNumAdd):
            fValuesNASG1[i1] = np.maximum(0.0,1.0 - np.exp(0.5-np.prod(aPointsNASG1[i1,:]+0.2)**(1.0/gridDim)))

        gridNASG.loadNeededPoints(fValuesNASG1)
        aNumLoad = gridNASG.getNumLoaded()

        aNumNASG1 = gridNASG.getNumPoints()
        #print(i0,aNum1)

        gridNASG.setSurplusRefinement(surplThreshold, -1, "classic")
        ilev += 1

    NASGPoints[i0] = aNumNASG1
    fTrueNASG = np.zeros(TT)
    fImpNASG = np.zeros(TT)

    for i1 in range(TT):

        fTrueNASG[i1] = np.maximum(0.0,1.0 - np.exp(0.5-np.prod(evalPoints[i1,:]+0.2)**(1.0/gridDim)))
        fImpNASG[i1] = gridNASG.evaluate(evalPoints[i1,:])[0]

    l1ErrorNASG[i0] = 1.0/TT*np.sum(np.abs(fTrueNASG-fImpNASG))
    l2ErrorNASG[i0] = (1.0/TT*np.sum(np.abs(fTrueNASG-fImpNASG)**2))**0.5
    lInfErrorOrder = np.sort(np.abs(fTrueNASG-fImpNASG),axis=None)
    lInfErrorNASG[i0] = np.amax(lInfErrorOrder[:-10])





fig, ax = plt.subplots(1,2,figsize=(14,8))
ax[0].scatter(NSGPoints,lInfErrorNSG, marker='o',label='$SG$')
ax[0].plot(NSGPoints,lInfErrorNSG)
ax[0].scatter(NASGPoints,lInfErrorNASG, marker = 'x', label='$ASG$')
ax[0].plot(NASGPoints,lInfErrorNASG,'--')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlabel("Points", fontsize=14)
ax[0].set_ylabel("$L_{\infty}$ error", fontsize=14)
ax[0].legend()

ax[1].scatter(NSGPoints,l1ErrorNSG, marker='o',label='$SG$')
ax[1].plot(NSGPoints,l1ErrorNSG)
ax[1].scatter(NASGPoints,l1ErrorNASG, marker = 'x', label='$ASG$')
ax[1].plot(NASGPoints,l1ErrorNASG,'--')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlabel("Points", fontsize=14)
ax[1].set_ylabel("$L_1$ error", fontsize=14)
ax[1].legend()
plt.show()
