import numpy as np

## Choose parameter values

gamma = 1.0
betta = 0.99
zeta = 0.36
delta = 0.01
rhoZ = 0.95
sigE = 0.01
kappa = 0.5

## Compute the steady state of the endogenous variables

k_ss = 1.0
A_tfp = (1.0-betta*(1.0-delta))/(zeta*betta)
pareto = A_tfp**(1.0/gamma)

## Technical parameters

maxiter = 100
dimRef = -1
tol = 0.001

iterRefStart = 25
tol_ti = 1e-4

nCountries = 2
nShocks = nCountries+1
kMin = 0.8
kMax = 1.2
aMin = -0.8*sigE/(1.0-rhoZ)
aMax = 0.8*sigE/(1.0-rhoZ)

# Adaptive dimensions
maxRef = 5
scaleCorr = np.zeros(nCountries*2+1)
scaleCorr[0:nCountries] = 1
