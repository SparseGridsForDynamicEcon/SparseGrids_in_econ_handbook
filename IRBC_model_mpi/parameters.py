"""
This file sets the economic parameters as well as the technical parameters that
govern the time iteration and the error computation process.

"""

import numpy as np

## Economic parameters


# Type of IRBC model ({'smooth','non-smooth'})
typeIRBC = 'non-smooth'
# Type of numerical integration: ({'GH-quadrature','monomials_2d','monomials_power'})
typeInt = 'monomials_power'

# Number of countries
nCountries = 3
# Number of shocks (Country-specific shocks + aggregate shock)
nShocks = nCountries+1
# Number of policies (nCountries+1 for smooth IRBC, nCountries*2+1 for nonsmooth)
if typeIRBC=='non-smooth':
    nPols = nCountries*2+1
elif typeIRBC=='smooth':
    nPols = nCountries+1
else:
    print('Error in determining the type of IRBC')
    exit()

# Intertemporal elasticity of substitution
a_eis = 0.25
b_eis = 1.0
gamma = np.zeros(nCountries)
for j1 in range(nCountries):
    gamma[j1] = a_eis+j1*(b_eis-a_eis)/(nCountries-1)
# Discount factor
betta = 0.99
# Capital share of income
zeta = 0.36
# Depreciation rate
delta = 0.01
# Persistence of TFP shocks
rhoZ = 0.95
# Standard deviation of TFP shocks
sigE = 0.01
# Intensity of capital adjustment costs
kappa = 0.5
# Steady state for capital
k_ss = 1.0
# Aggregate productivity
A_tfp = (1.0-betta*(1.0-delta))/(zeta*betta)
# Welfare weight
pareto = A_tfp**(1.0/gamma)


# Lower bound for capital
kMin = 0.8
# Upper bound for capital
kMax = 1.2
# Lower bound for TFP
aMin = -0.8*sigE/(1.0-rhoZ)
# Upper bound for TFP
aMax = 0.8*sigE/(1.0-rhoZ)


## Technical parameters

# Iteration to start at (start from scratch -- numstart =0; restart: numstart >0)
numstart = 0
# Maximum number of iterations
maxiter = 30
# Iteration at which the refinement starts
iterRefStart = 25
# Convergence criterion in time iteration
tol_ti = 1e-4
# Number of random draws for the error computation
TT = 10000
# Number of burn-in periods for SIMULATION
burnin = 1000

# Location where data is stored
if typeIRBC == "non-smooth":
    data_location = "data/data_nonsmooth/"
else:
    data_location = "data/data_smooth/"
# Frequency of saving grid
savefreq = 100
