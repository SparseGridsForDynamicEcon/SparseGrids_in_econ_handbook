import Tasmanian
import numpy as np
from parameters import *

## Construct grid

gridDim = nCountries*2
gridOut = nCountries*2+1
gridDepth = 2
gridOrder = 1
gridRule = "localp"

gridDomain = np.empty((gridDim,2))

gridDomain[0:nCountries,0] = kMin
gridDomain[0:nCountries,1] = kMax
gridDomain[nCountries:,0] = aMin
gridDomain[nCountries:,1] = aMax

grid0 = Tasmanian.makeLocalPolynomialGrid(gridDim,gridOut,gridDepth,gridOrder,gridRule)
grid0.setDomainTransform(gridDomain)
aPoints = grid0.getPoints()
aNum = grid0.getNumPoints()
