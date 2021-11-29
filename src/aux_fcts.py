from parameters import *
import numpy as np


# Auxilliary functions

def F(capital,sh):

    FF = A_tfp * np.exp(sh)*np.maximum(capital,1e-6)**zeta

    return FF


#--------------------------------------------------------------------

def Fk(capital,sh):

    F_k =  A_tfp * zeta*np.exp(sh)*np.maximum(capital,1e-6)**(zeta-1.0)

    return F_k


#--------------------------------------------------------------------

def AdjCost(ktod,ktom):

    captod = ktod
    captom = ktom

    if (captod < 1e-6):
        print("Error - Today's Capital Stock is negative in AdjCost")
        captod = 1e-6
    elif (captom < 1e-6):
        print("Error - Tomorrow's Capital Stock is negative in AdjCost")
        captom = 1e-6
    else:
        j = captom/captod - 1.0
        Adj_cost = 0.5 * kappa * j * j * captod

    return Adj_cost

#----------------------------------------------------------------------

def AdjCost_k(ktod,ktom):

    captod = np.maximum(ktod,1e-6)
    captom = np.maximum(ktom,1e-6)

    j = captom/captod - 1.0
    j1 = captom/captod + 1.0
    AdjCostk = (-0.5)*kappa*j*j1

    return AdjCostk


#-------------------------------------------------------------------------

def AdjCost_ktom(ktod,ktom):

    captod = ktod
    captom = ktom

    if (captod < 1e-6):
        print("Error - Today's Capital Stock is negative in AdjCost_ktom")
        captod = 1e-6
    elif (captom < 1e-6):
        print("Error - Tomorrow's Capital Stock is negative in AdjCost_ktom")
        captom = 1e-6
    else:
        j = captom/captod - 1.0
        AdjCostktom = kappa * j
    return AdjCostktom
