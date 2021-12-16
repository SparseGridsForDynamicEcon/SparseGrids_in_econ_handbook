"""
This file contains auxilliary functions that are used within the code.

"""


from parameters import *
import numpy as np


################################################################################
#                           Production function                                #
################################################################################

def F(capital,sh):

    FF = A_tfp * np.exp(sh)*np.maximum(capital,1e-6)**zeta

    return FF


################################################################################
#                      Marginal product of capital                             #
################################################################################

def Fk(capital,sh):

    F_k =  A_tfp * zeta*np.exp(sh)*np.maximum(capital,1e-6)**(zeta-1.0)

    return F_k


################################################################################
#                          Capital adjustment cost                             #
################################################################################

def AdjCost(ktod,ktom):

    captod = np.maximum(ktod,1e-6)
    captom = np.maximum(ktom,1e-6)

    j = captom/captod - 1.0
    Adj_cost = 0.5 * kappa * j * j * captod

    return Adj_cost


################################################################################
#        Derivative of capital adjustment cost w.r.t today's cap stock         #
################################################################################

def AdjCost_k(ktod,ktom):

    captod = np.maximum(ktod,1e-6)
    captom = np.maximum(ktom,1e-6)

    j = captom/captod - 1.0
    j1 = captom/captod + 1.0
    AdjCostk = (-0.5)*kappa*j*j1

    return AdjCostk


################################################################################
#      Derivative of capital adjustment cost w.r.t tomorrows's cap stock       #
################################################################################

def AdjCost_ktom(ktod,ktom):

    captod = np.maximum(ktod,1e-6)
    captom = np.maximum(ktom,1e-6)

    j = captom/captod - 1.0
    AdjCostktom = kappa * j


    return AdjCostktom
