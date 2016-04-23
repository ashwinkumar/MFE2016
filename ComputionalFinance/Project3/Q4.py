import math
import EulerDiscretization as ed
import numpy as np
import timeit
## Question 4

'''
Implementation of Heston model
'''
rho = -0.6
r = 0.03
S0 = 48
K = 50
v0 = 0.05
sigma = 0.42
alpha = 5.8
beta = 0.0625
T = 0.5
nsims = 1000
ndiv = 100


# We will initialize the Euler Discretization (vectorized) with the initial values of the series we are trying to simulate. T
ED0 = np.array([S0,v0])
const_dt = np.array([alpha,beta])
const_dw = np.array([sigma,rho])

def coeff_EDdw_FT(ed,y,t,z,cdw):
    s = ed[0]
    v = ed[1]
    sigma = cdw[0]
    rho = cdw[1]
    z1 = z[0]
    z2 = z[0]*rho + math.sqrt(1 - rho*rho)*z[1]
    x1 = s * math.sqrt(max(v,0))*z1

    x2 = sigma*math.sqrt(max(v,0))*z2

    res = np.array([x1, x2])
    return res

def coeff_EDdt_FT(ed,y,t,cdt):
    x1 = ed[0]*y
    alpha = cdt[0]
    beta = cdt[1]
    x2= alpha*( beta- max(ed[1],0))
    res = np.array([x1,x2])
    return res

def coeff_EDdw_PT(ed,y,t,z,cdw):
    s = ed[0]
    v = ed[1]
    sigma = cdw[0]
    rho = cdw[1]
    z1 = z[0]
    z2 = z[0]*rho + math.sqrt(1 - rho*rho)*z[1]
    x1 = s * math.sqrt(max(v,0))*z1

    x2 = sigma*math.sqrt(max(v,0))*z2

    res = np.array([x1, x2])
    return res

def coeff_EDdt_PT(ed,y,t,cdt):
    x1 = ed[0]*y
    alpha = cdt[0]
    beta = cdt[1]
    x2= alpha*( beta- ed[1])
    res = np.array([x1,x2])
    return res

def coeff_EDdw_RM(ed,y,t,z,cdw):
    s = ed[0]
    v = ed[1]
    sigma = cdw[0]
    rho = cdw[1]
    z1 = z[0]
    z2 = z[0]*rho + math.sqrt(1 - rho*rho)*z[1]
    x1 = s * math.sqrt(abs(v))*z1

    x2 = sigma*math.sqrt(abs(v))*z2

    res = np.array([x1, x2])
    return res

def coeff_EDdt_RM(ed,y,t,cdt):
    x1 = ed[0]*y
    alpha = cdt[0]
    beta = cdt[1]
    x2= alpha*( beta- abs(ed[1]))
    res = np.array([x1,x2])
    return res

def ECall(S,X,r,t):
    array_zeros = np.zeros(len(S))
    payoff = np.maximum(S-X,array_zeros)
    disc = math.exp(-r*t)
    return disc*np.mean(payoff), disc*disc*np.var(payoff, ddof=1)/len(S)


# ED will be 3-D
ED_FT = np.zeros([nsims, 2,ndiv+1])
ED_PT = np.zeros([nsims, 2,ndiv+1])
ED_RM = np.zeros([nsims, 2,ndiv+1])


start_time1 = timeit.default_timer()
for i in range(0,nsims):
    ED_FT[i,:] = ed.EulerDiscretization(ED0, const_dt, const_dw, coeff_EDdt_FT, r, coeff_EDdw_FT, [], T, ndiv, 2)
    ED_PT[i,:] = ed.EulerDiscretization(ED0, const_dt,const_dw,coeff_EDdt_PT, r ,coeff_EDdw_PT ,[],T, ndiv,2)
    ED_RM[i,:] = ed.EulerDiscretization(ED0, const_dt, const_dw, coeff_EDdt_RM, r, coeff_EDdw_RM, [], T, ndiv, 2)


callprice_FT, var_callprice_FT = ECall(ED_FT[:,0,ndiv],K,r,T) # Full Truncation
callprice_PT, var_callprice_PT = ECall(ED_PT[:,0,ndiv],K,r,T) # Partial Truncation
callprice_RM, var_callprice_RM = ECall(ED_RM[:,0,ndiv],K,r,T) # Reflection Methods




elapsed1 = timeit.default_timer() - start_time1


print('Q4. European Call Option Prices and respective variances of MC Simulation  using stochastic volatalities :  ****[%f sec]' %elapsed1)
print('    a) By Full Truncation = %f with var of MC = %f' %(callprice_FT,var_callprice_FT))
print('    b) By Partial Truncation = %f with var of MC = %f' %(callprice_FT,var_callprice_FT))
print('    c) By Reflection Method = %f with var of MC = %f' %(callprice_FT,var_callprice_FT))
