import math
import numpy as np
import mynormal as nl

def myWienerMC1(z,SQRT_T):
    done = object()
    n1 = next(z, done)
    while (n1 is not done) :
        wt = n1 * SQRT_T
        tmp = wt*wt + math.sin(wt)
        yield tmp
        n1 = next(z, done)



def myWienerMC1_vr(z, SQRT_T):
    done = object()
    n1 = next(z, done)
    while (n1 is not done) :
        wt = n1 * SQRT_T
        tmp1 = wt*wt + math.sin(wt)
        yield tmp1
        tmp2 = wt*wt + math.sin(-wt)
        yield tmp2
        n1 = next(z, done)



def myCos(z, sqrt_t):
    done = object()
    n1 = next(z, done)
    while (n1 is not done):
        yield np.cos(n1 * sqrt_t)
        n1 = next(z, done)




def myCos_vr(z, gamma, sqrt_t):
    done = object()
    n1 = next(z, done)
    while (n1 is not done):
        tmp1 = np.cos(n1 * sqrt_t)
        tmp2 = np.multiply(n1*n1,sqrt_t) - sqrt_t
        yield tmp1 -gamma*tmp2
        n1 = next(z, done)


def mc_norm_payoff(const_drift, const_sigt, const_x_s, z):

    done = object()
    n1 = next(z, done)
    num = 0
    while (n1 is not done):
        #tmp = math.exp(const_drift + const_sigt * n1 - const_x_s)
        yield max(0, math.exp(const_drift + const_sigt * n1) - const_x_s)
        n1 = next(z, done)
        num += 1

def mc_norm_payoff_variancered(const_drift, const_sigt, const_x_s, z):

    done = object()
    n1 = next(z, done)
    num = 0
    while (n1 is not done):
        yield max(0, math.exp(const_drift + const_sigt * n1) - const_x_s)
        yield max(0, math.exp(const_drift + const_sigt * -n1) - const_x_s)
        n1 = next(z, done)
        num += 1

def mc_callprice(S0, X, r, sigma, T, z, bool = False):
    const_x_s = X / S0
    disc = math.exp(-r*T)
    const_drift = (r - sigma * sigma / 2.0) * T
    const_sigt = sigma * math.sqrt(T)
    if bool == False:
        l_payoff = list(mc_norm_payoff(const_drift, const_sigt, const_x_s, z))
    else:
        l_payoff = list(mc_norm_payoff_variancered(const_drift, const_sigt, const_x_s, z))
    return disc*S0*(np.mean(l_payoff)) , (disc*S0*np.std(l_payoff, ddof=1)/ math.sqrt(len(l_payoff)))


def mc_stockprice(S0, r, sigma, T, nsims):
    normal = nl.NormalDistribution()
    for i in range(1,T+1):
        normal.generateNormalDistribution(nsims)
        z = normal.getRdmNumbers()
        const_drift = (r - sigma * sigma / 2.0) * i
        const_sigt = sigma * math.sqrt(i)
        l_payoff = list(mc_norm_payoff_variancered(const_drift, const_sigt, 0, z))
        yield S0 *(np.mean(l_payoff))

def mc_calc_exp(const_drift, const_sigt, z, S0):
    done = object()
    n1 = next(z, done)
    while (n1 is not done):
        S = S0*(math.exp(const_drift + const_sigt* n1))
        S0 = S
        yield S
        n1 = next(z, done)


def mc_stockpath(S0, r, sigma, T, ndiv, seed=14):
    ndiv = (int)(ndiv)
    normal = nl.NormalDistribution("Polar-Marsaglia",seed)
    normal.generateNormalDistribution(ndiv)

    dt = T /ndiv
    const_drift = (r - sigma * sigma / 2.0) * dt
    const_sigt = sigma * math.sqrt(dt)
    z = normal.getRdmNumbers()
    return mc_calc_exp(const_drift, const_sigt,z, S0)


