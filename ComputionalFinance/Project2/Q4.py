import mynormal as nl
import MonteCarlo as mc
import math
import numpy as np
import BlackScholesOptionPrice as bs
## Question 4

'''
a) Estimate the price c of a European Call option on the stock following geometric brownian motion
'''

nsim = 10000
sigma = 0.2
r = 0.04
S0 = 88
T = 5
X = 100


nldist = nl.NormalDistribution("Polar-Marsaglia")
nldist.generateNormalDistribution(nsim)
z = nldist.getRdmNumbers()

mc_call = mc.mc_callprice(S0, X, r, sigma, T,  z, False)
print('Q4. a) European Call Option Price by MC Simulation = %f with std dev = %f' % (mc_call[0], mc_call[1]))

'''
Compute the exact value of the option c by the Black-Scholes formula. Now use variance reduction techniques (whichever you want)
to estimate the price in part (a) again. Did the accuracy improve? Comment.
'''
bsCallPrice = bs.EuropeanCallPrice(S0, X, r, sigma, T)
print('Q4. b) Black Scholes Price = %f' % bsCallPrice)

# We generate only half random numbers. nsim = nsim/2
nldist = nl.NormalDistribution("Polar-Marsaglia")
nldist.generateNormalDistribution(nsim/2)
z = nldist.getRdmNumbers()
mc_call_var = mc.mc_callprice(S0, X, r, sigma, T,  z, True)
print('Q4. b) European Call Option Price by MC Simulation (reduced var) = %f with std dev = %f' % (mc_call_var[0], mc_call_var[1]))




