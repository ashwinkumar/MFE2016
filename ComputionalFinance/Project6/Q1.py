from ExoticOption import LookBackOption as lbo
import numpy as np
import matplotlib.pyplot as plt

#### Question 1 - Fixed Strike LookBack Call and Put

S0 = 98
K = 100
r = 0.03
T = 1
ndiv = 100
dT = T/ndiv
nsims = 10000
start_vol = 0.12
end_vol = 0.48
step_vol = 0.04
sigma = np.arange(start_vol,end_vol+step_vol,step_vol)
l_it = len(sigma)
lbfsCall = np.zeros(l_it)
lbfsPut = np.zeros(l_it)
for i in range(0,l_it):
    lbo_class = lbo(S0,K,r,sigma[i],T,dT,nsims)
    lbo_class.mcStockTree()
    lbfsCall[i],_ = lbo_class.getFixedStrikeCall()
    lbfsPut[i],_ = lbo_class.getFixedStrikePut()
    del lbo_class

plt.plot(sigma, lbfsCall)
plt.ylabel("Option Price")
plt.xlabel("Vol")
plt.title("Q1 a) Fixed Strike Lookback Call MC Simulation")
plt.savefig("Proj6_1a.png")
plt.clf()
plt.plot(sigma, lbfsPut)
plt.ylabel("Option Price")
plt.xlabel("Vol")
plt.title("Q1 a) Fixed Strike Lookback Put MC Simulation")
plt.savefig("Proj6_1b.png")
