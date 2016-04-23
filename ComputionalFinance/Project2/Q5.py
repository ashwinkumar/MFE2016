import math
import MonteCarlo as mc
import numpy as np
import matplotlib.pyplot as plt
import random as rdm

# Queston 5
'''
(a) For each integer number 𝑛 from 1 to 10, use 1000 simulations of 𝑆𝑛 to estimate 𝐸𝑆𝑛, where 𝑆𝑡 is a Geometric Brownian Motion process:
Plot all of the above 𝐸(𝑆𝑛), for 𝑛 ranging from 1 to 10, in one graph.
'''


S0 = 88
sigma = 0.18
r = 0.04
nsims = 1000

T = 10
st= list(mc.mc_stockprice(S0, r, sigma, T, nsims))
time_linspace = np.linspace(1, len(st), len(st))
plt.plot(time_linspace, st)
plt.ylabel("E(t) for Brownian Stock Price")
plt.xlabel("Time in years")
plt.title("Q5 a) E(St)")
plt.show()
print('Q5. a) Last value of E(St) = %f' % st[T-1])
'''
(b) Now simulate 6 paths of 𝑆𝑡 for 0≤𝑡≤10 (defined in part (a)) by dividing up the interval [0, 10] into 1,000 equal parts.
'''
T = 10
ndiv = 1000
paths = 6
time_linspace1 = np.linspace(1, T, ndiv)
for i in range(paths):
    ST = list(mc.mc_stockpath(S0, r, sigma, T, ndiv, rdm.randint(10,50)))
    plt.plot(time_linspace1, ST)

plt.plot(time_linspace, st)
plt.ylabel("St for Brownian Stock Price")
plt.xlabel("Time in years")
plt.title("Q5 c) Stock Price Simulation")
plt.show()

'''
(c) What would happen to the 𝐸𝑆𝑛 graph if you increased 𝜎 from 18% to 35%?
What would happen to the 6 plots of 𝑆𝑡 for 0≤𝑡≤10, if you increased 𝜎 from 18% to 35%?
'''



sigma = 0.35
for i in range(paths):
    ST = list(mc.mc_stockpath(S0, r, sigma, T, ndiv, rdm.randint(10,50)))
    plt.plot(time_linspace1, ST)

st = list(mc.mc_stockprice(S0, r, sigma, T, nsims))

plt.plot(time_linspace, st)
plt.ylabel("St for Brownian Stock Price")
plt.xlabel("Time in years")
plt.title("Q5 c) Stock Price Simulation")
plt.show()
print('Q5. d) Last value of E(St) = %f' %st[T-1])

