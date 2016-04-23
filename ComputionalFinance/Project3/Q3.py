import math
import mynormal as nl
import BlackScholesOptionPrice as bp
import numpy as np
import matplotlib.pyplot as plt

## Question 3

'''
Q3. a) European Call Option Prices via MC
'''
def mc_norm_payoff(const_drift, const_sigt, const_x_s, z):

    done = object()
    n1 = next(z, done)
    num = 0
    while (n1 is not done):
        yield max(0, math.exp(const_drift + const_sigt * n1) - const_x_s)
        yield max(0, math.exp(const_drift + const_sigt * -n1) - const_x_s)
        n1 = next(z, done)
        num += 1

def mc_callprice(S0, X, r, sigma, T):
    nldist = nl.NormalDistribution("Polar-Marsaglia")
    nsim = 1000
    nldist.generateNormalDistribution(nsim)
    z = nldist.getRdmNumbers()
    const_x_s = X / S0
    disc = math.exp(-r*T)
    const_drift = (r - sigma * sigma / 2.0) * T
    const_sigt = sigma * math.sqrt(T)
    l_payoff = list(mc_norm_payoff(const_drift, const_sigt, const_x_s, z))
    return disc*S0*(np.mean(l_payoff)),(disc*S0*np.std(l_payoff, ddof=1)/ math.sqrt(len(l_payoff)))


mc_callprice = mc_callprice(100,80,0.04,0.3,5)
upper_callprice = mc_callprice[0]+ 1.96*mc_callprice[1]
lower_callprice = mc_callprice[0]- 1.96*mc_callprice[1]
print('Q3. a) MC Call Price for (S,k,r,sigma,t) =(100,80,0.04,0.30,5) = %f with 95%% confidence in (%f,%f)' % (mc_callprice[0], lower_callprice,upper_callprice))


bpprice = bp.EuropeanCallPrice(100,80,0.04,0.3,5)
print('Q3. b) Black Scholes European Call Price for (S,k,r,sigma,t) =(100,80,0.04,0.30,5) = %f' % bpprice)



'''
Q3. c) Hedging parameters S0 = [15,25], X= 20 , sigma= 0.25 , r=0.04, T=0.5
'''

S0 = list(range(15,26))
X = 20
sigma = 0.25
r = 0.04
T = 0.5

delta = np.zeros(11)
gamma = np.zeros(11)
theta = np.zeros(11)
vega = np.zeros(11)
rho = np.zeros(11)

for i in range(0,11):
    delta[i] = bp.getDelta(S0[i],X,r,sigma,T)
    gamma[i] = bp.getGamma(S0[i],X,r,sigma,T)
    theta[i] = bp.getTheta(S0[i],X,r,sigma,T)
    vega[i] = bp.getVega(S0[i],X,r,sigma,T)
    rho[i] = bp.getRho(S0[i],X,r,sigma,T)


plt.plot(S0,delta)
plt.ylabel("Delta")
plt.xlabel("S")
plt.title("Delta vs S")
plt.show()


plt.plot(S0,gamma)
plt.ylabel("Gamma")
plt.xlabel("S")
plt.title("Gamma vs S")
plt.show()


plt.plot(S0,theta)
plt.ylabel("Theta")
plt.xlabel("S")
plt.title("Theta vs S")
plt.show()

plt.plot(S0,vega)
plt.ylabel("Vega")
plt.xlabel("S")
plt.title("Vega vs S")
plt.show()

plt.plot(S0,rho)
plt.ylabel("Rho")
plt.xlabel("S")
plt.title("Rho vs S")
plt.show()

