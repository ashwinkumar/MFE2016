import InterestRateModel as IRM
import PrePaymentModel as PPM
import numpy as np
import timeit
import matplotlib.pyplot as plt

import math

T = 30
WAC = 0.08
PV0 = 100000

### CIR Model
r0 = 0.078
kappa = 0.6
rbar = 0.08
sigma = 0.12
dT = 1/365
nsims = 10000


'''
(a) Compute the price of the MBS using this model for prepayments. The code should be generic: the user is prompted for
inputs and the program runs and gives the output.
'''
start_time1 = timeit.default_timer()
cir = IRM.CIR(r0= r0,r_bar=rbar,T=T+10, dT= dT, sigma= sigma, nsims = nsims, kappa= kappa)
r_sim1 = cir.EulerDiscretize_R()

numerixPPM = PPM.NumerixPrepayment(PV0,WAC,r_sim1,T)
MBSPrice = numerixPPM.MBSPrice()
del cir
del numerixPPM
elapsed1 = timeit.default_timer() - start_time1
print('a) MBS = %f      ****[%f sec]'%(MBSPrice,elapsed1))

'''
(b) Compute the price of the MBS for the following ranges of the parameters: 𝑘 in 0.3 to 0.9 (in increments of 0.1)
and draw the graph of the price vs. 𝑘.
'''
kappa_start = 0.3
kappa_end = 0.9
dk = 0.1
kappa_array = np.arange(kappa_start,kappa_end+dk,dk)
MBSPrice = np.zeros(len(kappa_array))
for i in range(len(kappa_array)):
    cir = IRM.CIR(r0=r0, r_bar=rbar, T=T + 10, dT=dT, sigma=sigma, nsims=nsims, kappa=kappa_array[i])
    r_sim1 = cir.EulerDiscretize_R()
    numerixPPM = PPM.NumerixPrepayment(PV0, WAC, r_sim1, T)
    MBSPrice[i] = numerixPPM.MBSPrice()
    del cir
    del numerixPPM

plt.xlabel("Kappa")
plt.ylabel("MBS Price")
plt.title("MBS Price vs Kappa of CIR Model")
plt.plot(kappa_array,MBSPrice)
plt.savefig("1bfig1")
print(MBSPrice)


'''
(c) Compute the price of the MBS for the following ranges of the parameters: 𝑟̅ in 0.03 to 0.09 (in increments of 0.01)
and draw the graph of the price vs. 𝑟̅.
'''

rbar_start = 0.03
rbar_end = 0.09
dr = 0.01
rbar_array = np.arange(rbar_start,rbar_end+dr,dr)
MBSPrice = np.zeros(len(rbar_array))
for i in range(len(rbar_array)):
    cir = IRM.CIR(r0=r0, r_bar=rbar_array[i], T=T + 10, dT=dT, sigma=sigma, nsims=nsims, kappa=kappa)
    r_sim1 = cir.EulerDiscretize_R()
    numerixPPM = PPM.NumerixPrepayment(PV0, WAC, r_sim1, T)
    MBSPrice[i] = numerixPPM.MBSPrice()
    del cir
    del numerixPPM
plt.clf()

plt.xlabel("rbar")
plt.ylabel("MBS Price")
plt.title("MBS Price vs rbar of CIR Model")
plt.plot(rbar_array,MBSPrice)
plt.savefig("1cfig2")