import InterestRateModel as IRM
import PrePaymentModel as PPM
import timeit
import numpy as np

T = 30
WAC = 0.08
PV0 = 100000

### CIR Model
r0 = 0.078
kappa = 0.6
rbar = 0.08
sigma = 0.12
dT = 1/360
nsims = 1000


### Question 3
'''
Compute the Option-Adjusted-Spread (OAS) for the Numerix-Prepayment model case with the Market Price of MBS being $110,000.
'''

MBSSecurity = 110000
start_time1 = timeit.default_timer()
cir = IRM.CIR(r0= r0,r_bar=rbar,T=T+10, dT= dT, sigma= sigma, nsims = nsims, kappa= kappa)
r_sim1 = cir.EulerDiscretize_R()

numerixPPM = PPM.NumerixPrepayment(PV0,WAC,r_sim1,T)
numerixPPM.setInterestRateModel(cir)
oasspread = numerixPPM.OASSpread(MBSSecurity)

del cir
del numerixPPM
elapsed1 = timeit.default_timer() - start_time1
print('a) OAS = %f      ****[%f sec]'%(oasspread,elapsed1))