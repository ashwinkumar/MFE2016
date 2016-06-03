import numpy as np
import InterestRateModel as IRM
import  PrePaymentModel as PPM
import timeit
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
ybps = 5

#### Question 4
'''
Compute the OAS-adjusted Duration and Convexity of the MBS, considered in the previous question.
'''

MBSSecurity = 110000
start_time1 = timeit.default_timer()
cir = IRM.CIR(r0= r0,r_bar=rbar,T=T+10, dT= dT, sigma= sigma, nsims = nsims, kappa= kappa)
r_sim1 = cir.EulerDiscretize_R()

numerixPPM = PPM.NumerixPrepayment(PV0,WAC,r_sim1,T)
numerixPPM.setInterestRateModel(cir)
duration,convexity = numerixPPM.getDurationConvexity(ybps,MBSSecurity)

del cir
del numerixPPM
elapsed1 = timeit.default_timer() - start_time1
print('a) Duration = %f Convexity = %f     ****[%f sec]'%(duration,convexity,elapsed1))