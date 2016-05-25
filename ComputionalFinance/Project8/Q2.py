import InterestRateModel as IRM
import BondPricing as bond
import numpy as np
import timeit
import FiniteDifferenceMethod as FDM

def EuropeanCallPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(St-X,array_zeros)


###  Question 2 CIR Model
r0 = 0.05
sigma = 0.12
kappa = 0.92
rbar = 0.055

'''
(a) Use Monte Carlo Simulation to find at time 𝑡=0 the price 𝑐(𝑡,𝑇,𝑆) of a European Call option, with strike price of 𝐾=$980,
 maturity of 𝑇=0.5 years on a Pure Discount Bond with Face Value of $1,000, that matures in 𝑆=1 year:
'''
start_time1 = timeit.default_timer()
T = 1
dT = 1/365
nsims1 = 1000
nsims2 = 500

FV=1000
C0 = 0
bond_dT = T
t = 0.5
K= 980
cir = IRM.CIR(r0= r0,r_bar=rbar,T=t, dT= dT, sigma= sigma, nsims = nsims1, kappa= kappa)
r_sim1 = cir.EulerDiscretize_R()

r_avg = np.mean(r_sim1,1)
del cir
BP = np.zeros(nsims1)
for i in range(nsims1):
    cir = IRM.CIR(r0=r_sim1[i,-1], r_bar=rbar, T=(T-t), dT=dT, sigma=sigma, nsims=nsims2, kappa=kappa)
    r = cir.EulerDiscretize_R()
    ZeroCouponBond = bond.Bond(C0=C0,FV=FV,T=T,dT=bond_dT, t0=t)
    ZeroCouponBond.setInterestRates(r,dT)
    BP[i]= ZeroCouponBond.BondPrice()
    del cir
    del ZeroCouponBond
    del r

E_Call_Coupon = np.mean(EuropeanCallPayOff(BP,K)*np.exp(-r_avg*t))

elapsed_time1 = timeit.default_timer() - start_time1
print("2.a) European Call Zero Coupon Bond = %f    ****[%f sec]"%(E_Call_Coupon,elapsed_time1))


'''
(b) Use the Implicit Finite-Difference Method to find at time 𝑡=0 the price 𝑐(𝑡,𝑇,𝑆) of a European Call option, with strike price of 𝐾=$980,
 maturity of 𝑇=0.5 years on a Pure Discount Bond with Face Value of $1,000, that matures in 𝑆=1 year.
'''
start_time2 = timeit.default_timer()
rmin = 0
rmax = 0.12

dT = 0.001
dr = 0.01
#dr= 0.0001

fdm_class = FDM.CIR_ImplicitFiniteDifferenceMethod(r0,rbar,kappa,sigma,t,dT)
r_fdm,ind = fdm_class.initRate(rmin,rmax,dr)

cir = IRM.CIR(r0= r0,r_bar=rbar,T=t, dT= dT, sigma= sigma, nsims = nsims1, kappa= kappa)
BP = FV*cir.BondPrice_ClosedForm(r_fdm,t,T)
del cir

P =fdm_class.getCallOptionPrice(BP,K)
elapsed_time2 = timeit.default_timer() - start_time2

print("2.b) Implicit Differential Method European Call Zero Coupon Bond = %f   ****[%f sec]"%(P[ind],elapsed_time2))

'''
(c) Compute the price 𝑐(𝑡,𝑇,𝑆) of the European Call option above using the explicit formula, and compare it to your findings in parts (a) and (b)
 and comment on your findings.
'''
start_time3 = timeit.default_timer()

nsims=10000

cir = IRM.CIR(r0= r0,r_bar=rbar,T=t, dT= dT, sigma= sigma, nsims = nsims, kappa= kappa)
ZeroCoupon = bond.Bond(C0=C0,FV=FV,T=T,dT=bond_dT, t0=t)
E_Call_Pure_Full_ClosedForm = ZeroCoupon.CIR_CallPrice_ClosedForm(cir,K/FV,r0)
elapsed_time3=  timeit.default_timer() - start_time3
print("2.c) European Call Zero Coupon Bond = %f   ****[%f sec]"%(E_Call_Pure_Full_ClosedForm,elapsed_time3))
