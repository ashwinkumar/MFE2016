import InterestRateModel as IRM
import BondPricing as bond
import numpy as np
import timeit
### Question 1

def EuropeanCallPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(St-X,array_zeros)

r0 = 0.05
sigma = 0.1
kappa = 0.82
rbar = 0.05

'''
(a) Use Monte Carlo Simulation (assume each time step is a day) to find the price of a pure discount bond, with Face Value of $1,000,
maturing in 𝑇=0.5 years (at time 𝑡=0):
'''
start_time1 = timeit.default_timer()

T = 0.5
dT = 1/365
nsims = 100000

FV=1000
bond_dT = 0.5
vasicek = IRM.Vasicek(r0= r0,r_bar=rbar,T=T, dT= dT, sigma= sigma, nsims = nsims, kappa= kappa)
r = vasicek.EulerDiscretize_R()
zeroCouponBond = bond.Bond(C0=0,FV=FV,T=T,dT=bond_dT, t0=0)
zeroCouponBond.setInterestRates(r,dT)
purediscountbondprice = zeroCouponBond.BondPrice()
elapsed1 = timeit.default_timer() - start_time1
print("1.a) Pure Discount Bond = %f      ****[%f sec]"%(purediscountbondprice,elapsed1))
del zeroCouponBond
del vasicek

'''
(b) Use Monte Carlo Simulation to find the price of a coupon paying bond, with Face Value of $1,000, paying semiannual
coupons of $30, maturing in 𝑇=4 years:
'''
start_time2 = timeit.default_timer()

T = 4
C0 = 30
vasicek = IRM.Vasicek(r0= r0,r_bar=rbar,T=T, dT= dT, sigma= sigma, nsims = nsims, kappa= kappa)
r = vasicek.EulerDiscretize_R()
CouponBond = bond.Bond(C0=C0,FV=FV,T=T,dT=bond_dT,t0=0)
CouponBond.setInterestRates(r,dT)
couponbondprice = CouponBond.BondPrice()
elapsed2 = timeit.default_timer() - start_time2
print("1.b) Coupon Bond = %f               ****[%f sec]"%(couponbondprice,elapsed2))
del CouponBond
del vasicek


'''
(c) Use Monte Carlo Simulation to find the price of a European Call option on the pure discount bond in part (a).
The option matures in 3 months and has a strike price of 𝐾=$980. Use the explicit formula for the underlying bond price (only for the bond price).
'''
start_time3 = timeit.default_timer()

T = 0.5
C0 = 0
nsims = 10000
t=0.25
K= 980

vasicek = IRM.Vasicek(r0= r0,r_bar=rbar,T=t, dT= dT, sigma= sigma, nsims = nsims, kappa= kappa)
r_sim1 = vasicek.EulerDiscretize_R()

vasicek = IRM.Vasicek(r0= r0,r_bar=rbar,T=T, dT= dT, sigma= sigma, nsims = nsims, kappa= kappa)
ZeroBond = bond.Bond(C0=C0,FV=FV,T=T,dT=bond_dT, t0=t)
T_array = ZeroBond.get_T_array()
A,B = vasicek.getA_B(t,T_array)

E_Call_Pure = ZeroBond.CallPrice_Semi_ClosedForm(A,B,r_sim1,K)
elapsed3 = timeit.default_timer() - start_time3

print("1.c) European Call Zero Coupon Bond = %f    ****[%f sec]"%(E_Call_Pure,elapsed3))
del vasicek
del ZeroBond


'''
(d) Use Monte Carlo Simulation to find the price of a European Call option on the coupon paying bond in part (b).
The option matures in 3 months and has a strike price of 𝐾=$980. Use Monte Carlo simulation for pricing the underlying bond.
'''
start_time4 = timeit.default_timer()
T = 4
C0 = 30
nsims1 = 500
nsims2 = 100
t=0.25
K= 980

vasicek = IRM.Vasicek(r0= r0,r_bar=rbar,T=t, dT= dT, sigma= sigma, nsims = nsims1, kappa= kappa)
r_sim1 = vasicek.EulerDiscretize_R()
r_avg = np.mean(r_sim1,1)
del vasicek
BP = np.zeros(nsims1)
for i in range(nsims1):
    vasicek = IRM.Vasicek(r0=r_sim1[i,-1], r_bar=rbar, T=(T-t), dT=dT, sigma=sigma, nsims=nsims2, kappa=kappa)
    r = vasicek.EulerDiscretize_R()
    CouponBond = bond.Bond(C0=C0,FV=FV,T=T,dT=bond_dT, t0=t)
    CouponBond.setInterestRates(r,dT)
    BP[i]= CouponBond.BondPrice()
    del vasicek
    del CouponBond
    del r
#print(np.mean(E_Call))
E_Call_Coupon = np.mean(EuropeanCallPayOff(BP,K)*np.exp(-r_avg*t))
elapsed4 = timeit.default_timer() - start_time4
print("1.d) European Call Coupon Bond = %f     ****[%f sec]"%(E_Call_Coupon,elapsed4))


'''
(e) Find the price of a European Call option of part (d) by using the explicit formula for the underlying bond price, and reconcile
 the findings with the ones of part (d).
'''
start_time5 = timeit.default_timer()
T = 4
C0 = 30
t=0.25
K= 980

vasicek = IRM.Vasicek(r0= r0,r_bar=rbar,T=t, dT= dT, sigma= sigma, nsims = nsims, kappa= kappa)
r_sim1 = vasicek.EulerDiscretize_R()
CouponBond = bond.Bond(C0=C0,FV=FV,T=T,dT=bond_dT, t0=t)
'''

E_Call_Pure_Full_ClosedForm = CouponBond.Vasicek_CallPrice_ClosedForm(vasicek,K,r0)

print("1.e) European Call Zero Coupon Bond(using r* method [Fully explicit] = %f"%E_Call_Pure_Full_ClosedForm)
'''
T_array = CouponBond.get_T_array()
A,B = vasicek.getA_B(t,T_array)
#A1,B1 = vasicek.getA_B()


E_Call_Pure_Semi_ClosedForm = CouponBond.CallPrice_Semi_ClosedForm(A,B,r_sim1,K)
elapsed5 = timeit.default_timer() - start_time5
print("1.e) European Call Zero Coupon Bond (explicit) = %f  ****[%f sec]"%(E_Call_Pure_Semi_ClosedForm,elapsed5))
