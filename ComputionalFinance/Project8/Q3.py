import InterestRateModel as IRM
import numpy as np
import BondPricing as bond
import timeit
#### Question 3 G2++ Model

def EuropeanPutPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(X-St,array_zeros)

x0 = 0
y0 = 0
phi0 = 0.03
phit = 0.03
r0  = 0.03
rho = 0.7
a  = 0.1
b = 0.3
sigma = 0.03
eta = 0.08
nsims1 = 500
nsims2 = 200




T = 1
dT = 1/365

FV=1000
C0 = 0
bond_dT = T
t = 0.5
K= 950

start_time1 = timeit.default_timer()
g2pp = IRM.G2pp(r0,t,dT,sigma,phi0,phit,rho,nsims1)
g2pp.initX(x0,y0)
g2pp.initY(y0,b,eta)
r_sim,x_sim,y_sim = g2pp.EulerDiscretize_R()

r_avg = np.mean(r_sim,1)
del g2pp

BP = np.zeros(nsims1)

for i in range(nsims1):
    g2pp = IRM.G2pp(r0,T-t,dT,sigma,phi0,phit,rho,nsims2)
    g2pp.initX(x_sim[i,-1],a)
    g2pp.initY(y_sim[i,-1],b,eta)
    r,_,_ = g2pp.EulerDiscretize_R()
    ZeroCouponBond = bond.Bond(C0=C0,FV=FV,T=T,dT=bond_dT, t0=t)
    ZeroCouponBond.setInterestRates(r,dT)
    BP[i]= ZeroCouponBond.BondPrice()
    del g2pp
    del ZeroCouponBond
    del r

E_Put_Coupon = np.mean(EuropeanPutPayOff(BP,K)*np.exp(-r_avg*t))
elapsed_time1 = timeit.default_timer()- start_time1
print("3.a) MC Simulation European Put Zero Coupon Bond = %f    ****[%f sec]"%(E_Put_Coupon,elapsed_time1))


#### Explicit method

start_time2 = timeit.default_timer()

nsims=100

g2pp = IRM.G2pp(r0,t,dT,sigma,phi0,phit,rho,nsims1)
g2pp.initX(x0, a)
g2pp.initY(y0, b, eta)
ZeroCoupon = bond.Bond(C0=C0,FV=FV,T=T,dT=bond_dT, t0=t)
E_Put_Pure_Full_ClosedForm = ZeroCoupon.G2pp_CallPrice_ClosedForm(g2pp,K/FV)
elapsed_time2=  timeit.default_timer() - start_time2
print("3.b) European Put Zero Coupon Bond = %f   ****[%f sec]"%(E_Put_Pure_Full_ClosedForm,elapsed_time2))


