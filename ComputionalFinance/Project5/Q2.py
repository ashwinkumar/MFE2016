from MonteCarlo import  MCEuropeanOption
import numpy as np
import math
import timeit
import LeastSquareMC as lsmc
#### Question 2


'''
def myEuropeanPayOff(St):
    if type(St).__module__ != np.__name__ :
        raise Exception('St is not numpy array')
    l = len(St)
    return np.max(St)-St[l-1]
'''
def putInTheMoney(St,X):
    return np.where(St <X)

def PutPayOff(St,X):
    if isinstance(St,np.float):
        return max(0,X-St)
    array_zeros = np.zeros(len(St))
    return np.maximum(X-St,array_zeros)
S0 = 65
K = 60
sigma = 0.2
r = 0.06
t = 0.2
T = 1
nsims = 100000

nsims1 = 100
nsims2 = (int)(nsims/nsims1)
ndiv = 100
dT_t =(T-t)/ndiv
dt= t/ndiv

'''
a) Forward start European Put Option
'''
start_time1 = timeit.default_timer()
myMCEuropean = MCEuropeanOption(S0,r,sigma,t,dt,nsims1)
myMCEuropean.mcStockTree()
X = myMCEuropean.getStrikePrice(t)
FS_EuropeanPutPrice_sim = np.zeros(nsims1)
#FS_Put_var = np.zeros(nsims1)
for i in range(nsims1):
    myMCEuropean_branch = MCEuropeanOption(X[i],r,sigma,T-t,dT_t,nsims2)
    myMCEuropean_branch.mcStockTree()
    myMCEuropean_branch.setPayOff(PutPayOff)
    FS_EuropeanPutPrice_sim[i],var= myMCEuropean_branch.getOptionPrice(X[i],0)
    del myMCEuropean_branch
elapsed_time1 = timeit.default_timer() - start_time1
FS_EuropeanPutPrice = np.mean(FS_EuropeanPutPrice_sim)*math.exp(-r*t)
FS_EuropeanPut_var = math.exp(-2*r*t)*np.var(FS_EuropeanPutPrice_sim, ddof=1)/nsims1
print('Forward-start European Put Price = %f  with var = %f   ****[%f sec]' % (FS_EuropeanPutPrice,FS_EuropeanPut_var,elapsed_time1))

'''
b) Forward start American Put Option
'''

k =2
start_time2 = timeit.default_timer()
FS_AmericanPutPrice_sim = np.zeros(nsims1)

for i in range(nsims1):
    lsclass = lsmc.LeastSquareMC(X[i], r, sigma, T-t, dT_t, nsims2, k, "Monomials")
    lsclass.setPayOff(PutPayOff)
    lsclass.setInMoney(putInTheMoney)
    lsclass.mcStockTree()
    X_array = X[i]*np.ones(nsims2)
    FS_AmericanPutPrice_sim[i] = lsclass.calcualteOptionPrice(X_array, 0)
    del lsclass

elapsed_time2 = timeit.default_timer() - start_time2
FS_AmericanPutPrice = np.mean(FS_AmericanPutPrice_sim)*math.exp(-r*t)
FS_AmericanPut_var = math.exp(-2*r*t)*np.var(FS_AmericanPutPrice_sim, ddof=1)/nsims1
print('Forward-start American Put Price = %f  with var = %f   ****[%f sec]' % (FS_AmericanPutPrice,FS_AmericanPut_var,elapsed_time2))
