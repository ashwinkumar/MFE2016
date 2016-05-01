from MonteCarlo import  MCEuropeanOption
import numpy as np
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
    if isinstance(St,np.float)  or isinstance(X,np.float):
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
ndiv = 100
dT =T/ndiv

'''
a) Forward start European Put Option
'''
start_time1 = timeit.default_timer()
myMCEuropean = MCEuropeanOption(S0,r,sigma,T,dT,nsims)
myMCEuropean.mcStockTree()
X = myMCEuropean.getStrikePrice(t)
myMCEuropean.setPayOff(PutPayOff)
FS_EuropeanPutPrice, FS_Put_var= myMCEuropean.getOptionPrice(X,t)
elapsed_time1 = timeit.default_timer() - start_time1
print('Forward-start European Put Price = %f  with var = %f   ****[%f sec]' % (FS_EuropeanPutPrice,FS_Put_var,elapsed_time1))

'''
b) Forward start American Put Option
'''
k =2
start_time2 = timeit.default_timer()
lsclass = lsmc.LeastSquareMC(S0, r, sigma, T, dT, nsims, k, "Monomials")
lsclass.setPayOff(PutPayOff)
lsclass.setInMoney(putInTheMoney)
lsclass.mcStockTree()
X =lsclass.getStrikePrice(t)
# elapsed1 = timeit.default_timer() - start_time1
# print(' Time elapsed in generate the stock tree %f seconds' %elapsed1)
# start_time2 = timeit.default_timer()
optionPrice = lsclass.calcualteOptionPrice(X,t)
elapsed_time2 = timeit.default_timer() - start_time2
print('Forward-start American Put Price = %f    ****[%f sec]' % (optionPrice,elapsed_time2))
