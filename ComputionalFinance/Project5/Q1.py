import LeastSquareMC as lsmc
import numpy as np
import timeit
#### Queston 1

def putInTheMoney(St,X):
    return np.where(St <X)

def PutPayOff(St,X):
    if isinstance(St,np.float)  or isinstance(X,np.float):
        return max(0,X-St)
    array_zeros = np.zeros(len(St))
    return np.maximum(X-St,array_zeros)


S0_array = np.array([36,40,44])
X = 40
#t_array = np.zeros(3) # Current time is t=0
T_array = np.array([0.5,1,2])
r = 0.06
sigma = 0.2
nsims = 100000
ndiv = 100
dT_array = T_array/ndiv
k_array = np.array([2,3,4], dtype= np.int) # k=2 is OLS

method_array = ["Laguerre","Hermite", "Monomials"]
sub = ['a','b','c']

def myRoutineQ1(S0,T,dT,k,method):
    start_time1 = timeit.default_timer()
    lsclass = lsmc.LeastSquareMC(S0,r,sigma,T,dT,nsims,k,method)
    lsclass.setPayOff(PutPayOff)
    lsclass.setInMoney(putInTheMoney)
    lsclass.mcStockTree()
    #elapsed1 = timeit.default_timer() - start_time1
    #print(' Time elapsed in generate the stock tree %f seconds' %elapsed1)
    #start_time2 = timeit.default_timer()
    optionPrice =lsclass.calcualteOptionPrice1(X)
    elapsed1 = timeit.default_timer() - start_time1
    print('American Put OptionPrice (S0,T)=[%f,%f]  = %f  [Method = %s]    ****[%f sec]' %(S0,T,optionPrice,method,elapsed1))

for j in range(len(method_array)):
    print('Q1.%s ' %sub[j])
    for i in range(len(S0_array)):
        myRoutineQ1(S0_array[i],T_array[i],dT_array[i],k_array[i],method_array[j])
