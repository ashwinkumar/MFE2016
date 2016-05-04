import LeastSquareMC as lsmc
import numpy as np
import timeit
import BinomialModel as bm
import  math
from  tabulate import  tabulate
#### Queston 1

def putInTheMoney(St,X):
    return np.where(St <X)

def PutPayOff(St,X):
    if isinstance(St,np.float)  or isinstance(X,np.float):
        return max(0,X-St)
    array_zeros = np.zeros(len(St))
    return np.maximum(X-St,array_zeros)


S0_array = np.array([36,40,44])
#t_array = np.zeros(3) # Current time is t=0
T_array = np.array([0.5,1,2])
r = 0.06
sigma = 0.2
nsims = 10000

ndiv = 100
dT_array = T_array/ndiv

X = 40

bin_put_price = np.zeros([len(S0_array),len(dT_array)])

for k in range(len(dT_array)):
    u, d, p = bm.getBinomialParameters(r, sigma, dT_array, 3)
    for j in range(len(S0_array)):
        bin_model = bm.BinomialModel(S0_array[j], r, T_array[k])
        bin_model.n = ndiv
        if bin_model.verify(u[k], d[k], p[k]):
            bin_put_price[j,k] = bin_model.getOptionPrice(PutPayOff, X, "A")
        del bin_model

X = 40*np.ones(nsims)

k_array = np.array([2,3,4], dtype= np.int) # k=2 is OLS

method_array = ["Laguerre","Hermite", "Monomials"]
sub = ['a','b','c']

def myRoutineQ1(S0,T,dT,k,method):
    #start_time1 = timeit.default_timer()
    lsclass = lsmc.LeastSquareMC(S0,r,sigma,T,dT,nsims,k,method)
    lsclass.setPayOff(PutPayOff)
    lsclass.setInMoney(putInTheMoney)
    lsclass.mcStockTree()
    #elapsed1 = timeit.default_timer() - start_time1
    #print(' Time elapsed in generate the stock tree %f seconds' %elapsed1)
    start_time2 = timeit.default_timer()
    optionPrice =lsclass.calcualteOptionPrice(X,0)
    #elapsed1 = timeit.default_timer() - start_time1
    #print('American Put OptionPrice (S0,T,K)=[%f,%f,%f]  = %f  [Method = %s]    ****[%f sec]' %(S0,T,k,optionPrice,method,elapsed1))
    return optionPrice

for i in range(len(method_array)):
    print('Q1.%s [method = %s]' %(sub[i],method_array[i]))
    for j in range(len(S0_array)):
        print(' S0 = %f' %S0_array[j])
        result = [[["T/K","nan"],["Binomial-using J-R Model","% deviation"],["K=2","% deviation"], ["K=3","% deviation"],["K=4","% deviation"]]]
        for k in range(len(dT_array)):
            sub_res = np.zeros((len(k_array)+2,2))
            sub_res[0,0]= T_array[k]
            sub_res[0,1]= None
            sub_res[1,0] = bin_put_price[j,k]
            sub_res[1,1] = 0

            for l in range(len(k_array)):
                sub_res[l+2,0]=myRoutineQ1(S0_array[j],T_array[k],dT_array[k],k_array[l],method_array[i])
                sub_res[l+2,1] = math.fabs(sub_res[l+2,0]-sub_res[1,0])*100/ sub_res[1,0]
            result.append(sub_res)
            del sub_res
        print(tabulate(result))
