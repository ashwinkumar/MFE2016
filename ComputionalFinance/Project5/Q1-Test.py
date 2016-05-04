import BinomialModel as bm
import math
import numpy as np
import matplotlib.pyplot as plt


def PutPayOff(St,X):
    if isinstance(St,np.float)  or isinstance(X,np.float):
        return max(0,X-St)
    array_zeros = np.zeros(len(St))
    return np.maximum(X-St,array_zeros)

S0_array = np.array([36,40,44])

T_array = np.array([0.5,1,2])
ndiv = 100
dT_array =  T_array/ndiv

r = 0.06
sigma = 0.2
nsims = 100000
X = 40

put_price = np.zeros([3,3])

for k in range(len(dT_array)):
    u, d, p = bm.getBinomialParameters(r, sigma, dT_array, 3)
    for j in range(len(S0_array)):
        bin_model = bm.BinomialModel(S0_array[j], r, T_array[k])
        bin_model.n = ndiv
        if bin_model.verify(u[k], d[k], p[k]):
            put_price[j,k] = bin_model.getOptionPrice(PutPayOff, X, "A")
        del bin_model

print('American Put_price for S0 = %f  t =(%s) = [%s] ' %(S0_array[0],T_array,put_price[0,:]))
print('American Put_price for S0 = %f  t =(%s) = [%s] ' %(S0_array[1],T_array,put_price[1,:]))
print('American Put_price for S0 = %f  t =(%s) = [%s] ' %(S0_array[2],T_array,put_price[2,:]))






