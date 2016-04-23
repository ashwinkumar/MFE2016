import numpy as np
import math
import random as rndom
import mynormal as nl

'''
Note: Normal Random number generator takes random seed
Generalized Euler Discretization
coeff_dt -> generic function which can depend on x,t,y1
coeff_dq -> generic function which can depend on x,t,y1,dw1
'''
def EulerDiscretization(x0,cdt,cdw ,coeff_dt, y1, coeff_dw,y2, T, ndiv, dim_z):
    num = 1
    dt = T/ ndiv
    dim_x = x0.shape[0]
    x = np.zeros((dim_x,ndiv+1))
    x[:,0]= x0

    mynormal = nl.NormalDistribution("Polar-Marsaglia", rndom.randint(10,10*ndiv))
    mynormal.generateNormalDistribution(ndiv*dim_z)
    ndist = mynormal.getRdmNumbers()
    sqrt_t = math.sqrt(dt)
    done = object()
    z = np.zeros(dim_z)
    n= next(ndist, done)
    while (n is not done):
        k=0
        while (k< dim_z) and (n is not done):
            z[k] = n
            k= k+1
            n = next(ndist, done)

        x[:,num] = x[:,num - 1] + coeff_dt(x[:,num-1], y1,dt*num, cdt) * dt + coeff_dw(x[:,num-1],y2,dt*num,z,cdw)*sqrt_t
        num = num+1
    return x


