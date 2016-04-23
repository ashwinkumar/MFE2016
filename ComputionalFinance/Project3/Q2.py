import EulerDiscretization as ed
import mynormal as nl
import numpy as np
import timeit
import mywrapper as wp
import math
##Question 2
'''
Estimate the following expected values:
    𝐸(1+𝑋3)13⁄,𝐸(1+𝑌3)13⁄
'''

X0 = np.array([1])
ndiv = 1000
T = 3
nsims = 100


def coeff_xdw(x,y,t,z,cdw):
    z1 = z[0]
    w1 = z[1]
    return (1/3)*x*z1 - (3/4)*x*w1;

def coeff_xdt(x,y,t,cdt):
    return x/4



X = np.zeros([nsims, ndiv+1])
start_time1 = timeit.default_timer()

for i in range(0,nsims):
    X[i,:] = ed.EulerDiscretization(X0,[],[], coeff_xdt, [],coeff_xdw ,[], T, ndiv,2)

res = np.mean(wp.mycuberoot(1 + X[:,ndiv ]))


elapsed1 = timeit.default_timer() - start_time1
print('Q2.  a) Expection of the expression = %f                                          ****[%f sec]' % (res,elapsed1))


'''
Q2. b) Y = exp(-0.08t + 1/3*Wt + 3/4*Zt)
'''


def myMC(n1,n2,t,sqr_t):
    done = object()
    z1 = next(n1, done)
    z2 = next(n2, done)
    while (z1 is not done) and (z2 is not done):
        yield (1 + math.exp(-0.08*t + sqr_t*(z1/3 + 3*z2/4)))
        z1 = next(n1, done)
        z2 = next(n1, done)


start_time2 = timeit.default_timer()

T = 3
sqrt_t = math.sqrt(T)
nsims = 10000

normal1  = nl.NormalDistribution("Polar-Marsaglia",14)
normal2 = nl.NormalDistribution("Polar-Marsaglia",18)
normal1.generateNormalDistribution(nsims)
normal2.generateNormalDistribution(nsims)
n1= normal1.getRdmNumbers()
n2= normal2.getRdmNumbers()


res= wp.mycuberoot(np.asarray(list(myMC(n1,n2,T,sqrt_t))))
res_mean = np.mean(res)
res_var = np.var(res, ddof=1)
elapsed2 = timeit.default_timer() - start_time2

print('     b) Expected value of the function E[(1+Y)^(1/3)] = %f with var =%f     ****[%f sec]' %(res_mean,res_var,elapsed2))



