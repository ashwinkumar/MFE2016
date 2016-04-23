
import mynormal as nl
import math
import statistics as st
import numpy as np
import MonteCarlo as mc
import mywrapper as wp

CONST_MAX = 1000
## Question 3
'''
(a) Estimate the following expected values by simulation where 𝑊𝑡 is a Standard Wiener Process. :
    1.  𝛦(𝑊52+𝑠𝑖𝑛(𝑊5))  2. 𝛦(𝑒𝑡2𝑐𝑜𝑠(𝑊𝑡))for 𝑡=0.5,3.2,6.5.
'''

nsim = 10000
nldist = nl.NormalDistribution("Polar-Marsaglia")
nldist.generateNormalDistribution(nsim)
z = nldist.getRdmNumbers()
t = 5
sqrtt = math.sqrt(t)
mc1 = list(mc.myWienerMC1(z, sqrtt))
res1 = st.mean(mc1)
var1 = st.variance(mc1)/len(mc1)
print('Q3. a) The value of expression by MC Simulaiton is %f with var = %f'  % (res1, var1))

list_t = np.array([0.5,3.2,6.5])
sqrt_t = np.sqrt(list_t)
exp_t = np.exp(list_t/2.0)
exp_t2 = np.exp(list_t)
nldist.generateNormalDistribution(nsim)
z = nldist.getRdmNumbers()

mc2 = list(mc.myCos(z, sqrt_t))
res2 = np.multiply(np.mean(mc2, axis = 0), exp_t)
var2 = np.multiply(np.var(mc2, axis=0), exp_t2)/(len(mc2))
print('       The value of individual expressions with t = (0.5,3.2,6.5) is %s with var = %s' % (res2, var2))
print('Q3. b) The E[cos(Wt)] = exp(-t/2). Hence the expression should evaluate to 1 for all values of t as nsims increses. This is also seen in our MC simulation where the expectation~ 1')

## For variance reduction we use +/- z
nldist.generateNormalDistribution(nsim/2)
z = nldist.getRdmNumbers()
mc1_vr = list(mc.myWienerMC1_vr(z, sqrtt))
res1_vr = st.mean(mc1_vr)
var1_vr = st.variance(mc1_vr) / len(mc1_vr)

print('Q3. c) The value of expression 1 by MC Simulaiton (reduced variance) is %f with var = %f'  % (res1_vr, var1_vr))
perred = (var1 /var1_vr -1)*100
print('       The variance reduction = %f%%' % perred)
nldist.generateNormalDistribution(nsim/2)
z = nldist.getRdmNumbers()

# Variance reduction using control variate method take Y = z1*z1 where z1(0,1)

def gencorrarrays(z,sqrt_t):
    done = object()
    n1 = next(z, done)
    while (n1 is not done):
        tmp1 = math.cos(n1 * sqrt_t)
        tmp2 = n1 * n1 *sqrt_t
        yield tmp1 ,tmp2
        n1 = next(z, done)

def getGamma(z, sqrt_t):
    mat = list(gencorrarrays(z, sqrt_t))
    X1 = np.asarray([item[0] for item in mat])
    Y1 = np.asarray([item[1] for item in mat])
    cov = wp.cov(X1,Y1)
    v = np.var(Y1)
    return cov/v

nsim_loc = 1000
nldist.generateNormalDistribution(nsim_loc)
z = nldist.getRdmNumbers()
gamma = getGamma(z, sqrtt)
print('       The gamma for the functions (cos(n1*sqrt(t),n1*n1) is %f ' %gamma)
nldist.generateNormalDistribution(nsim)
z = nldist.getRdmNumbers()

mc2_vr = list(mc.myCos_vr(z,gamma , sqrt_t))

res2_vr = np.multiply(np.mean(mc2_vr, axis = 0), exp_t)
var2_vr = np.multiply(np.var(mc2_vr, axis=0), exp_t2)/len(mc2_vr)

print('Q3. c) The value of expression 2 with t =  (0.5,3.2,6.5)  by MC Simulaiton (reduced variance) is %s with var = %s'  % (res2_vr, var2_vr))

perred2 = (1- np.multiply(var2_vr,1/var2))*100
print('       The variance reduction in %%', perred2)