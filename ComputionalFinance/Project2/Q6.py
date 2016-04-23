import EulerIntegral as ei
import myrandgen as rnd
import math
import statistics as st
from scipy.optimize import fmin

import numpy as np

###Question 1
'''
a) Euler integral for given function
'''

def f(x):
    return 4* math.sqrt(1 - x*x)

def MC_pi(u):
    done = object()
    n1 = next(u, done)
    while (n1 is not done):
        yield f(n1)
        n1 = next(u, done)

n = 10000
e1 = ei.EulerIntegral(f, 0, 1, n)
print('Q6. a)The value of pi by Euler Integration %f' % e1)

'''
b) Monte Carlo Simulation for evaluating pi
'''
uniform = rnd.LGMRandomGenerator(14)
u = uniform.generateRdmNumberByGenerator(n)
mc1= list(MC_pi(u))
e2 = st.mean(mc1)
stde2 = st.stdev(mc1) / math.sqrt(n)
print('Q6. b)The value of pi from MC = %f with std dev %f'  % (e2 ,stde2))

'''
c) Importance sampling
'''

def g(x,alpha):
    return (1- alpha * x * x) / (1 - alpha/3)

def h(x,alpha):
    return f(x)/g(x,alpha)
'''
def minvar(a):
    tanh = math.atanh(math.pow(a,-0.5))
    return  (1- a/3)*(1/a - (1-a)*tanh/math.pow(a,1.5))

def mc_min_alpha():
    return fmin(minvar, 0.74)
'''


alpha_min = 0.74

gmax = 1/ (1- alpha_min/3)

def generate_from_g(n, urand, g , alpha_min,gmax):
    num = 0
    while num < n:
        Y = urand.getNextRdmNumber()
        u = urand.getNextRdmNumber()
        if Y >1 or Y < 0:
            print('Uniform (0,1) generating out of range. Check myrandge')
            raise AttributeError
        x = g(Y,alpha_min)/gmax
        if u <= x :
            yield Y
            num = num+1

def pi_varred(n, urand, f,g,alpha_min, gmax):
    sum = 0
    X = list(generate_from_g(n,urand,g,alpha_min,gmax))
    res = np.zeros(n, dtype = float)
    for i in range(0,len(X)):
        res[i] = h(X[i],alpha_min)
    return np.mean(res), np.std(res,ddof=1)/math.sqrt(n)

urand = rnd.LGMRandomGenerator(14)
pi_var_rd = pi_varred(n, urand, f, g , alpha_min, gmax)

print('Q6. c)The value of pi from Importance Sampling = %f with std dev = %f'  % (pi_var_rd[0], pi_var_rd[1]))
perred =  math.fabs(math.pow((pi_var_rd[1])/stde2,2) -1)*100
print('      We see that the variance has reduced by %f %% using Importance Sampling '  % perred)
















