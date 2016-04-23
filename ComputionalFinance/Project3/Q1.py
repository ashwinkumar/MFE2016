import EulerDiscretization as ed
import numpy as np
import timeit
import mywrapper as wp
import mynormal as nl
## Question 1

'''
Evaluate the expected values and probabilities
'''

X0 = np.array([1])
Y0 = np.array([3/4])
ndiv = 1000
T = 4
nsims = 1000


def coeff_xdw(x,x1,t,z,cdw):
    return 2/3*z

def coeff_xdt(x,x1,t,cdt):
    return 1/5 - x/2


def coeff_ydw(y,y1,t,z,cdw):
    return (1+ t*t*t)*z/3

def coeff_ydt(y,y1,t,cdt):
    return ((y*2/(1+t)) + (1 + t*t*t)/3);

def getProb(Y,t,p):
    x = Y[:,t]
    s = sum(x >p)
    return s/len(x)

X = np.zeros([nsims, ndiv+1])
Y = np.zeros([nsims,ndiv+1])
start_time1 = timeit.default_timer()

for i in range(0,nsims):
    X[i,:] = ed.EulerDiscretization(X0,[],[], coeff_xdt, [] ,coeff_xdw ,[],T, ndiv,1)
    Y[i,:] = ed.EulerDiscretization(Y0,[],[], coeff_ydt, [],coeff_ydw,[],T, ndiv,1)

t1 = 2
t2 = 3
dt1 = (int) (ndiv/T*t1)
dt2 = (int) (ndiv/T*t2)
vl = 5
p1 = getProb(Y, dt1,vl)
elapsed1 = timeit.default_timer() - start_time1
print('Q1. a) P(Y2 > 5) = %f                        ****[%f sec]' % (p1,elapsed1))


start_time2 = timeit.default_timer()
exp1 = np.mean(wp.mycuberoot(X[:,dt1]))
elapsed2 = timeit.default_timer() - start_time2

print('Q1. b) Expected value of X2^(1/3) = %f       ****[%f sec]' % (exp1, elapsed2))

start_time3 = timeit.default_timer()
exp2 = np.mean(Y[:,dt2])
elapsed3 = timeit.default_timer() - start_time3
print('Q1. c) Expected value of Y3 = %f             ****[%f sec]' % (exp2,elapsed3))

start_time4 = timeit.default_timer()
X1 = X[:,dt1]
ncount = sum(X1 >=1)
X1[X1 <=1] = 0
Y1 = Y[:,dt1]
exp3 = np.mean(X1*Y1)
exp4 = sum(X1*Y1)/ncount
elapsed4 = timeit.default_timer() - start_time4
print('Q1. d) Expected value of X2*Y2 , X2>1 = %f    ****[%f sec]' % (exp3, elapsed4))

print('Q1. d) Expected value (conditional) of X2*Y2 | X2>1 = %f    ****[%f sec]' % (exp4, elapsed4))
