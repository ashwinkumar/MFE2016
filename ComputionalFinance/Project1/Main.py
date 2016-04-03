import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom
from array import *


import myrandgen as rndgen
import mybernoulli as bl
import mybinomial as bn
import mywrapper as wp


## Question 1

'''
(a) Using LGM method generate 10,000 Uniformly distributed random numbers on [0,1] and
compute the empirical mean and the standard deviation of the sequence.
'''
lgm = rndgen.LGMRandomGenerator(14)
lgm.generateRdmNumber(10000)
print('Q1. a) LGM Mean %f' % lgm.meanRdmNumber())
print('       LGM Std Dev %f' % lgm.stdDeviationRdmNumber())

'''
(b) Now use built-in functions of whatever software you are using to do the same thing as in (a).
'''
x_array = np.random.uniform (0, 1, 10000)
print('Q1. b) Built in Fn Mean %f' % x_array.mean())
print('       Built in Fn Std Dev %f' % x_array.std())



## Question 2

'''
(a)
Generate 10,000 random numbers with the following distribution:
'''
p = np.array((0.3, 0.35, 0.2, 0.15 ))
X = np.array((-1, 0, 1, 2))
bernouli = bl.BernoulliDistribution(p, X)
bernouli.generateBernoulliDistribution(10000)

'''
(b)
Draw the histogram and compute the empirical mean and the standard deviation of the sequence of 10,000 numbers generated above in part (a).
'''

bdist = bernouli.getBernoulDistNumbers()
plt.hist(bdist)
plt.title("Histogram")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()
print('Q2 b) Mean of given dist %f' % bdist.meanRdmNumber())
print('Q2 b) Std deviation of dist %f' % bdist.stdDeviationRdmNumber())


### Question 3

'''
(a)
Generate 1,000 random numbers with Binomial distribution with 𝑛=44 and 𝑝=0.64.
'''

binomial = bn.BinomialDistribution(0.64, 44)
binomial.generateBinomialDistribution(1000)
print('Q3. a) Binomial Mean %f' % binomial.meanRdmNumber())
print('       Binomial Std  %f' % binomial.stdDeviationRdmNumber())
'''
(b) (b) Draw the histogram. Compute the probability that the random variable X that has Binomial (44, 0.64) distribution, is at least 40: 𝑃(𝑋≥40).
Use any statistics textbook or online resources for the exact number for the above probability and compare it with your finding and comment.
'''


binomialdist = binomial.getBinomialDistribution()

'''
plt.hist(binomialdist)
plt.title("Histogram")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()
'''
cdf_binfunction = binomial.getCumulativeDistribution(40)
print('Q3. d) Given Binomial distribution %f' % (1-cdf_binfunction))

print('Q3. d) Built in fn %f ' % binom.sf(40, 44, 0.64))

### Question 4
'''
(a) Generate 10,000 Exponentially distributed random numbers with parameter 𝜆=1.5.
'''
lgm1 = rndgen.LGMRandomGenerator(14)
lgm1.generateRdmNumber(10000)
exponentialRdmNum = -1.5* np.log(lgm1.getRdmNumber())
'''
(b) Compute 𝑃(𝑋≥1) and 𝑃(𝑋≥4).
'''
cdf_expfun = wp.getCumulativeDistribution(1, exponentialRdmNum)
print('Q4. b) CDF for given [ P >= 1] exponential is %f' % (1- cdf_expfun) )
cdf_expfun2 = wp.getCumulativeDistribution(4, exponentialRdmNum)
print('Q4. b) CDF for given [ P >= 4] exponential is %f' % (1 - cdf_expfun2))

'''
(c) Compute the empirical mean and the standard deviation of the sequence of 10,000 numbers generated above in part (a).
 Draw the histogram by using the 10,000 numbers of part (a).
'''

print('Q4. c) Mean %f' % exponentialRdmNum.mean())
print('Q4. c) Std dev %f' % exponentialRdmNum.std())

'''
plt.hist(exponentialRdmNum)
plt.title("Histogram")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()
'''

### Question 5

'''
(a) Generate 5,000 Uniformly distributed random numbers on [0,1].
'''
lgm2 = rndgen.LGMRandomGenerator(14)
n = 5000
lgm2.generateRdmNumber(n)

'''
(b) Generate 5,000 Normally distributed random numbers with mean 0 and variance 1, by Box- Muller Method.
'''
U = lgm2.getRdmNumber().reshape(2, n/2)
myln = np.sqrt(-np.log(U[0]))
mycoses = np.cos(2.*np.pi* U[1])
mysines = np.sin(2.*np.pi* U[1])
z1= np.multiply(myln, mycoses) * np.sqrt(2)
z2= np.multiply(myln, mysines) * np.sqrt(2)
z = np.concatenate([z1, z2])

'''
(c) Compute the empirical mean and the standard deviation of the sequence of numbers generated above of part (b).
'''
print(z.mean())
print(z.std())

'''
(d) Now use the Polar-Marsaglia method to do the same as in (b).
'''
v1 = 2* U[0] -1
v2 = 2* U[1] -1
W = v1*v1 + v2*v2
z_p =  array('f', [])
for i in range(0, len(W)):
    if W[i] <=1:
        z_p1 = v1[i] * np.sqrt(-2 * np.log(W[i]) / W[i])
        z_p2 = v2[i] * np.sqrt(-2 * np.log(W[i]) / W[i])
        z_p.append(z_p1)
        z_p.append(z_p2)

print(wp.meanfn(z_p))
print(wp.stdDeviationfn(z_p))