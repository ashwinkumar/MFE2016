import myrandgen as rndgen
import numpy as np
import mywrapper as wp
import matplotlib.pyplot as plt


### Question 4
'''
(a) Generate 10,000 Exponentially distributed random numbers with parameter 𝜆=1.5.
'''
n = 10000
lgm1 = rndgen.LGMRandomGenerator(14)
lgm1.generateRdmNumber(n)
exponentialRdmNum = -1.5* np.log(lgm1.getRdmNumber())
'''
(b) Compute 𝑃(𝑋≥1) and 𝑃(𝑋≥4).
'''
cdf_expfun = wp.getCumulativeDistribution(1, exponentialRdmNum)
print('Q4. b) CDF for given [ P >= 1] exponential is %f' % (1- cdf_expfun) )
cdf_expfun2 = wp.getCumulativeDistribution(4, exponentialRdmNum)
print('       CDF for given [ P >= 4] exponential is %f' % (1 - cdf_expfun2))

'''
(c) Compute the empirical mean and the standard deviation of the sequence of 10,000 numbers generated above in part (a).
 Draw the histogram by using the 10,000 numbers of part (a).
'''

print('Q4. c) Mean %f' % exponentialRdmNum.mean())
print('       Std dev %f' % exponentialRdmNum.std())

bins = np.linspace(0, 10, 100)
plt.hist(exponentialRdmNum, bins, alpha=0.5)
plt.title("Histogram of exponential distribution")
plt.xlabel("Random Number")
plt.ylabel("Frequency")
plt.show()