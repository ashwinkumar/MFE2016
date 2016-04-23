import numpy as np
import mybernoulli as bl
import matplotlib.pyplot as plt

## Question 2

'''
(a)
Generate 10,000 random numbers with the following distribution:
'''
n = 10000
p = np.array((0.3, 0.35, 0.2, 0.15 ))
X = np.array((-1, 0, 1, 2))
bernouli = bl.BernoulliDistribution(p, X)
bernouli.generateBernoulliDistribution(n)

'''
(b)
Draw the histogram and compute the empirical mean and the standard deviation of the sequence of 10,000 numbers generated above in part (a).
'''

bdist = bernouli.getBernoulDistNumbers()
plt.hist(bdist, alpha=0.5)
plt.title("Histogram of Distribution with p = (0.3,0.35,0.2,0.15)")
plt.xlabel("Random Number")
plt.ylabel("Frequency")
plt.show()

print('Q2. b) Mean of given dist %f' % bernouli.meanRdmNumber())
print('Q2. b) Std deviation of dist %f' % bernouli.stdDeviationRdmNumber())
