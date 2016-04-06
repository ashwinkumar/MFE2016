import mybinomial as bn
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binom

### Question 3

'''
(a)
Generate 1,000 random numbers with Binomial distribution with 𝑛=44 and 𝑝=0.64.
'''
n = 1000
binomial = bn.BinomialDistribution(0.64, 44)
binomial.generateBinomialDistribution(n)
print('Q3. a) Binomial Mean %f' % binomial.meanRdmNumber())
print('       Binomial Std  %f' % binomial.stdDeviationRdmNumber())
'''
(b) (b) Draw the histogram. Compute the probability that the random variable X that has Binomial (44, 0.64) distribution, is at least 40: 𝑃(𝑋≥40).
Use any statistics textbook or online resources for the exact number for the above probability and compare it with your finding and comment.
'''


binomialdist = binomial.getBinomialDistribution()
bins = np.linspace(0, 44, 45)
plt.hist(binomialdist, bins, alpha=0.5)
plt.title("Histogram of Binomial distribution(0.64,44)")
plt.xlabel("Random Number")
plt.ylabel("Frequency")
plt.show()

cdf_binfunction = binomial.getCumulativeDistribution(40)
print('Q3. b) My Binomial P(X>= 40) %f' % (1-cdf_binfunction))

print('       Built in Binomial P(X>=40) %f ' % binom.sf(40, 44, 0.64))
