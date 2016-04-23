import sys
sys.path.append("/Users/akumar/Python/com/ashwin/computationalmethodsinfinance")
import numpy as np
import matplotlib.pyplot as plt
from array import *


import myrandgen as rndgen



## Question 1

'''
(a) Using LGM method generate 10,000 Uniformly distributed random numbers on [0,1] and
compute the empirical mean and the standard deviation of the sequence.
'''
n = 10000
lgm = rndgen.LGMRandomGenerator(14)
lgm.generateRdmNumber(n)
myrdn = lgm.getRdmNumber()

bins = np.linspace(0, 1, 20)


plt.hist(myrdn, bins, alpha=0.5)
plt.title("Histogram of Uniform Random Number (0,1)")
plt.xlabel("Random Number")
plt.ylabel("Frequency")
plt.draw()

print('Q1. a) LGM Mean %f' % lgm.meanRdmNumber())
print('       LGM Std Dev %f' % lgm.stdDeviationRdmNumber())

'''
(b) Now use built-in functions of whatever software you are using to do the same thing as in (a).
'''
x_array = np.random.uniform (0, 1, n)
print('Q1. b) Built in Fn Mean %f' % x_array.mean())
print('       Built in Fn Std Dev %f' % x_array.std())


'''
(c) Compare your findings in (a) and (b) and comment (be short but precise).
'''
print(' The mean and standard deviation values of computed uniform distribution is close to built-in distribution. Also it is approximately equal to actual values of (0.5,0.288)')
plt.show()
