import myrandgen as rndgen
import mynormal as nl
import mywrapper as wp
import timeit
import numpy as np
import matplotlib.pyplot as plt

### Question 5
np.seterr(all='raise')

'''
(a) Generate 5,000 Uniformly distributed random numbers on [0,1].
'''
n = 5000
lgm2 = rndgen.LGMRandomGenerator(14)
lgm2.generateRdmNumber(n)

'''
(b) Generate 5,000 Normally distributed random numbers with mean 0 and variance 1, by Box- Muller Method.
'''
boxmullernormal = nl.NormalDistribution("Box-Muller")
boxmullernormal.generateNormalDistribution(n)
nrm1 = list(boxmullernormal.getRdmNumbers())

bins = np.linspace(-10, 10, 1000)
plt.hist(nrm1, bins, alpha=0.5)
plt.title("Histogram of Normal distribution : Box-Muller Method")
plt.xlabel("Random Number")
plt.ylabel("Frequency")
#plt.show()

'''
(c) Compute the empirical mean and the standard deviation of the sequence of numbers generated above of part (b).
'''
#print('Q5.  b) Mean %f' % nrm1.mean())
#print('     c) Std Deviation %f' % nrm1.std())
print('Q5.  b) Mean %f' % wp.meanfn(nrm1))
print('     c) Std Deviation %f' % wp.stdDeviationfn(nrm1))


'''
(d) Now use the Polar-Marsaglia method to do the same as in (b).
'''

polarmarsaglia = nl.NormalDistribution("Polar-Marsaglia")
polarmarsaglia.generateNormalDistribution(n)
nrm2 = list(polarmarsaglia.getRdmNumbers())

bins = np.linspace(-10, 10, 1000)
plt.hist(nrm2, bins, alpha=0.5)
plt.title("Histogram of Normal distribution : Polar-Marsaglia Method")
plt.xlabel("Random Number")
plt.ylabel("Frequency")
#plt.show()
'''
(e) Compute the empirical mean and the standard deviation of the sequence of numbers generated above of part
'''
print('Q5.  b) Mean %f' % wp.meanfn(nrm2))
print('     c) Std Deviation %f' % wp.stdDeviationfn(nrm2))


'''
(f) Now compare the efficiencies of the two above-algorithms, by comparing the execution times to generate 5,000 normally distributed random numbers by the two methods. Which one is more efficient?
'''

#N = 1000000000
N = 10000
start_time1 = timeit.default_timer()
bmnormal = nl.NormalDistribution("Box-Muller")
bmnormal.generateNormalDistribution(N)
n1 = list(bmnormal.getRdmNumbers())
elapsed1 = timeit.default_timer() - start_time1


start_time2 = timeit.default_timer()

pnormal = nl.NormalDistribution("Polar-Marsaglia")
pnormal.generateNormalDistribution(N)
n2 = list(pnormal.getRdmNumbers())

elapsed2 = timeit.default_timer() - start_time2

start_time3 = timeit.default_timer()
uniformrd = rndgen.LGMRandomGenerator(14)
u1 = list(uniformrd.generateRdmNumberByGenerator(N))
elapsed3 = timeit.default_timer() - start_time3

print('     f) Run time for Box-Muller Normal Distribution %f' % elapsed1)
print('        Run time for Polar-Marsaglia Normal Distribution %f' % elapsed2)
print('        Polar-Marsaglia method is faster when n= %d' % N)
print('        Run time for Uniform Random Distribution %f' % elapsed3)
