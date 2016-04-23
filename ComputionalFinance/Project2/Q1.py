__author__= 'ashwin'

import mybivariatenormal as bv
## Question 1
n = 10000
a = -0.7


binormal = bv.BiVariateNormalDistribution(a)
binormal.generateDistribution(n)
binormal.getRdmNumbers()
p = binormal.getCorr()

print('Q1. The value of p from simulation is %f' % p)





