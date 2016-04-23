__author__= 'ashwin'

import mybivariatenormal as bv
import math
import  statistics as st
## Question 2

'''
(a) Evaluate the following expected values by using Monte Carlo simulation:  Ε [max (0,(𝑋3+sin(𝑌)+𝑋2𝑌))]
where X and Y have 𝑁(0,1) distribution and a correlation of 𝜌=0.6.
'''

def myMC(n, X1, Y1):
    num = 0
    while num < n:
        xsq = X1[num]*X1[num]
        yield max(0, xsq * X1[num] + Y1[num]*xsq + math.sin(Y1[num]))
        num = num +1

n = 10000
p = 0.6


binormal = bv.BiVariateNormalDistribution(p)
binormal.generateDistribution(2*n)
Z = list(binormal.getRdmNumbers())

X1 = [item[0] for item in Z]
Y1 = [item[1] for item in Z]
sim = list(myMC(n, X1, Y1))
res = st.mean(sim)
var = st.variance(sim)/len(sim)
print('Q2. The value of expression by MC Simulaiton =  %f with var = %f' % (res,var) )


