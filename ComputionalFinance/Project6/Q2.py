from ExoticOption import DefaultOption as do
import numpy as np
import matplotlib.pyplot as plt
from  tabulate import  tabulate

#### Question 2- Default Option Pricing using Jump Diffusion process

V0 = 20000
L0 = 22000
mu = -0.1
sigma = 0.2
gamma = -0.4

l1_start = 0.05
l1_end  =0.4
l1_step = 0.05
lambda1 = np.arange(l1_start,l1_end+l1_step,l1_step)
l_l1 = len(lambda1)

l2_start = 0
l2_end = 0.8
l2_step =0.1
lambda2 = np.arange(l2_start,l2_end+ l2_step, l2_step)
l_l2 = len(lambda2)

alpha = 0.7
ephsilon = 0.95

T_start = 3
T_end = 8
T = np.arange(T_start,T_end+1,1)
l_T = len(T)

r0 = 0.02
delta = 0.25
ndiv = 12 * T  # 12 months per year
nsims = 10000
defaultOptionPrice = np.zeros([l_T,l_l1,l_l2])
defaultProb = np.zeros([l_T,l_l1,l_l2])
defaultExerciseTime = np.zeros([l_T,l_l1,l_l2])

for i in range(l_T):
    def_option= do(T[i],ndiv[i],nsims)
    for j in range(l_l1):
        def_option.initCollateral(V0, mu, sigma, gamma, lambda1[j])
        def_option.simulateCollateral()
        for k in range(l_l2):
            def_option.initLoan(L0, r0, delta, lambda2[k])
            def_option.setThreshold(alpha, ephsilon)
            def_option.simulateLoanValue()
            defaultOptionPrice[i,j,k],_ = def_option.getDefaultOptionPrice(r0)
            defaultProb[i,j,k] = def_option.getDefaultProb()
            defaultExerciseTime[i,j,k] = def_option.getExpectedExerciseTime()
    del def_option

print('Q1. a) Default Option Prices with T = 5 , lambda1 = 0.2, lambda2 = 0.4  = %f' %defaultOptionPrice[2,5,5])
print('Q1. b) Default Probability  with T = 5 , lambda1 = 0.2, lambda2 = 0.4  = %f' %defaultProb[2,5,5])
print('Q1. c) Expected Exercise Time  with T = 5 , lambda1 = 0.2, lambda2 = 0.4  = %f' %defaultExerciseTime[2,5,5])


### Default Option Price
## Fixing lambda1 = 0.2
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)
ax = plt.subplot(211)
plt.title('Default Option Price for $\lambda1$= 0.2')
plt.xlabel("Time in Years")
plt.ylabel('Option Price')
for i in range(l_l2):
    ax.plot(T,defaultOptionPrice[:,5,i],label = "$\lambda2 = %0.2f$"%lambda2[i])



box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.show()
#plt.savefig("Proj6_2a.png")


## Fixing lambda2 = 0.4
#fig = plt.figure()
ax = plt.subplot(212)

plt.title('Default Option Price for $\lambda2$ = 0.4')
plt.xlabel("Time in Years")
plt.ylabel('Option Price')
for i in range(l_l1):
    ax.plot(T,defaultOptionPrice[:,i,5],label = "$\lambda1 = %0.2f$"%lambda1[i])



box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.show()
plt.savefig("Proj6_2a.png")




###  Default Probability
## Fixing lambda1 = 0.2
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)
ax = plt.subplot(211)
plt.title('Default Probability for $\lambda1$= 0.2')
plt.xlabel("Time in Years")
plt.ylabel('Probability')
for i in range(l_l2):
    ax.plot(T,defaultProb[:,5,i],label = "$\lambda2 = %0.2f$"%lambda2[i])



box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.show()
#plt.savefig("Proj6_2a.png")


## Fixing lambda2 = 0.4
#fig = plt.figure()
ax = plt.subplot(212)

plt.title('Default Probability for $\lambda2$ = 0.4')
plt.xlabel("Time in Years")
plt.ylabel('Probability')
for i in range(l_l1):
    ax.plot(T,defaultProb[:,i,5],label = "$\lambda1 = %0.2f$"%lambda1[i])



box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.show()
plt.savefig("Proj6_2b.png")


###  Expected Time in Years
## Fixing lambda1 = 0.2
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)
ax = plt.subplot(211)
plt.title('Expected Exercise Time for $\lambda1$= 0.2')
plt.xlabel("Time in Years")
plt.ylabel('Expected Exercise Time ')
for i in range(l_l2):
    ax.plot(T,defaultExerciseTime[:,5,i],label = "$\lambda2 = %0.2f$"%lambda2[i])



box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.show()
#plt.savefig("Proj6_2a.png")


## Fixing lambda2 = 0.4
#fig = plt.figure()
ax = plt.subplot(212)

plt.title('Expected Exercise Time for $\lambda2$ = 0.4')
plt.xlabel("Time in Years")
plt.ylabel('Expected Exercise Time ')
for i in range(l_l1):
    ax.plot(T,defaultExerciseTime[:,i,5],label = "$\lambda1 = %0.2f$"%lambda1[i])



box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.show()
plt.savefig("Proj6_2c.png")