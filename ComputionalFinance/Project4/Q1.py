import BinomialModel as bm
import math
import numpy as np
import timeit
import matplotlib.pyplot as plt
import BlackScholesOptionPrice as bp
###Question 1

def EuropeanCallPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(St-X,array_zeros)

ax = plt.subplot(1, 1, 1)
plt.xlabel("No of intervals")
plt.ylabel("Call Option Price")
t= 0.5
r = 0.05
sigma = 0.24
S0 = 32
X = 30
n_periods = np.array([10,20,40,80,100,200,500])
l_n = len(n_periods)
call_price = np.zeros(l_n)
dt = t/n_periods #[t/i for i in n_periods]

'''
Q1. a)   u=1/d
'''
u,d,p = bm.getBinomialParameters(r,sigma,dt,1)
bin_model = bm.BinomialModel(S0,r, t)

start_time1 = timeit.default_timer()

for i in range(0,l_n):
    bin_model.n= n_periods[i]
    if bin_model.verify(u[i], d[i],p[i]):
        call_price[i] = bin_model.getOptionPrice(EuropeanCallPayOff,X,"E")
    bin_model.resetParams()


elapsed1 = timeit.default_timer() - start_time1


print('Q1. a) Method 1: u=1/d    :  ****[%f sec]' % elapsed1)
print(call_price)

plt.plot(n_periods,call_price,label='Method 1 [u=1/d]')

'''
plt.plot(n_periods, call_price)
plt.ylabel("Call Price")
plt.xlabel("No of intervals")
plt.title("European Call Price Vs No. of intervals  [u=1/d]")
plt.show()
'''

'''
Q1. b)  p=1/2
'''
u,d,p = bm.getBinomialParameters(r,sigma,dt,2)
bin_model.resetParams()

start_time2 = timeit.default_timer()

for i in range(0,l_n):
    bin_model.n= n_periods[i]
    if bin_model.verify(u[i], d[i],p[i]):
        call_price[i] = bin_model.getOptionPrice(EuropeanCallPayOff,X,"E")
    bin_model.resetParams()

elapsed2 = timeit.default_timer() - start_time2


print('Q1. b) Method 2: p=1/2    :  ****[%f sec]' % elapsed2)
print(call_price)
plt.plot(n_periods,call_price, label='Method 2 [p=1/2]')
'''
plt.plot(n_periods, call_price)
plt.ylabel("Call Price")
plt.xlabel("No of intervals")
plt.title("European Call Price Vs No. of intervals  [p=1/2]")
plt.show()
'''

'''
Q1.  c) p=1/2, Jarrow-Rudd Model
'''

u,d,p = bm.getBinomialParameters(r,sigma,dt,3)
bin_model.resetParams()

start_time3 = timeit.default_timer()
for i in range(0,l_n):
    bin_model.n= n_periods[i]
    if bin_model.verify(u[i], d[i],p[i]):
        call_price[i] = bin_model.getOptionPrice(EuropeanCallPayOff,X,"E")
    bin_model.resetParams()

elapsed3 = timeit.default_timer() - start_time3


print('Q1. c) Method 3: J-R Model    :  ****[%f sec]' % elapsed3)
print(call_price)
plt.plot(n_periods,call_price,label='Method 3 [p=1/2, J-R Model]')

'''
plt.plot(n_periods, call_price)
plt.ylabel("Call Price")
plt.xlabel("No of intervals")
plt.title("European Call Price Vs No. of intervals  [J-R Model]")
plt.show()
'''

u,d,p = bm.getBinomialParameters(r,sigma,dt,4)
bin_model.resetParams()

start_time4 = timeit.default_timer()
for i in range(0,l_n):
    bin_model.n= n_periods[i]
    if bin_model.verify(u[i], d[i],p[i]):
        call_price[i] = bin_model.getOptionPrice(EuropeanCallPayOff,X,"E")
    bin_model.resetParams()

elapsed4 = timeit.default_timer() - start_time4
print('Q1. d) Method 4: CRR Model    :  ****[%f sec]' % elapsed4)
print(call_price)

plt.plot(n_periods,call_price, label='Method 4 [p=1/2, CRR Model]')

'''
plt.plot(n_periods, call_price)
plt.ylabel("Call Price")
plt.xlabel("No of intervals")
plt.title("European Call Price Vs No. of intervals  [CRR Model]")
plt.show()
'''

plt.legend(loc="upper left", bbox_to_anchor=[0, 1],
           ncol=2, shadow=True, title="Legend", fancybox=True)
ax.get_legend().get_title().set_color("red")

black_scholes = bp.EuropeanCallPrice(S0,X,r,sigma,t)
b = black_scholes*np.ones(l_n)
plt.plot(n_periods,b,"k--",label='BS Price')

plt.show()
