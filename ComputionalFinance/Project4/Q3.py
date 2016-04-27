import numpy as np
import BinomialModel as bm
import timeit
import matplotlib.pyplot as plt

#### Question 3
def EuropeanCallPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(St-X,array_zeros)

start_time = timeit.default_timer()
ephsilon_perc = 0.01 # Use in 1 % percentage
epshilon = 0.01

fig = plt.figure()
fig.suptitle('Greeks for European Call [CRR Model for Stock Price]')

'''
(i) Delta of the call option as a function of 𝑆0, for 𝑆0ranging from $20 to $80, in increments of $2.
'''
start = 20
end =80
step =2
S0 = np.arange(start,end+step, step)
S0_ephsilon_perc = S0*(1 + ephsilon_perc)
l_it = len(S0)

#delta_st_1 = np.zeros(l_it)  # Using the tree nodes
#call_1 = np.zeros(l_it)
#call_1ep = np.zeros(l_it)
delta_st = np.zeros(l_it)
X = 50
r = 0.03
sigma = 0.2
t = 0.3846
mu = 0.14

n = 100
dt = t/n
u,d,p = bm.getBinomialParameters(r,sigma,dt,4)

start_time1 = timeit.default_timer()

for i in range(0,l_it):
    bin_model = bm.BinomialModel(S0[i], r, t)
    bin_model1 = bm.BinomialModel(S0_ephsilon_perc[i],r,t)
    bin_model.n= n
    bin_model1.n=n
    if bin_model.verify(u, d,p) and bin_model1.verify(u, d, p):
        #delta_st_1[i] = bin_model.getDelta(EuropeanCallPayOff,X)
        call_1 = bin_model.getOptionPrice(EuropeanCallPayOff,X,"E")
        call_1ep= bin_model1.getOptionPrice(EuropeanCallPayOff, X,"E")
        delta_st[i] = (call_1ep- call_1)/(ephsilon_perc*S0[i])
    del bin_model
    del bin_model1

ax1 = fig.add_subplot(321)
ax1.plot(S0, delta_st, "k--", color = "red")
ax1.set_ylabel("Delta of Call Price")
ax1.set_xlabel("Stock Price")
'''
(ii) Delta of the call option, as a function of T (time to expiration), from 0 to 0.3846 in increments of 0.01.
'''
start = 0
end =0.3846
step =0.01
S0 = 49
S0_ephsilon_perc = S0*(1+ephsilon_perc)
t = np.arange(start+step,end+step,step)
l_it = len(t)
delta_t = np.zeros(l_it)

dt = t/n
u,d,p = bm.getBinomialParameters(r,sigma,dt,4)
start_time2 = timeit.default_timer()

for i in range(0,l_it):
    bin_model = bm.BinomialModel(S0, r, t[i])
    bin_model1 = bm.BinomialModel(S0_ephsilon_perc, r, t[i])
    bin_model.n= n
    bin_model1.n = n
    if bin_model.verify(u[i], d[i],p[i]) and bin_model1.verify(u[i], d[i], p[i]):
        call_1 = bin_model.getOptionPrice(EuropeanCallPayOff,X,"E")
        call_1ep = bin_model1.getOptionPrice(EuropeanCallPayOff, X,"E")
        delta_t[i] = (call_1ep- call_1)/(ephsilon_perc*S0)
    del bin_model
    del bin_model1
ax2 = fig.add_subplot(322)
ax2.set_ylabel("Delta of Call Price")
ax2.set_xlabel("Time to expiration")
ax2.plot(t,delta_t,"k--",color= "blue")


'''
(iii) Theta of the call option, as a function of 𝑆0, for 𝑆0ranging from $20 to $80 in increments of $2.
'''

start = 20
end =80
step =2
S0 = np.arange(start,end+step, step)
l_it = len(S0)
theta_st = np.zeros(l_it)
t = 0.3846
t_epsilon = t*(1-ephsilon_perc)

dt = t/n
dt_e = t_epsilon/n
u,d,p = bm.getBinomialParameters(r,sigma,dt,4)
u_e,d_e,p_e = bm.getBinomialParameters(r,sigma,dt_e,4)


for i in range(0,l_it):
    bin_model = bm.BinomialModel(S0[i], r, t)
    bin_model_e = bm.BinomialModel(S0[i], r, t_epsilon)
    bin_model.n= n
    bin_model_e.n= n
    if bin_model.verify(u, d,p) and bin_model_e.verify(u_e,d_e,p_e):
        #theta_st[i] = bin_model.getTheta(EuropeanCallPayOff,X)
        call_1 = bin_model.getOptionPrice(EuropeanCallPayOff, X, "E")
        call_1ep = bin_model_e.getOptionPrice(EuropeanCallPayOff, X, "E")
        theta_st[i] = (call_1ep - call_1) / (ephsilon_perc * t)  #ephsilon_perc in percentage
    del bin_model
    del bin_model_e

#plt.plot(S0, theta_st,label = "Theta vs Stock Price")
ax3 = fig.add_subplot(323)
ax3.set_ylabel("Theta of Call Price")
ax3.set_xlabel("Stock Price")
ax3.plot(S0,theta_st, "k--",color = "black")

'''
(iv) Gamma of the call option, as a function of 𝑆0, for 𝑆0ranging from $20 to $80 in increments of $2.
'''

start = 20
end =80
step =2
ephsilon_perc = 0.25
S0 = np.arange(start,end+step, step)
S0_ephsilon_perc1 = S0*(1 + ephsilon_perc)
S0_ephsilon_perc2 = S0*(1 + 2*ephsilon_perc)
l_it = len(S0)
gamma = np.zeros(l_it)

u,d,p = bm.getBinomialParameters(r,sigma,dt,4)

start_time4 = timeit.default_timer()

for i in range(0,l_it):
    bin_model = bm.BinomialModel(S0[i], r, t)
    bin_model1 = bm.BinomialModel(S0_ephsilon_perc1[i],r,t)
    bin_model2 = bm.BinomialModel(S0_ephsilon_perc2[i], r, t)
    bin_model.n = n
    bin_model1.n = n
    bin_model2.n = n
    if bin_model.verify(u, d,p) and bin_model1.verify(u, d, p) and bin_model2.verify(u,d,p):
        call_1 = bin_model.getOptionPrice(EuropeanCallPayOff,X,"E")
        call_1ep= bin_model1.getOptionPrice(EuropeanCallPayOff, X,"E")
        call_2ep= bin_model2.getOptionPrice(EuropeanCallPayOff, X,"E")
        gamma[i] = (call_2ep- 2*call_1ep + call_1 )/(ephsilon_perc*S0[i])**2 # epsilon in percentage
    del bin_model
    del bin_model1
    del bin_model2

ax4 = fig.add_subplot(324)
ax4.set_ylabel("Gamma of Call Price")
ax4.set_xlabel("Stock Price")
ax4.plot(S0, gamma ,"k--", color = "green")

'''
(v) Vega of the call option, as a function of 𝑆0, for 𝑆0ranging from $20 to $80 in increments of $2.
'''

vega = np.zeros(l_it)
sigma_epsilon = sigma*(1+ephsilon_perc)
u,d,p = bm.getBinomialParameters(r,sigma,dt,4)
u_1,d_1,p_1 = bm.getBinomialParameters(r,sigma_epsilon,dt,4)

for i in range(0,l_it):
    bin_model = bm.BinomialModel(S0[i], r, t)
    bin_model1 = bm.BinomialModel(S0[i],r,t)
    bin_model.n = n
    bin_model1.n = n
    if bin_model.verify(u, d,p) and bin_model1.verify(u_1, d_1, p_1):
        call_1 = bin_model.getOptionPrice(EuropeanCallPayOff,X,"E")
        call_1ep= bin_model1.getOptionPrice(EuropeanCallPayOff, X,"E")
        vega[i] = (call_1ep- call_1 )/(ephsilon_perc*sigma) # epsilon is in % of sigma
    del bin_model
    del bin_model1


ax5 = fig.add_subplot(325)
ax5.set_ylabel("Vega of Call Price")
ax5.set_xlabel("Stock Price")
ax5.plot(S0, vega,"k--", color = "violet")

'''
(1v) Rho of the call option, as a function of 𝑆0, for 𝑆0ranging from $20 to $80 in increments of $2.
'''

rho = np.zeros(l_it)
r_epsilon = r*(1+ephsilon_perc)
u,d,p = bm.getBinomialParameters(r,sigma,dt,4)
u_1,d_1,p_1 = bm.getBinomialParameters(r_epsilon,sigma,dt,4)

for i in range(0,l_it):
    bin_model = bm.BinomialModel(S0[i], r, t)
    bin_model1 = bm.BinomialModel(S0[i],r_epsilon,t)
    bin_model.n = n
    bin_model1.n = n
    if bin_model.verify(u, d,p) and bin_model1.verify(u_1, d_1, p_1):
        call_1 = bin_model.getOptionPrice(EuropeanCallPayOff,X,"E")
        call_1ep= bin_model1.getOptionPrice(EuropeanCallPayOff, X,"E")
        rho[i] = (call_1ep- call_1 )/(ephsilon_perc*r) # epsilon is in % of sigma
    del bin_model
    del bin_model1

#plt.plot(S0, rho,label = "Vega vs Stock Price")

ax6 = fig.add_subplot(326)
ax6.set_ylabel("Rho of Call Price")
ax6.set_xlabel("Stock Price")
ax6.plot(S0,rho,"k--", color = "brown")



elapsed = timeit.default_timer() - start_time

print('Q3. 6 Graphs plotted.      :  ****[%f sec]' % elapsed)
plt.show()
