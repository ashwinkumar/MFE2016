import TrinomialModel as tm
import numpy as np
import matplotlib.pyplot as plt
import math
import BlackScholesOptionPrice as bp
#### Question 5 - Trinomial Model
def EuropeanCallPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(St-X,array_zeros)

ax = plt.subplot(1, 1, 1)

r = 0.05
sigma = 0.24
S0 = 32
X = 30
n = np.array([10,15,20,40,70,80,100,200,500])
l_n = len(n)
call_price = np.zeros(l_n)
t = 0.5
dt = t/n
trimodel = tm.TrinomialModel(S0,r,t)
'''
(a) Use the trinomial method applied to the stock price-process (𝑆𝑡) in which
'''
u,d,pu,pm,pd = tm.getBinomialParameters1(r,sigma,dt)

for i in range(0,l_n):
    trimodel.n= n[i]
    if trimodel.verify(u[i],d[i],pu[i],pm[i],pd[i]):
        call_price[i] = trimodel.getOptionPrice(EuropeanCallPayOff,X,"E")
    trimodel.resetParams()
plt.plot(n,call_price,label='Method 1 [St]')


'''
(b) Use the trinomial method applied to the Log-stock price-process (𝑋𝑡) in which
'''
trimodel1 = tm.TrinomialModel(math.log(S0),r,t)

u,d,pu,pm,pd = tm.getBinomialParameters2(r,sigma,dt)

for i in range(0,l_n):
    trimodel1.n= n[i]
    if trimodel1.verifylog(u[i],d[i],pu[i],pm[i],pd[i]):
        call_price[i] = trimodel1.getOptionPricebyLog(EuropeanCallPayOff,X)
    trimodel1.resetParams()

plt.plot(n,call_price,label='Method 2 [Ln St]')


# Compute BS option price for comparing
black_scholes_price = bp.EuropeanCallPrice(S0,X,r,sigma,t)
b= black_scholes_price* np.ones(l_n)
plt.plot(n,b,"k--",label='BS Price')


plt.legend(loc="upper left", bbox_to_anchor=[0, 1], ncol=2, shadow=True, title="Legend", fancybox=True)
ax.get_legend().get_title().set_color("red")

plt.show()
