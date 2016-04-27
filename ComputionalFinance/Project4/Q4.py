import numpy as np
import timeit
import BinomialModel as bm
import matplotlib.pyplot as plt

#### Question 4

def EuropeanPutPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(X-St,array_zeros)

r = 0.05
sigma = 0.3
X= 100
start = 80
end = 120
step = 4
S0 = np.arange(start,end+step,step)
t=1
n=100

l_it = len(S0)
put_price_european =  np.zeros(l_it)
put_price_american = np.zeros(l_it)

dt = t/n
u,d,p = bm.getBinomialParameters(r,sigma,dt,4)
start_time1 = timeit.default_timer()

for i in range(0,l_it):
    bin_model = bm.BinomialModel(S0[i], r, t)
    bin_model.n= n
    if bin_model.verify(u, d,p):
        put_price_american[i] = bin_model.getOptionPrice(EuropeanPutPayOff,X,"A")
        put_price_european[i] = bin_model.getOptionPrice(EuropeanPutPayOff,X,"E")

    del bin_model

elapsed1 = timeit.default_timer() - start_time1
plt.plot(S0,put_price_american,label="American")
plt.plot(S0,put_price_european, label="European")
plt.ylabel("Put Option Price")
plt.xlabel("Stock Price")
plt.title("Put Price Vs Stock  [CRR Method]")
plt.legend(loc="upper left", bbox_to_anchor=[0, 1], ncol=2, shadow=True, title="Legend", fancybox=True)
plt.show()

plt.show()