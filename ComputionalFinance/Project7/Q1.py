import FiniteDifferenceMethod as FDM
import BlackScholesOptionPrice as bs
import math
import numpy as np
import matplotlib.pyplot as plt
from  tabulate import  tabulate

#### Question 1

def EuropeanCallPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(St-X,array_zeros)

def EuropeanPutPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(X- St,array_zeros)


def getStockIndices(S,S1,S2):
    ind = []
    iS = S1
    for i in range(len(S)):
        if S[i] > iS:
            iS = iS+1
            ind.append(i)
    return np.asarray(ind)

S0 = 10
Smin = 4
Smax = 16
X0 = math.log(S0)
Xmin = math.log(Smin)
Xmax = math.log(Smax)

sigma = 0.20
r = 0.04
dT = 0.002
T = 0.5
K = 10
dX = np.array([1, math.sqrt(3),2])* sigma*math.sqrt(dT)
l_dx = len(dX)

for i in range(l_dx):
    filename = "Proj7_1a_dx_%d"%i
    title = "$\Delta X = %0.4f$"%dX[i]
    fdm_class = FDM.LogFiniteDifferenceMethod(X0,r,sigma,T,dT)
    S,indS0 = fdm_class.initLogStock(Xmin, Xmax, dX[i])
    fdm_class.setProbabilities("EFD")
    optionPrice_EFM = fdm_class.getOptionPrice(K,EuropeanPutPayOff,"EFD")
    fdm_class.setProbabilities("IFD")
    optionPrice_IFM = fdm_class.getOptionPrice(K,EuropeanPutPayOff,"IFD")
    fdm_class.setProbabilities("CNFD")
    optionPrice_CNFD = fdm_class.getOptionPrice(K,EuropeanPutPayOff,"CNFD")

    blsOptionPrice = bs.EuropeanPutPrice(S,K,r,sigma,T)
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.xlabel("Stock Price")
    plt.ylabel("Euopean Put Option Price")
    plt.title(title)
    ax.plot(S,blsOptionPrice,label= "BS")
    ax.plot(S, optionPrice_EFM, label= "EFM")
    ax.plot(S,optionPrice_IFM,"k--", label = "IFM")
    ax.plot(S,optionPrice_CNFD,"k:", label = "CNFM")
    plt.xlim(4,16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


    #ax.legend(loc="upper left", bbox_to_anchor=[0, 1], ncol=4, shadow=True, title="Legend", fancybox=True,prop={'size':8})
    ax.legend(loc="center left", bbox_to_anchor=[1, 0.5], ncol=1, shadow=True, title="Legend", fancybox=True)
    plt.savefig(filename+".png")
    plt.clf()
    print("Q1.a Put Option Price for S= %f, dX = %f"%(S[indS0],dX[i]))
    print("     Price = %f [Black Scholes]" % (blsOptionPrice[indS0]))
    print("     Price = %f [Explicit Method]"%(optionPrice_EFM[indS0]))
    print("     Price = %f [Implicit Method]"%(optionPrice_IFM[indS0]))
    print("     Price = %f [Crank-Nicolson Method]"%(optionPrice_CNFD[indS0]))

    ind = getStockIndices(S,4,16)
    S = S[ind]
    blsOptionPrice = blsOptionPrice[ind]
    optionPrice_EFM = optionPrice_EFM[ind]
    optionPrice_IFM = optionPrice_IFM[ind]
    optionPrice_CNFD =optionPrice_CNFD[ind]
    result = np.zeros([len(S),5,2])
    result[:,0,0] = S
    result[:,1,0] = blsOptionPrice
    result[:,2,0] = optionPrice_EFM
    result[:,2,1]=  ((optionPrice_EFM/ blsOptionPrice) -1)*100
    result[:,3,0]= optionPrice_IFM
    result[:,3,1]= ((optionPrice_IFM/blsOptionPrice) -1)*100
    result[:,4,0] = optionPrice_CNFD
    result[:,4,1] = ((optionPrice_CNFD/blsOptionPrice)-1)*100

    print(tabulate(result, headers=["S,na","Black-Scholes, % Error ","Explicit Method, % Error","Implicit Method, %Error","Crank-Nicolson Method, %Error"]))
    del fdm_class

