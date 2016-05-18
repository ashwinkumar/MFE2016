import numpy as np
import BlackScholesOptionPrice as bs
import math
import matplotlib.pyplot as plt
import FiniteDifferenceMethod as FDM
from  tabulate import  tabulate



def CallPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(St-X,array_zeros)

def PutPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(X- St,array_zeros)


S0 = 10
sigma = 0.2
r = 0.04
K = 10
T = 0.5
dT = 0.002
Smin = 0
Smax = 16
dS = np.array([0.5,1,1.5])
l_dS = len(dS)

for i in range(l_dS):
    filename1 = "Proj7_2_Putdx_%d"%i
    filename2 = "Proj7_2_Calldx_%d" % i
    title = "$\Delta X = %0.4f$"%dS[i]
    fdm_class = FDM.FiniteDifferenceMethod(S0,r,sigma,T,dT)
    S,indS0 =fdm_class.initStock(Smin,Smax,dS[i])
    putoptionPrice_EFM = fdm_class.getOptionPrice(K, PutPayOff, alpha=1,type= "P")
    putoptionPrice_IFM = fdm_class.getOptionPrice(K, PutPayOff, alpha=0,type="P")
    putoptionPrice_CNFM = fdm_class.getOptionPrice(K, PutPayOff, alpha=0.5, type ="P")

    calloptionPrice_EFM = fdm_class.getOptionPrice(K, CallPayOff, alpha=1,type = "C")
    calloptionPrice_IFM = fdm_class.getOptionPrice(K, CallPayOff, alpha=0,type = "C")
    calloptionPrice_CNFM = fdm_class.getOptionPrice(K, CallPayOff, alpha=0.5,type = "C")

    #blsPutOptionPrice = bs.EuropeanPutPrice(S, K, r, sigma, T)
    #blsCallOptionPrice = bs.EuropeanCallPrice(S, K, r, sigma, T)
    if( dS[i]==1):
        fig = plt.figure()
        ax = plt.subplot(111)

        plt.xlabel("Stock Price")
        plt.ylabel("American Put Option Price")
        plt.xlim(4,16)
        plt.ylim(0,7)
        plt.title(title)
        #plt.plot(S,blsPutOptionPrice,label="Black Scholes European")
        ax.plot(S,putoptionPrice_EFM,label = "EFM")
        ax.plot(S,putoptionPrice_IFM,label = "IFM")
        ax.plot(S, putoptionPrice_CNFM,"k--", label="CNFM")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc="center left", bbox_to_anchor=[1, 0.5], ncol=1, shadow=True, title="Legend", fancybox=True)
        #plt.legend(loc="upper left", bbox_to_anchor=[0, 1], ncol=2, shadow=True, title="Legend", fancybox=True)
        plt.savefig(filename1 + ".png")
        plt.clf()

        fig = plt.figure()
        ax = plt.subplot(111)
        plt.xlabel("Stock Price")
        plt.ylabel("American Call Option Price")
        plt.xlim(4,16)
        #plt.ylim(0,7)
        plt.title(title)
        #plt.plot(S,blsCallOptionPrice,label="Black Scholes European")
        ax.plot(S, calloptionPrice_EFM, label="EFM")
        ax.plot(S, calloptionPrice_IFM, label="IFM")
        ax.plot(S, calloptionPrice_CNFM,"k--", label="CNFM")
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc="center left", bbox_to_anchor=[1, 0.5], ncol=1, shadow=True, title="Legend", fancybox=True)
        #plt.legend(loc="upper left", bbox_to_anchor=[0, 1], ncol=2, shadow=True, title="Legend", fancybox=True)
        plt.savefig(filename2 + ".png")
        plt.clf()



    print("Q2.a American Put Option Price for S= %f, dX = %f"%(S[indS0],dS[i]))
    #print("     Price = %f [BS Method]" % (blsPutOptionPrice[indS0]))
    print("     Price = %f [Explicit Method]"%(putoptionPrice_EFM[indS0]))
    print("     Price = %f [Implicit Method]"%(putoptionPrice_IFM[indS0]))
    print("     Price = %f [Crank-Nicolson Method]"%(putoptionPrice_CNFM[indS0]))

    print("Q2.a American Call Option Price for S= %f, dX = %f" % (S[indS0],dS[i]))
    #print("     Price = %f [BS Method]" % (blsCallOptionPrice[indS0]))
    print("     Price = %f [Explicit Method]" % (calloptionPrice_EFM[indS0]))
    print("     Price = %f [Implicit Method]" % (calloptionPrice_IFM[indS0]))
    print("     Price = %f [Crank-Nicolson Method]" % (calloptionPrice_CNFM[indS0]))
    del fdm_class

