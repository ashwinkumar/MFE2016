import myHalton as hl
import numpy as np
import math
import statistics as st
import timeit
import BlackScholesOptionPrice as bp
#### Question 6

CONST_MIN = 0.00001

def read(type):
    while True:
        if type =="S0":
            x = float(input("Enter S0: "))
            if x <=0 :
                raise Exception('S0 should be positive number')
            return x
        elif type == "K":
            x = float(input("Enter K: "))
            if x <= 0:
                raise Exception('K should be positive number')
            return x
        elif type == "T":
            x = float(input("Enter T: "))
            if x <= 0:
                raise Exception('T should be positive number')
            return x
        elif type == "r":
            x = float(input("Enter r: "))
            return x
        elif type == "sigma":
            x = float(input("Enter sigma: "))
            if x <= 0:
                raise Exception('sigma should be positive number')
            return x
        elif type == "N":
            x = int(input("Enter N of points: "))
            if x <= 0:
                raise Exception('N should be positive integer')
            return x
        elif type == "b1":
            x = int(input("Enter b1: "))
            if x <= 0:
                raise Exception('b1 should be positive integer')
            return x
        elif type == "b2":
            x = int(input("Enter b2: "))
            if x <= 0:
                raise Exception('b2 should be positive integer')
            return x

def generateNormalbyBoxMuller(halton1, halton2,n):
    num = 0
    while num <n:
        u1 = max(CONST_MIN, halton1.generateNextRdmNumber())
        u2 = max(CONST_MIN, halton2.generateNextRdmNumber())
        try:
            n1 = np.sqrt(-2 * np.log(u1)) * np.cos(2. * np.pi * u2)
            n2 = np.sqrt(-2 * np.log(u1)) * np.sin(2. * np.pi * u2)
        except ZeroDivisionError:
            print(u1, u2)

        yield n1
        num += 1
        if num <n:
            yield n2
        num +=1

def payOff(S0,K,r,sigma,T,z):
    tmp = S0*math.exp((r - sigma*sigma/2)*T)
    sqrt = math.sqrt(T)
    done = object()
    n1 = next(z, done)

    while (n1 is not done):
        rdm = math.exp(sigma*n1*sqrt)
        rdm2 = math.exp(sigma*(-n1)*sqrt)
        yield max(0, tmp*rdm - K)
        yield max(0,tmp*rdm2 - K)
        n1 = next(z, done)


def MC_EuropeanCall(S0,K,r,sigma,T,z):
    c_t = list(payOff(S0,K,r,sigma,T,z))
    disc = math.exp(-r*T)
    return disc*st.mean(c_t), disc*disc*st.variance(c_t)/len(c_t)



S0 = read("S0")
K = read("K")
T = read("T")
r = read("r")
sigma = read("sigma")
N = read("N")
b1 = read("b1")
b2 = read("b2")

start_time1 = timeit.default_timer()
'''
S0 = 100
K = 105
T =2
r = 0.05
sigma = 0.3
N = 100000
b1 = 2
b2 = 7
'''
halton1 = hl.Halton(b1)
halton2 = hl.Halton(b2)
z = generateNormalbyBoxMuller(halton1,halton2,N)
call, var_call = MC_EuropeanCall(S0,K,r,sigma,T,z)
elapsed1 = timeit.default_timer() - start_time1

bs = bp.EuropeanCallPrice(S0,K,r,sigma,T)
print('Q6. Price of European Call Option [MC using Halton seq] = %f with variance = %f for %d simulations *****[%f sec]' %(call,var_call,N, elapsed1))
print('    Using BS Model %f', bs)
