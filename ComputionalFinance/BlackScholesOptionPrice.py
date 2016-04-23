import  math
from scipy.stats import norm



#(CONST_D1, CONST_D2, CONST_D3, CONST_D4, CONST_D5, CONST_D6) = (0.0498673470, 0.0211410061, 0.0032776263, 0.0000380036, 0.0000488906, 0.0000053830)

CONST_D1 = 0.0498673470
CONST_D2 = 0.0211410061
CONST_D3 = 0.0032776263
CONST_D4 = 0.0000380036
CONST_D5 = 0.0000488906
CONST_D6 = 0.0000053830

def N(d):
    d_2 = d*d
    d_3 = d_2*d
    d_4 = d_3*d
    d_5 = d_4*d
    d_6 = d_5*d
    if d >= 0 :
        return (1- 0.5*math.pow((1+ CONST_D1*d + CONST_D2*d_2 + CONST_D3*d_3 + CONST_D4*d_4 + CONST_D5*d_5 + CONST_D6*d_6),-16))
    if d <0 :
        return(1- N(-d))
def n(d):
    return 1/math.sqrt(2*math.pi)*math.exp(-d*d/2)

def getd1(S0, X, r, sigma, t):
    return (math.log(S0/X) + (r + sigma*sigma/2) * t )/ (sigma * math.sqrt(t))

def getd2(S0, X, r, sigma, t):
    return getd1(S0, X, r , sigma, t) - sigma* math.sqrt(t)

def EuropeanCallPrice(S0, X, r , sigma, t):
    d1 = getd1(S0, X, r, sigma, t)
    d2 = getd2(S0, X, r, sigma, t)
    #return S0* norm.cdf(d1) - X* math.exp(-r*t)* norm.cdf(d2)
    return S0* N(d1) - X* math.exp(-r*t)* N(d2)

def getDelta(S0,X,r,sigma,t):
    return N(getd1(S0,X,r,sigma,t))

def getGamma(S0,X,r,sigma,t):
    d1 = getd1(S0,X,r,sigma,t)
    return n(d1)/(S0*sigma*math.sqrt(t))

def getTheta(S0,X,r,sigma,t):
    d1 = getd1(S0,X,r,sigma,t)
    d2 = getd2(S0,X,r,sigma,t)
    return -S0*sigma*n(d1)/(2*math.sqrt(t))  - r*X* math.exp(-r*t)*N(d2)

def getVega(S0,X,r,sigma,t):
    d1 = getd1(S0, X, r, sigma, t)
    return S0*math.sqrt(t)*n(d1)

def getRho(S0,X,r,sigma,t):
    d2 = getd2(S0,X,r,sigma,t)
    return X*t*math.exp(-r*t)*N(d2)




