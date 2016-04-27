import math
import numpy as np

'''
                               S
            /                  |                    \
         u*S                   S                    d*S
    /     |    \          /    |    \         /      |       \
u^2*S    u*S    u*d*S   u*S    S    d*S    u*d*s    d*S    d^2*S

'''
def getBinomialParameters(r,sigma,dt, method=1):
    if method == 1:
        return getBinomialParameters1(r,sigma,dt)
    elif method == 2:
        return getBinomialParameters2(r,sigma,dt)
    else:
        raise Exception('Unknown method. Please give one of the following options for method [1,2]')

def getBinomialParameters1(r, sigma, dt):

    u = np.exp(sigma*np.sqrt(3*dt))
    d = 1/u
    pd = (r*dt*(1-u) + r*dt*r*dt + sigma*sigma*dt)/((u-d)*(1-d))

    pu = (r*dt*(1-d) + r*dt*r*dt + sigma*sigma*dt)/((u-d)*(u-1))
    pm = 1- pu - pd
    return u, d, pu,pm,pd


def getBinomialParameters2(r, sigma, dt):

    u = sigma*np.sqrt(3*dt)
    d = -u

    drift = (r- sigma*sigma/2)*dt
    tmp1= (sigma*sigma*dt + drift*drift)/(u*u)
    tmp2 = drift/u
    pd = (tmp1- tmp2)/2
    pu = (tmp1+tmp2)/2
    pm = 1- pu - pd
    return u, d, pu,pm,pd

class TrinomialModel:
    def __init__(self, s0,r,t):
        self.s0= s0
        self.r = r
        self.t = t
        self.resetParams()

    @property
    def s0(self):
        return self.__s0

    @s0.setter
    def s0(self, s0):
        if s0 < 0 :
            raise AttributeError('S0 cannot be negative')
        self.__s0 = s0

    @property
    def t(self):
        return self.__t

    @t.setter
    def t(self, t):
        if t <= 0:
            raise AttributeError('Time period cannot be negative')
        self.__t = t

    @property
    def n(self):
        return self.__n

    @n.setter
    def n(self, n):
        if n < 0:
            raise AttributeError('N of intervals cannot be negative')
        '''
        if isinstance(n, int):
            raise AttributeError('N of intervals has to be integer')
        '''
        if self.__t < 0:
            raise AttributeError('Time period cannot be negative')
        self.__n = n
        self.__dt = self.__t/self.__n
        if self.__r >0:
            self.__disc = math.exp(self.__r* self.__dt)

    def disc(self):
        return self.__disc

    @property
    def r(self):
        return self.__r

    @r.setter
    def r(self, r):
        if r < 0:
            raise AttributeError('r cannot be negative')
        self.__r = r

    def isParamsSet(self):
        params = [self.__u, self.__d,self.__pu,self.__pm, self.__pd, self.__n, self.__dt, self.__disc]
        return all(v is not None for v in params)

    def verify(self,u,d,pu,pm,pd):

        if u < 0 or d < 0:
            raise AttributeError('u & d cannot be negative')
        if u*d != 1:
            raise AttributeError('u*d should be 1')
        if d >u:
            raise AttributeError('d should be greater than u')
        if (d >= self.__disc) or (u <= self.__disc):
            raise AttributeError('Combination of u,r,d cannot be used.')
        if pu < 0 or pu > 1 or pm <0 or pm >1 or pd <0 or pd >1:
            raise AttributeError('p s not between 0 and 1')
        self.__u = u
        self.__d = d
        self.__pu = pu
        self.__pm = pm
        self.__pd = pd
        return True

    def verifylog(self, u, d, pu, pm, pd):

        if u <0 or d>0 or u<d:
            raise AttributeError('Incorrect u & d')
        if u != -d:
            raise AttributeError('u+d should be 0')
        if (1+d >= self.__disc) or (1+ u <= self.__disc):
            raise AttributeError('Combination of u,r,d cannot be used.')
        if pu < 0 or pu > 1 or pm < 0 or pm > 1 or pd < 0 or pd > 1:
            raise AttributeError('p s not between 0 and 1')
        self.__u = u
        self.__d = d
        self.__pu = pu
        self.__pm = pm
        self.__pd = pd
        return True


    def preCalculate(self):
        self.__upow = np.ones(self.__n+1, dtype=np.float)
        self.__dpow = np.ones(self.__n+1, dtype=np.float)

        for i in range(1,self.__n+1):
            self.__upow[i]= self.__upow[i-1]* self.__u
            self.__dpow[i] = self.__dpow[i-1]*self.__d

    def generateLogTree(self):
        for i in range(0, self.__n + 1):
            for j in range(0, 2 * i + 1):
                if i > j:
                    yield self.__s0 +(i-j)*self.__u
                elif j > i:
                    yield self.__s0 + (j-i)*self.__d
                else:
                    yield self.__s0

    def generateTree(self):
        self.preCalculate()
        for i in range(0,self.__n+1):
            for j in range(0,2*i+1):
                if i>j:
                    yield self.__s0*self.__upow[i-j]
                elif j>i:
                    yield self.__s0*self.__dpow[j-i]
                else:
                    yield self.__s0


    def resetParams(self):
        self.__u = None
        self.__d = None
        self.__pd = None
        self.__pu = None
        self.__pm = None
        self.__n = None
        self.__dt = None
        self.__disc = None


    def getOptionPrice(self, payOff_Function,X, type="E"):
        if not self.isParamsSet():
            raise Exception('Unset parameters for Binomial Model')
        st = np.asarray(list(self.generateTree()))
        option = np.zeros(len(st))
        intrinsic = np.zeros(len(st))
        num = 0
        s1 = (int)(self.__n*self.__n)
        e1 = (int)((self.__n+1)*(self.__n+1)-1)
        option[s1:(e1+1)]=payOff_Function(st[s1:(e1+1)],X)
        intrinsic = payOff_Function(st,X)
        while( num < self.__n):
            i = (self.__n) - (num+1)
            s2 = (int)(i*i)
            e2 = (int)((i+1)*(i+1)-1)
            k=0
            for j in range(s2,e2+1):
                option[j] = (option[s1 + k] * self.__pu + option[s1 + k + 1] *self.__pm  + (option[s1+k+2])*self.__pd) / self.__disc
                if(type == "A"):
                    option[j] = max(intrinsic[j], option[j])
                k= k+1
            s1 = s2
            num = num+1

        return option[0]

    def getOptionPricebyLog(self, payOff_Function, X):
        if not self.isParamsSet():
            raise Exception('Unset parameters for Binomial Model')
        xt = np.asarray(list(self.generateLogTree()))
        option = np.zeros(len(xt))
        num = 0
        s1 = (int)(self.__n * self.__n)
        e1 = (int)((self.__n + 1) * (self.__n + 1) - 1)
        st = np.exp(xt[s1:(e1+1)])
        option[s1:(e1 + 1)] = payOff_Function(st, X)
        while (num < self.__n):
            i = (self.__n) - (num + 1)
            s2 = (int)(i * i)
            e2 = (int)((i + 1) * (i + 1) - 1)
            k = 0
            for j in range(s2, e2 + 1):
                option[j] = (option[s1 + k] * self.__pu + option[s1 + k + 1] * self.__pm + (
                option[s1 + k + 2]) * self.__pd) / self.__disc
                k = k + 1
            s1 = s2
            num = num + 1
        return option[0]
