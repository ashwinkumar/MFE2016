import math
import numpy as np

'''
                S
            /      \
          u*S      d*S
        /    \    /    \
     u^2*S   u*d*S    d^2*S
'''

class Node:
    def __init__(self, s0):
        self.val =s0
        self.no_ex = -1

def getBinomialParameters(r,sigma,dt, method=1):
    if method == 1:
        return getBinomialParameters1(r,sigma,dt)
    elif method == 2:
        return getBinomialParameters2(r,sigma,dt)
    elif method == 3:
        return getBinomialParameters3(r,sigma,dt)
    elif method == 4:
        return getBinomialParameters4(r, sigma, dt)
    else:
        raise Exception('Unknown method. Please give one of the following options for method [1,2,3,4]')

def getBinomialParameters1(r, sigma, dt):
    g = (np.exp(-r * dt) + np.exp((r + sigma * sigma) * dt)) / 2
    l = np.sqrt(g * g - 1)
    u = g + l
    d = g - l
    p = (np.exp(r * dt) - d) / (u - d)
    return u, d, p

def getBinomialParameters2(r, sigma, dt):
    tmp = np.sqrt(np.exp(sigma * sigma * dt) - 1)
    disc = np.exp(r * dt)
    u = disc * (1 + tmp)
    d = disc * (1 - tmp)
    p = 0.5 * np.ones(len(u))
    return u, d, p

def getBinomialParameters3(r, sigma, dt):
    drift = np.exp((r - sigma * sigma / 2) * dt)
    rd = np.exp(sigma * np.sqrt(dt))
    u = drift * rd
    d = drift / rd
    p = 0.5 * np.ones(len(u))
    return u, d, p

def getBinomialParameters4(r, sigma, dt):
    p = 1 / 2 + np.sqrt(dt) * (r - sigma * sigma / 2) / (2 * sigma)
    rd = np.exp(sigma * np.sqrt(dt))
    u = rd
    d = 1 / rd
    return u, d, p


class BinomialModel:
    def __init__(self, s0,r,t): # u, d, r,t,n):
        self.s0= s0
        self.r = r
        self.t = t
        self.resetParams()
        '''
        self.n = n
        self.verify(u,d)
        '''

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
        params = [self.__u, self.__d,self.__p, self.__n, self.__dt, self.__disc]
        return all(v is not None for v in params)

    def verify(self,u,d,p):

        if u < 0 or d < 0:
            raise AttributeError('u & d cannot be negative')
        if (d >= self.__disc) or (u <= self.__disc):
            raise AttributeError('Combination of u,r,d cannot be used.')
        if p < 0 or p > 1:
            raise AttributeError('p is not between 0 and 1')
        self.__u = u
        self.__d = d
        self.__p = p
        return True

    def preCalculate(self):
        self.__upow = np.ones(self.__n+1, dtype=np.float)
        self.__dpow = np.ones(self.__n+1, dtype=np.float)

        for i in range(1,self.__n+1):
            self.__upow[i]= self.__upow[i-1]* self.__u
            self.__dpow[i] = self.__dpow[i-1]*self.__d


    def generateTree(self):
        self.preCalculate()
        for i in range(0,self.__n+1):
            for j in range(0,i+1):
                yield self.__s0*self.__upow[i-j]*self.__dpow[j]


    def resetParams(self):
        self.__u = None
        self.__d = None
        self.__p = None
        self.__n = None
        self.__dt = None
        self.__disc = None

    '''
    def getSt(self,t):
        num = 0
        tmp = self.__S0*math.pow(self.__u,t)
        while num<=t:
            yield tmp
            tmp = (tmp/self.__u)*self.__d
            num=num+1
    '''

    def getOptionPrice(self, payOff_Function,X, type="E"):
        if not self.isParamsSet():
            raise Exception('Unset parameters for Binomial Model')
        st = np.asarray(list(self.generateTree()))
        option = np.zeros(len(st))
        intrinsic = np.zeros(len(st))
        num = 0
        s1 = (int)((self.__n)*(self.__n+1)/2)
        e1 = (int)(s1+ self.__n)
        option[s1:(e1+1)]=payOff_Function(st[s1:(e1+1)],X)
        intrinsic = payOff_Function(st,X)
        while( num < self.__n):
            i = (self.__n) - (num+1)
            s2 = (int)(i* (i+1)/2)
            e2 = (int)(s2+ i)
            k=0
            for j in range(s2,e2+1):
                option[j] = (option[s1 + k] * self.__p + option[s1 + k + 1] * (1 - self.__p)) / self.__disc
                if(type == "A"):
                    option[j] = max(intrinsic[j], option[j])
                k= k+1
            s1 = s2
            num = num+1

        return option[0]


    def getDelta(self,payOff_Function,X):
        if not self.isParamsSet():
            raise Exception('Unset parameters for Binomial Model')
        st = np.asarray(list(self.generateTree()))
        option = np.zeros(len(st))
        delta = np.zeros(len(st) - self.__n-1)
        num = 0
        s1 = (int)((self.__n) * (self.__n + 1) / 2)
        e1 = (int)(s1 + self.__n)
        option[s1:(e1+1)] = payOff_Function(st[s1:(e1+1)], X)

        ## calculate u-d
        delta_ud = self.__s0*(self.__u - self.__d)
        while (num < self.__n):
            i = (self.__n) - (num + 1)
            s2 = (int)(i * (i + 1) / 2)
            e2 = (int)(s2 + i)
            k = 0
            for j in range(s2, e2 + 1):
                option[j] = (option[s1 + k] * self.__p + option[s1 + k + 1] * (1 - self.__p)) / self.__disc
                # delta = (Cu - Cd)/(St*(u-d))
                delta[j] = (option[s1 + k] - option[s1 + k + 1]) /(st[j]*delta_ud)
                k = k + 1
            s1 = s2
            num = num + 1

        return delta[0]


    '''
            theta = (C(+-) - C )/2∆t
    '''
    def getTheta(self, payOff_Function,X):
        if not self.isParamsSet():
            raise Exception('Unset parameters for Binomial Model')
        st = np.asarray(list(self.generateTree()))
        option = np.zeros(len(st))
        theta = np.zeros(len(st)) #  We will need less indices (2 steps)
        num = 0
        s1 = (int)((self.__n) * (self.__n + 1) / 2)
        e1 = (int)(s1 + self.__n)
        '''
        s2 = (int)( (self.__n)*(self.__n -1)/2)
        e2 = (int)(s2 + self.__n-1)
        '''
        option[s1:(e1+1)] = payOff_Function(st[s1:e1], X)

        while (num < self.__n):
            i = (self.__n) - (num + 1)
            i2 = (self.__n) - (num - 1)

            s2 = (int)(i * (i + 1) / 2)
            e2 = (int)(s2 + i)

            s3 = (int)(i2*(i2+1)/2)
            e3 = (int)(s3 + i2)

            k = 0
            for j in range(s2, e2 + 1):
                option[j] = (option[s1 + k] * self.__p + option[s1 + k + 1] * (1 - self.__p)) / self.__disc
                # theta = (Cud - C)/2dt
                if s3 <= len(st) and e3 <= len(st):
                    theta[j] = (option[s3 + k+1]- option[j])/(2*self.__dt)
                k = k + 1
            s1 = s2
            num = num + 1

        return theta[0]


