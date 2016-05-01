import math
import numpy as np
import mynormal as nl

def myWienerMC1(z,SQRT_T):
    done = object()
    n1 = next(z, done)
    while (n1 is not done) :
        wt = n1 * SQRT_T
        tmp = wt*wt + math.sin(wt)
        yield tmp
        n1 = next(z, done)



def myWienerMC1_vr(z, SQRT_T):
    done = object()
    n1 = next(z, done)
    while (n1 is not done) :
        wt = n1 * SQRT_T
        tmp1 = wt*wt + math.sin(wt)
        yield tmp1
        tmp2 = wt*wt + math.sin(-wt)
        yield tmp2
        n1 = next(z, done)



def myCos(z, sqrt_t):
    done = object()
    n1 = next(z, done)
    while (n1 is not done):
        yield np.cos(n1 * sqrt_t)
        n1 = next(z, done)




def myCos_vr(z, gamma, sqrt_t):
    done = object()
    n1 = next(z, done)
    while (n1 is not done):
        tmp1 = np.cos(n1 * sqrt_t)
        tmp2 = np.multiply(n1*n1,sqrt_t) - sqrt_t
        yield tmp1 -gamma*tmp2
        n1 = next(z, done)


def mc_norm_payoff(const_drift, const_sigt, const_x_s, z):

    done = object()
    n1 = next(z, done)
    num = 0
    while (n1 is not done):
        #tmp = math.exp(const_drift + const_sigt * n1 - const_x_s)
        yield max(0, math.exp(const_drift + const_sigt * n1) - const_x_s)
        n1 = next(z, done)
        num += 1

def mc_norm_payoff_variancered(const_drift, const_sigt, const_x_s, z):

    done = object()
    n1 = next(z, done)
    num = 0
    while (n1 is not done):
        yield max(0, math.exp(const_drift + const_sigt * n1) - const_x_s)
        yield max(0, math.exp(const_drift + const_sigt * -n1) - const_x_s)
        n1 = next(z, done)
        num += 1

def mc_callprice(S0, X, r, sigma, T, z, bool = False):
    const_x_s = X / S0
    disc = math.exp(-r*T)
    const_drift = (r - sigma * sigma / 2.0) * T
    const_sigt = sigma * math.sqrt(T)
    if bool == False:
        l_payoff = list(mc_norm_payoff(const_drift, const_sigt, const_x_s, z))
    else:
        l_payoff = list(mc_norm_payoff_variancered(const_drift, const_sigt, const_x_s, z))
    return disc*S0*(np.mean(l_payoff)) , (disc*S0*np.std(l_payoff, ddof=1)/ math.sqrt(len(l_payoff)))


def mc_stockprice(S0, r, sigma, T, nsims):
    normal = nl.NormalDistribution()
    for i in range(1,T+1):
        normal.generateNormalDistribution(nsims)
        z = normal.getRdmNumbers()
        const_drift = (r - sigma * sigma / 2.0) * i
        const_sigt = sigma * math.sqrt(i)
        l_payoff = list(mc_norm_payoff_variancered(const_drift, const_sigt, 0, z))
        yield S0 *(np.mean(l_payoff))

def mc_calc_exp(const_drift, const_sigt, z, S0):
    done = object()
    n1 = next(z, done)
    while (n1 is not done):
        S = S0*(math.exp(const_drift + const_sigt* n1))
        S0 = S
        yield S
        n1 = next(z, done)


def mc_stockpath(S0, r, sigma, T, ndiv, seed=14):
    ndiv = (int)(ndiv)
    normal = nl.NormalDistribution("Polar-Marsaglia",seed)
    normal.generateNormalDistribution(ndiv)

    dt = T /ndiv
    const_drift = (r - sigma * sigma / 2.0) * dt
    const_sigt = sigma * math.sqrt(dt)
    z = normal.getRdmNumbers()
    return mc_calc_exp(const_drift, const_sigt,z, S0)


class MCEuropeanOption:
    def __init__(self, S0, r, sigma,T, dT, nsims):
        self.S0 = S0
        self.__r = r
        self.__sigma = sigma
        self.T = T
        #self.t = t
        self.dT = dT
        self.nsims = nsims
        ndiv = (int)(T / dT)
        self.ndiv = ndiv
        self.__tree = np.zeros([nsims, ndiv + 1])
        #self.__treemaximum = np.zeros([nsims,ndiv+1])
        self.__disc = np.ones(ndiv+1)

    @property
    def S0(self):
        return self.__S0

    @S0.setter
    def S0(self, S0):
        if S0 < 0:
            raise AttributeError('S0 cannot be negative')
        else:
            self.__S0 = S0
    '''
    @property
    def t(self):
        return self.__t

    @t.setter
    def t(self, t):
        if t < 0 or t >self.__T:
            raise AttributeError('t cannot be negative or greater than expiry time T')
        else:
            self.__t = t
    '''
    @property
    def T(self):
        return self.__T

    @T.setter
    def T(self, T):
        if T < 0:
            raise AttributeError('T cannot be negative')
        else:
            self.__T = T

    @property
    def dT(self):
        return self.__dT

    @dT.setter
    def dT(self, dT):
        if dT < 0 or dT > self.__T:
            raise AttributeError('dT cannot be negative or grater than t')
        else:
            self.__dT = dT

    @property
    def nsims(self):
        return self.__nsims

    @nsims.setter
    def nsims(self, nsims):
        if nsims < 0 or not isinstance(nsims, int):
            raise AttributeError('n is not positive integer')
        else:
            self.__nsims = nsims

    @property
    def ndiv(self):
        return self.__ndiv

    @ndiv.setter
    def ndiv(self, ndiv):
        if ndiv <= 0 or not isinstance(ndiv, int):
            raise AttributeError('ndiv is not positive integer')
        else:
            self.__ndiv = ndiv

    def setPayOff(self, payoff):
        self.__payoff = payoff

    def isParamsSet(self):
        params = [self.__T, self.__dT, self.__r, self.__sigma]
        return all(v is not None for v in params)

    def preCal_MC(self):
        if not self.isParamsSet():
            raise Exception('Required params not set. Initialize correctly')
        self.__edrift = math.exp((self.__r - self.__sigma ** 2 / 2) * self.__dT)
        self.__wt = self.__sigma * math.sqrt(self.__dT)

        ert = math.exp(-self.__r * self.__dT)
        for i in range(self.__ndiv):
            self.__disc[i + 1] = self.__disc[i] * ert

    def mcStock_t(self, t):
        z = np.random.normal(0, 1, self.__nsims)
        yield np.multiply(self.__tree[:, t], np.exp(self.__wt * z)) * self.__edrift

    def mcStockTree(self):
        self.preCal_MC()
        self.__tree[:, 0] = self.__S0 * np.ones(self.__nsims)
        #self.__treemaximum[:,0] = self.__S0*np.ones(self.__nsims)
        for i in range(1, (self.__ndiv + 1)):
            self.__tree[:, i] = np.asarray(list(self.mcStock_t(i - 1)))
            #self.__treemaximum[:,i] = np.maximum(self.__treemaximum[:,(i-1)],self.__tree[:,i])

    def getStrikePrice(self, t):
        if t > self.__T:
            raise Exception('t should be less than the Time to expiry T')
        ind = (int)(t / self.__T * self.__ndiv)
        return self.__tree[:, ind]

    def getOptionPrice(self,X,t):
        payOff=self.__payoff(self.__tree[:,self.__ndiv],X)
        ind = (int)((1- t/self.__T)*self.__ndiv)
        disc= self.__disc[ind]
        return  disc*np.mean(payOff), disc*disc*np.var(payOff,ddof=1)/self.__nsims


