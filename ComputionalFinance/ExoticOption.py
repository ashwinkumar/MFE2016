import numpy as np
import math
import PoissonJump as jump

import matplotlib.pyplot as plt
class LookBackOption:


    def __init__(self, S0, K, r, sigma, T, dT, nsims):
        self.S0 = S0
        self.K = K
        self.__r = r
        self.__sigma = sigma
        self.T = T
        self.dT = dT
        self.nsims = nsims
        ndiv = (int)(T/dT)
        self.ndiv = ndiv
        self.__tree = np.zeros([nsims,ndiv+1])
        self.__maximum = np.zeros(nsims)
        self.__minimum = np.zeros(nsims)

        self.__disc = np.ones(ndiv+1)
        #self.__treemaximum = np.zeros([nsims,ndiv+1])


    @property
    def S0(self):
        return self.__S0

    @S0.setter
    def S0(self, S0):
        if S0 < 0:
            raise AttributeError('S0 cannot be negative')
        else:
            self.__S0 = S0

    @property
    def K(self):
        return self.__K

    @K.setter
    def K(self, K):
        if K < 0:
            raise AttributeError('K cannot be negative')
        else:
            self.__K = K


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
        z1 = np.random.normal(0, 1, self.__nsims / 2)
        z2 = -z1
        z = np.concatenate((z1, z2))
        yield np.multiply(self.__tree[:, t], np.exp(self.__wt * z)) * self.__edrift

    def mcStockTree(self):
        self.preCal_MC()
        self.__tree[:, 0] = self.__S0 * np.ones(self.__nsims)
        for i in range(1, (self.__ndiv + 1)):
            self.__tree[:, i] = np.asarray(list(self.mcStock_t(i - 1)))

    def getFixedStrikeCall(self):
        self.__maximum= np.max(self.__tree,1)
        payOff = (self.__maximum - self.K )* self.__disc[self.__ndiv+1]
        return  np.mean(payOff), np.var(payOff,ddof=1)/(len(payOff))

    def getFixedStrikeCall(self):
        self.__maximum = np.max(self.__tree, 1)
        payOff = (self.__maximum - self.K) * self.__disc[self.__ndiv]
        return np.mean(payOff), np.var(payOff, ddof=1) / (len(payOff))


    def getFixedStrikePut(self):
        self.__minimum = np.min(self.__tree, 1)
        payOff = (self.K- self.__minimum) * self.__disc[self.__ndiv]
        return np.mean(payOff), np.var(payOff, ddof=1) / (len(payOff))


class DefaultOption:
    def __init__(self, T,ndiv, nsims):
        self.T = T
        self.ndiv = ndiv
        self.dT = T / ndiv
        self.__V = np.zeros((nsims,ndiv+1))
        self.__L = np.zeros(ndiv+1)
        self.__q = np.zeros(ndiv+1)
        self.nsims = nsims

    def initCollateral(self,V0,mu,sigma,gamma,lambda1):
        self.V0 = V0
        self.__mu = mu
        self.__sigma = sigma
        self.__gamma = gamma
        self.lambda1 = lambda1

    def initLoan(self,L0,r0,delta,lambda2):
        self.L0 = L0
        if(r0 < 0):
            raise Exception(" R0 cannot be negative")

        R = r0 + delta*lambda2
        if( self.__ndiv is None or self.__T is None):
            raise Exception("Contract rate period is not defined")

        r = R/self.__ndiv
        n = self.__T*self.__ndiv
        compd = (1+r)**n
        PMT = L0*r/(1 - 1/compd)

        self.a = PMT/r
        self.b = PMT/(r*compd)
        self.c = 1+r
        self.lambda2 = lambda2



    @property
    def L0(self):
        return self.__L0

    @L0.setter
    def L0(self, L0):
        if L0 < 0:
            raise AttributeError('L0 cannot be negative')
        else:
            self.__L0 = L0

    @property
    def V0(self):
        return self.__V0

    @V0.setter
    def V0(self, V0):
        if V0 < 0:
            raise AttributeError('V0 cannot be negative')
        else:
            self.__V0 = V0

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
        if ndiv <= 0 or not isinstance(ndiv, np.int64):
            raise AttributeError('ndiv is not positive integer')
        else:
            self.__ndiv = ndiv


    @property
    def lambda1(self):
        return self.__lambda1

    @lambda1.setter
    def lambda1(self, lambda1):
        if lambda1< 0:
            raise AttributeError('lambda1 cannot be negative')
        else:
            self.__lambda1 = lambda1

    @property
    def lambda2(self):
        return self.__lambda2

    @lambda2.setter
    def lambda2(self, lambda2):
        if lambda2 < 0:
            raise AttributeError('lambda2 cannot be negative')
        else:
            self.__lambda2 = lambda2

    @property
    def gamma(self):
        return self.__gamma

    @gamma.setter
    def gamma(self, gamma):
        if gamma > 0:
            raise AttributeError('gamma cannot be positive')
        else:
            self.__gamma = gamma

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, a):
        if a < 0:
            raise AttributeError('a cannot be negative')
        else:
            self.__a = a

    @property
    def b(self):
        return self.__b

    @b.setter
    def b(self, b):
        if b < 0:
            raise AttributeError('b cannot be negative')
        else:
            self.__b = b

    @property
    def c(self):
        return self.__c

    @c.setter
    def c(self, c):
        if c < 1:
            raise AttributeError('c cannot be less than 1')
        else:
            self.__c = c

    @property
    def alpha(self):
        return self.__alpha

    @alpha.setter
    def alpha(self, alpha):
        if self.__V0 is None or self.__L0 is None:
            raise Exception('V0 and L0 is not set. Initialize error')
        if alpha > (self.__V0/self.__L0):
            raise AttributeError('alpha cannot be greater than V0/L0')
        else:
            self.__alpha =alpha

    @property
    def beta(self):
        return self.__beta

    @beta.setter
    def beta(self, beta):
        if beta< 0:
            raise AttributeError('beta cannot be less than 0')
        else:
            self.__beta = beta

    def isCollateralParamsSet(self):
        params = [self.__T, self.__dT, self.__mu, self.__sigma, self.__gamma]
        return all(v is not None for v in params)

    def isLoanParamsSet(self):
        params = [self.__L0,self.__a, self.__b, self.__c]
        return all(v is not None for v in params)


    def simulateCollateral(self):
        if not self.isCollateralParamsSet() :
            raise Exception('Required params not set. Initialize correctly')
        self.__V[:,0] = self.__V0
        jdt = np.zeros((self.__nsims,self.__ndiv+1))
        poissonjump = jump.PoissonJump(self.__lambda1,self.__T,self.__ndiv)
        for i in range(self.__nsims):
            jdt[i,:] =np.asarray(list(poissonjump.getJt()))

        for i in range(1,(self.__ndiv+1)):
            wdt = self.__dT*np.random.normal(0, 1, self.__nsims)
            self.__V[:,i] = np.multiply(self.__V[:,i-1],(1+self.__mu*self.__dT + self.__sigma*wdt + self.__gamma*jdt[:,i]))

        '''
        t = np.arange(0, self.__T + self.__dT / 2, self.__dT)
        for i in range(self.__nsims):
            plt.plot(t,self.__V[i,:])
        plt.show()
        '''
    def simulateLoanValue(self):
        if not self.isLoanParamsSet():
            raise Exception('Required params not set. Initialize correctly')
        self.__L[0] = self.__L0
        for i in range(1,self.__ndiv+1):
            n = self.__T*i
            self.__L[i] = self.a - self.b*(self.c**n)

    def setThreshold(self,alpha,ephsilon):
        self.alpha = alpha
        self.__ephsilon = ephsilon
        self.__beta = (ephsilon - alpha) / self.__T

        t = np.arange(0,self.__T+self.__dT/2, self.__dT)
        self.__q = self.__alpha + self.__beta*t
    '''
    Defualt time Q is the first time at which Vt <= qt*Lt
    Default time S is the first time at which some adverse event happens (losing job etc)
    Default tim E is the minimum of Q,S
    Default_type 0 - No default , 1 - Default when Q , 2- Default when S
    '''
    def getDefaultTime(self):
        Q_t = -1*np.ones(self.__nsims)
        S_t = -1 * np.ones(self.__nsims)
        E_t = -1 * np.ones(self.__nsims, dtype=np.int)
        type = np.zeros(self.__nsims, dtype=np.int)
        poissonjump = jump.PoissonJump(self.__lambda2, self.__T, self.__ndiv)


        for i in range(self.__nsims):
            val = self.__V[i,:] - self.__L* self.__q
            S_t[i] = poissonjump.firstJumpTime()
            ind = np.argmax(val <0)
            if ind ==0:
                Q_t[i] = self.__ndiv+1 # No exercise
            else:
                Q_t[i] = ind
            if(S_t[i] < Q_t[i]):
                E_t[i] = S_t[i]
                type[i] =2
            elif (Q_t[i] < S_t[i]) :
                E_t[i] = Q_t[i]
                type[i] = 1
            else:
                E_t[i] = -1
                type[i] =0

        return E_t,type

    def getDefaultOptionPrice(self,r0):
        exercise_time,exercise_type = self.getDefaultTime()
        payOff = np.zeros(self.__nsims)
        for i in range(self.__nsims):
            t = exercise_time[i]
            disc = math.exp(-r0*t)
            type = exercise_type[i]
            if(t <0 ):
                payOff[i] =0
            if(type == 1):
                payOff[i]= max(0,self.__L[t]- self.__ephsilon *self.__V[i,t])*disc
            else:
                payOff[i]= math.fabs(self.__L[t]- self.__ephsilon *self.__V[i,t])*disc
        return np.mean(payOff), np.var(payOff, ddof=1)/len(payOff)


    def getDefaultProb(self):
        exercise_time,_ = self.getDefaultTime()
        return ((exercise_time > 0).sum()) / self.__nsims


    def getExpectedExerciseTime(self):
        exercise_time,_ = self.getDefaultTime()
        ind  = np.asarray(np.where(exercise_time >0 )).flatten()
        sum = 0
        for i in range(len(ind)):
            sum += exercise_time[ind[i]]
        return (self.__dT*sum) / len(ind)
