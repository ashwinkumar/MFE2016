import numpy as np
#import mynormal as nl
import random as rndom
import math



def Monomials(k,x):
    num = 0
    l = len(x)
    p = np.ones(l)
    while num <k:
        yield p
        p = p*x
        num = num+1

def Hermite(k,x):
    num = 0
    l = len(x)
    f = np.ones(l)
    s = 0*np.ones(l)
    p = np.ones(l)
    yield p

    while num <(k-1):
        if num ==1:
            p = 2*f*x - 2
        elif num ==2:
            p = 2*f*x - 2*(num)*s
        else:
            p = 2 * f * x - 2 * (num+1) * s
        s = f
        f = p
        yield p
        num = num+1

def Laguerre(k,x):
    num = 0
    if k >4:
        raise Exception('k should be less than equal 4')
    ex_p= np.exp(-x/2)
    while num < k:
        if num == 0:
            yield ex_p
        elif num== 1:
            yield np.multiply(ex_p,(1-x))
        elif num == 2:
            yield np.multiply(ex_p, (1 - 2*x + x*x/2))
        elif num == 3:
            x_sq = x*x
            yield np.multiply(ex_p,(1- 3*x + 3*x_sq/2 - x_sq*x/6))
        num = num + 1


class LeastSquareMC:
    list_methods = ["Hermite","Laguerre", "Monomials"]

    def __init__(self, S0,r, sigma, T, dT, nsims, k, method = "Monomials"):
        self.S0 = S0
        self.__r = r
        self.__sigma = sigma
        self.T = T
        #self.t = t
        self.dT = dT
        self.nsims = nsims
        self.k = k
        self.method = method
        ndiv = (int)(T/dT)
        self.ndiv = ndiv
        self.__tree = np.zeros([nsims,ndiv+1])
        self.__exercise = -1*np.ones(nsims, dtype= np.int)  # The time of exercise is stored. Intiated -1
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
    def T(self):
        return self.__T

    @T.setter
    def T(self, T):
        if T < 0:
            raise AttributeError('T cannot be negative')
        else:
            self.__T = T
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
        if nsims < 0 or not isinstance(nsims,int):
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

    @property
    def k(self):
        return self.__k

    @k.setter
    def k(self, k):
        if k <= 0 or k> 4 :#or (not isinstance(k,np.int64)):
            raise AttributeError('k should be integer between [1,4]')
        else:
            self.__k = k

    @property
    def method(self):
        return self.__method

    @method.setter
    def method(self, method):
        if method in self.list_methods:
            self.__method = method
        else:
            raise AttributeError('unknown orthogonal function families')


    def setPayOff(self,OptionPayOff):
        self.__payoff = OptionPayOff

    def setInMoney(self, InTheMoney):
        self.__inthemoney = InTheMoney

    def isParamsSet(self):
        params = [self.__T, self.__dT, self.__r, self.__sigma]
        return all(v is not None for v in params)

    def preCal_MC(self):
        if not self.isParamsSet():
            raise Exception('Required params not set. Initialize correctly')
        self.__edrift = math.exp((self.__r - self.__sigma**2/2)* self.__dT)
        self.__wt = self.__sigma*math.sqrt(self.__dT)

        ert = math.exp(-self.__r*self.__dT)
        for i in range(self.__ndiv):
            self.__disc[i+1]= self.__disc[i]*ert

    def mcStock_t(self, t):
        #if len() != self.__nsims:
        #    raise Exception('Fatal error. S_t should be of length of nsims')
        #mynormal = nl.NormalDistribution("Polar-Marsaglia", rndom.randint(10, 10 * self.__ndiv))
        #mynormal.generateNormalDistribution(self.__nsims)
        #z = np.asarray(list(mynormal.getRdmNumbers()))
        z1 = np.random.normal(0,1,self.__nsims/2)
        z2 = -z1
        z= np.concatenate((z1,z2))
        #done = object()
        #z = next(ndist, done)
        #while (z is not done):

        yield  np.multiply(self.__tree[:,t],np.exp(self.__wt*z))*self.__edrift
        #    z = next(ndist, done)

    def mcStockTree(self):
        self.preCal_MC()
        self.__tree[:,0] = self.__S0*np.ones(self.__nsims)
        #self.__treemaximum[:,0] = self.__S0*np.ones(self.__nsims)
        for i in range(1,(self.__ndiv+1)):
            self.__tree[:,i] = np.asarray(list(self.mcStock_t(i-1)))
            #self.__treemaximum[:,i] = np.maximum(self.__treemaximum[:,(i-1)],self.__tree[:,i])

    def getStrikePrice(self,t):
        if t > self.__T:
            raise Exception('t should be less than the Time to expiry T')
        ind = (int)(t/self.__T* self.__ndiv)
        return self.__tree[:,ind]

    def getorthogalFamily(self,x):
        if self.__method == "Monomials":
            return Monomials(self.__k,x)
        elif self.__method == "Hermite":
            return Hermite(self.__k,x)
        elif self.__method == "Laguerre":
            return Laguerre(self.__k, x)
        else:
            raise Exception('Unknown method')

    def getECV(self,x,Y):
        l = list(self.getorthogalFamily(x))
        sizeofX = len(x)
        A = np.zeros((self.__k,self.__k))
        b = np.zeros((self.__k,1))
        ECV = np.zeros(sizeofX)
        for i in range(self.__k):
            b[i] = np.dot(l[i],Y)
            for j in range(self.__k):
                A[i][j]= np.dot(l[i],l[j])
        if (np.linalg.det(A) ==0):
            return ECV
        a = np.dot(np.linalg.inv(A), b)

        for i in range(self.__k):
            ECV = ECV+ a[i]*l[i]
        return ECV



    def getStockInMoney(self,t):
        if not isinstance(t, int) or t < 0 or t > self.__nsims:
            raise Exception('Index not match')
        indices = self.__inthemoney(self.__tree[:,(t+1)])
        np.take(self.__tree[:,t],indices)

    '''
    def getInd(self,S,X,ECV):
        temp = self.__payoff(self.__tree[:,t],X)
        ind = np.where(self.__payoff(self.__tree[:,t],X)>ECV)
        return np.asarray(ind).flatten()
    '''

    def calcualteOptionPrice(self,X,t):

        self.__exercise[self.__inthemoney(self.__tree[:,self.__ndiv],X)] = self.__ndiv
        ECV = np.zeros(self.__nsims)
        curr_time_interval = (int)(t/self.__T*self.__ndiv)
        for i in range(self.__ndiv-1, curr_time_interval,-1):
            itm_ind = np.asarray(self.__inthemoney(self.__tree[:,i],X)).flatten()
            x = np.take(self.__tree[:, i],itm_ind)
            exercise_t = np.take( self.__exercise,itm_ind)
            y = np.zeros(len(itm_ind))
            for j in range(len(itm_ind)):
                if exercise_t[j] == -1:
                    y[j]=0
                else:
                    y[j] = self.__payoff(self.__tree[itm_ind[j],exercise_t[j]],X[itm_ind[j]])*self.__disc[exercise_t[j]-i]
            if(len(x)==0):
                continue

            ECV = self.getECV(x,y)
            for j in range(len(itm_ind)):
                if self.__payoff(x[j],X[itm_ind[j]]) > ECV[j]:
                    self.__exercise[itm_ind[j]]=i

        sum = 0
        for i in range(self.__nsims):
            ind= self.__exercise[i]
            if(ind != -1):
                sum = sum+ self.__payoff(self.__tree[i,ind],X[i])* self.__disc[ind]

        return sum/self.__nsims
    '''
    def calcualteOptionPrice(self):

        self.__exercise.fill(self.__ndiv)
        ECV = np.zeros(self.__nsims)
        #curr_time_interval = (int)(self.__t / self.__T * self.__ndiv)
        for i in range(self.__ndiv - 1, 0, -1):
            ## All the option prices are in the money
            #itm_ind = np.asarray(self.__inthemoney(self.__tree[:, i], X)).flatten()
            x = self.__tree[:, i]
            y = np.zeros(self.__nsims)
            for j in range(self.__nsims):
                if self.__exercise[j] == -1:
                    y[j] = 0
                else:
                    y[j] = self.__payoff(self.__treemaximum[j, self.__exercise[j]] , self.__tree[j,self.__exercise[j]])* self.__disc[self.__exercise[j] - i]
            #y = self.__payoff(self.__treemaximum[:,self.__exercise], self.__tree[:,self.__exercise])*self.__disc[self.__exercise - i]

            ECV = self.getECV(x, y)
            for j in range(self.__nsims):
                if (self.__treemaximum[j,i]-self.__tree[j,i]) > ECV[j]:
                    self.__exercise[j] = i
        sum = 0
        for i in range(self.__nsims):
            ind = self.__exercise[i]
            if (ind != -1):
                sum = sum + (self.__treemaximum[i,ind]- self.__tree[i,ind]) * self.__disc[ind]

        return sum / self.__nsims
    '''