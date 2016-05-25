import numpy as np
import math
## log Normal prices
class LogFiniteDifferenceMethod:
    def __init__(self, X0, r, sigma, T,dT):
        self.X0 = X0
        self.__r = r
        self.__sigma = sigma
        self.T = T
        self.dT = dT

    def initLogStock(self, Xmin,Xmax,dX):
        self.Xmin = Xmin
        self.Xmax = Xmax
        self.dX = dX
        Smax = np.exp(Xmax)
        S0= np.exp(self.__X0)
        N1 = math.ceil((Xmax -  self.__X0)/dX)
        N2 = math.ceil((self.__X0- Xmin) / dX)
        self.__S = np.zeros(N1+N2+1)
        for i in range(N1+1):
           self.__S[N2+i] = self.__X0 + i*dX

        for i in range(1,N2+1):
            self.__S[N2 - i] = self.__X0 - i*dX
        self.__S = np.exp(self.__S)
        return self.__S, N2
        #self.__S = np.exp(np.arange(Xmin,Xmax, dX))


    def isStockTreeCreated(self):
        params = [self.__Xmin, self.__Xmax, self.__X0, self.__dX , self.__S]
        return all(v is not None for v in params)

    def isProbSet(self):
        params = [self.__pu, self.__pm, self.__pd]
        return (all(v is not None for v in params) & ((self.__pu + self.__pm + self.__pd) == 1))

    @property
    def X0(self):
        return self.__X0


    @X0.setter
    def X0(self, X0):
        if X0 < 0:
            raise AttributeError('X0 cannot be negative')
        else:
            self.__X0 = X0


    @property
    def Xmin(self):
        return self.__Xmin


    @Xmin.setter
    def Xmin(self, Xmin):
        if Xmin < 0 or Xmin > self.__X0:
            raise AttributeError('Incorrect value of Xmin')
        else:
            self.__Xmin = Xmin

    @property
    def Xmax(self):
        return self.__Xmax

    @Xmax.setter
    def Xmax(self, Xmax):
        if Xmax < self.__Xmin or Xmax < self.__X0:
            raise AttributeError('Incorrect value of Xmax')
        else:
            self.__Xmax = Xmax


    @property
    def dX(self):
        return self.__dX


    @dX.setter
    def dX(self, dX):
        if dX < 0:
            raise AttributeError('Positive value of dX expected')
        else:
            self.__dX = dX

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

    def setProbabilities(self,method):
        if method == "EFD":
            return self.setProbabilitiesforEFD()
        elif method == "IFD":
            return self.setProbabilitiesforIFD()
        elif method == "CNFD":
            return self.setProbabilitiesforCNFD()
        else:
            raise Exception('Unknown method.')

    def setProbabilitiesforEFD(self):
        if not self.isStockTreeCreated():
            raise Exception('Required params not set. Initialize correctly')
        t1 = (self.__sigma / self.__dX) ** 2
        t2 = (self.__r - self.__sigma ** 2 / 2) / self.__dX
        self.__pu = self.__dT / 2 * (t1 + t2)
        self.__pm = 1 - self.__dT * t1 - self.__dT * self.__r
        self.__pd = self.__dT / 2 * (t1 - t2)

    def setProbabilitiesforIFD(self):
        if not self.isStockTreeCreated():
            raise Exception('Required params not set. Initialize correctly')
        t1 = (self.__sigma / self.__dX) ** 2
        t2 = (self.__r - self.__sigma ** 2 / 2) / self.__dX
        self.__pu = -self.__dT / 2 * (t1 + t2)
        self.__pm = 1 +  self.__dT * t1  + self.__dT * self.__r
        self.__pd = -self.__dT / 2 * (t1 - t2)

    def setProbabilitiesforCNFD(self):
        if not self.isStockTreeCreated():
            raise Exception('Required params not set. Initialize correctly')
        t1 = (self.__sigma / self.__dX) ** 2
        t2 = (self.__r - self.__sigma ** 2 / 2) / self.__dX
        self.__pu = -self.__dT / 4 * (t1 + t2)
        self.__pm = 1 + self.__dT/2 * t1 + self.__dT/2 * self.__r
        self.__pd = -self.__dT / 4 * (t1 - t2)

    def getOptionPrice(self,K,payoff,method):
        if method == "EFD":
            return self.getOptionPricenyEFD(K,payoff)
        elif method == "IFD":
            return self.getOptionPricebyIFD(K,payoff)
        elif method == "CNFD":
            return self.getOptionPricebyCNFD(K,payoff)
        else:
            raise Exception('Unknown method.')

    def getOptionPricenyEFD(self,K,payOff):
        num = 0
        ndiv = (int)(self.__T / self.__dT)
        N = len(self.__S)
        A = np.zeros((N, N), dtype=np.float64)

        A[0,0] = self.__pu
        A[0,1] = self.__pm
        A[0,2] = self.__pd
        A[(N-1), (N-3)] = self.__pu
        A[(N-1), (N-2)] = self.__pm
        A[(N-1), (N-1)] = self.__pd

        for i, v in enumerate((self.__pu, self.__pm, self.__pd)):
            np.fill_diagonal(A[1:(N-1), i:], v)

        B = np.zeros(N)
        F = payOff(self.__S,K)[::-1] # reverse the array
        B[(N-1)] = self.__S[0] - self.__S[1]
        while num <= ndiv:
            F = np.dot(A,F) + B
            num += 1
        return F[::-1]


    def getOptionPricebyIFD(self, K, payOff):
        num = 0
        ndiv = (int)(self.__T / self.__dT)
        N = len(self.__S)
        A = np.zeros((N, N), dtype=np.float64)

        A[0, 0] = 1
        A[0, 1] = -1
        A[(N - 1), (N - 2)] = 1
        A[(N - 1), (N - 1)] = -1
        '''
        A[(N-1), 0] = 1
        A[(N-1), 1] = -1
        A[0, (N - 2)] = 1
        A[0, (N - 1)] = -1
        '''
        for i, v in enumerate((self.__pu, self.__pm, self.__pd)):
            np.fill_diagonal(A[1:(N - 1), i:], v)


        #B = np.zeros(N)
        B = payOff(self.__S, K)[::-1]
        B[0] = 0
        B[(N-1)] = self.__S[0] - self.__S[1]


        while num <= ndiv:
            F = np.dot(np.linalg.inv(A),B)
            B = F
            B[0] = 0
            B[(N-1)] =self.__S[0] - self.__S[1]
            num += 1
        return F[::-1]

    def getOptionPricebyCNFD(self, K, payOff):
        num = 0
        ndiv = (int)(self.__T / self.__dT)
        N = len(self.__S)
        A = np.zeros((N, N), dtype=np.float64)
        A[0, 0] = 1
        A[0, 1] = -1

        for i, v in enumerate((self.__pu, self.__pm, self.__pd)):
            np.fill_diagonal(A[1:(N - 1), i:], v)
        A[(N - 1), (N - 2)] = 1
        A[(N - 1), (N - 1)] = -1

        Z = np.zeros(N)
        '''
        B = payOff(self.__S,K)
        for i in range(1,N-1):
            Z[i]  = -self.__pu*B[N- i] - (self.__pm-2)*B[N-i-1] - self.__pd*B[N-i-2]
        Z[0] = self.__S[(N - 1)] - self.__S[(N - 2)]
        '''
        F = payOff(self.__S, K)[::-1]
        while num <= ndiv:
            for i in range(1, N - 1):
                Z[i] = -self.__pu * F[i-1] - (self.__pm - 2) * F[i] - self.__pd * F[i + 1]
                #Z[i] = -self.__pu * F[N - i] - (self.__pm - 2) * F[N - i - 1] - self.__pd * F[N - i -2]
            Z[0] = 0
            Z[-1]= self.__S[0] - self.__S[1]
            F = np.dot(np.linalg.inv(A), Z)
            num += 1
        return F[::-1]


class FiniteDifferenceMethod:
    def __init__(self, S0, r, sigma, T, dT):
        self.S0 = S0
        self.__r = r
        self.__sigma = sigma
        self.T = T
        self.dT = dT

    def initStock(self, Smin, Smax, dS):
        self.Smin = Smin
        self.Smax = Smax
        self.dS = dS
        N1 = math.ceil((Smax - self.__S0)/dS)
        N2 = math.ceil((self.__S0 - Smin)/dS)

        self.__S = np.zeros(N1 + N2 + 1)
        for i in range(N1 + 1):
            self.__S[N2 + i] = self.__S0 + i * dS

        for i in range(1, N2 + 1):
            self.__S[N2 - i] = self.__S0 - i * dS
        return self.__S, N2

    def isStockTreeCreated(self):
        params = [self.__Smin, self.__Smax, self.__S0, self.__dS, self.__S]
        return all(v is not None for v in params)

    def isProbSet(self):
        params = [self.__pu, self.__pm, self.__pd]
        return (all(v is not None for v in params) & ((self.__pu + self.__pm + self.__pd) == 1))

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
    def Smin(self):
        return self.__Smin

    @Smin.setter
    def Smin(self, Smin):
        if Smin < 0 or Smin > self.__S0:
            raise AttributeError('Incorrect value of Smin')
        else:
            self.__Smin = Smin

    @property
    def Smax(self):
        return self.__Smax

    @Smax.setter
    def Smax(self, Smax):
        if Smax < self.__Smin or Smax < self.__S0:
            raise AttributeError('Incorrect value of Smax')
        else:
            self.__Smax = Smax

    @property
    def dS(self):
        return self.__dS

    @dS.setter
    def dS(self, dS):
        if dS < 0:
            raise AttributeError('Positive value of dS expected')
        else:
            self.__dS = dS

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

    def getOptionPrice(self, K, payOff,alpha, type):
        '''
        if alpha != 1 & alpha!=0.5 & alpha !=0:
            raise Exception('Unknown method.')
        '''
        num = 0
        ndiv = (int)(self.__T / self.__dT)
        N = len(self.__S)
        A = np.zeros((N, N), dtype=np.float64)
        B = np.zeros((N, N), dtype=np.float64)

        for j in range(N-2):
            for i in range(3):
                if i ==0:
                    A[j+1][i+j] = (self.__sigma**2*((j+1)**2) - self.__r*(j+1))*(1-alpha)/2
                    if j!=0:
                        B[j][i+j-1] = -1*(self.__sigma**2*(j**2) - self.__r*j)*alpha/2
                elif i==1:
                    A[j+1][i+j] = -1/ self.__dT - (self.__sigma**2*((j+1)**2) + self.__r)*(1-alpha)
                    if j!=0:
                        B[j][i+j-1] = -1* (1/ self.__dT - (self.__sigma**2*(j**2) + self.__r)*alpha)
                elif i == 2:
                    A[j+1][i+j] = (self.__sigma**2*((j+1)**2) + self.__r*(j+1))*(1-alpha)/2
                    if j!=0:
                        B[j][i+j-1] = -1*(self.__sigma**2*(j**2) + self.__r*j)*alpha/2

        A[0,0] = 1
        A[0,1] = -1
        A[0,2] = 0
        A[(N-1),(N-1)] = -1
        A[(N-1),(N-2)] = 1
        '''
        B[0,:] = 0
        B[(N-1),:] = 0
        '''
        Intrinsic = payOff(self.__S,K)
        F = payOff(self.__S, K)
        while num <= ndiv:
            D = np.dot(B, F)
            if type == "P":
                D[-1] = 0
                D[0] = self.__S[1] - self.__S[0]
            elif type == "C":
                D[0] = 0
                D[-1] = self.__S[-1] - self.__S[-2]

            F = np.maximum(np.dot(np.linalg.inv(A),D), Intrinsic)
            #F = np.dot(np.linalg.inv(A), D)
            num += 1
        return F


def EuropeanCallPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(St-X,array_zeros)

class CIR_ImplicitFiniteDifferenceMethod:
    def __init__(self, r0, rbar, kappa, sigma, T, dT):
        self.__r0 = r0
        self.__rbar = rbar
        self.__kappa = kappa
        self.__sigma = sigma
        self.T = T
        self.dT = dT

    def initRate(self, rmin, rmax, dr):
        self.__rmin = rmin
        self.__rmax = rmax
        self.__dr = dr
        N1 = math.ceil((rmax - self.__r0)/dr)
        N2 = math.ceil((self.__r0 - rmin)/dr)

        self.__r = np.zeros(N1 + N2 + 1)
        for i in range(N1 + 1):
            self.__r[N2 + i] = self.__r0 + i * dr

        for i in range(1, N2 + 1):
            self.__r[N2 - i] = self.__r0 - i * dr
        return self.__r, N2

    def isRateTreeCreated(self):
        params = [self.__rmin, self.__rmax, self.__r0, self.__dr, self.__r]
        return all(v is not None for v in params)

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

    def getCallOptionPrice(self, BP,K):
        num = 0
        ndiv = (int)(self.__T / self.__dT)
        N = len(self.__r)
        A = np.zeros((N, N), dtype=np.float64)
        B = np.zeros((N, N), dtype=np.float64)

        for j in range(N-2):
            for i in range(3):
                if i ==0:
                    A[j+1][i+j] = - self.__kappa* self.__rbar/(2*self.__dr) + self.__kappa*(j+1)/2
                    if j!=0:
                        B[j][i+j-1] = -1*(self.__sigma**2)*j/(2*self.__dr)
                        #B[j][i+j-1] = -1*((self.__sigma*j)**2)/(2*self.__dr**2)
                elif i==1:
                    A[j+1][i+j] = - self.__dr*(j+1) - 1/self.__dT
                    if j!=0:
                        B[j][i+j-1] = -1/self.__dT + (self.__sigma**2)*j/self.__dr
                        #B[j][i+j-1] = -1/self.__dT + ((self.__sigma*j)**2)/(self.__dr**2)
                elif i == 2:
                    A[j+1][i+j] = self.__kappa* self.__rbar/(2*self.__dr) - self.__kappa*(j+1)/2
                    if j!=0:
                        B[j][i+j-1] = -1*(self.__sigma**2)*j/(2*self.__dr)
                        #B[j][i+j-1] = -1*((self.__sigma*j)**2)/(2*self.__dr**2)

        A[0,0] = -1 # was 1 earlier
        A[0,1] = 1   # was -1 earlier
        A[0,2] = 0
        A[(N-1),(N-1)] = 1 # was -1
        A[(N-1),(N-2)] = -1  # was 1
        '''
        B[0,:] = 0
        B[(N-1),:] = 0
        '''
        #Intrinsic = EuropeanCallPayOff(BP,K)
        F = EuropeanCallPayOff(BP, K)
        while num <= ndiv:
            D = np.dot(B, F)

            D[-1] = 0
            D[0] = BP[1] - BP[0]
            '''
            D[0] = 0
            D[-1] = BP[-1] - BP[-2]
            '''
            F = np.dot(np.linalg.inv(A),D)
            #F = np.dot(np.linalg.inv(A), D)
            num += 1
        return F