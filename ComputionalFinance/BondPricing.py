import numpy as np
import math
from scipy.optimize import fsolve
from scipy.stats import norm

class Bond:
    def __init__(self, C0,FV, T, dT,t0):
        self.T = T
        self.__t0 = t0
        self.dT = dT
        nperiod = (int)(T/dT)
        self.__nperiods = nperiod
        self.__C = C0*np.ones(nperiod)
        self.__C[nperiod-1] += FV
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

    '''
    r using vasicek model                                <r0,r1,r2,r3,r4,r5,r6>
    Bond under consideration might have coupon payments  <C0,      C0,      C0+FV>
    t - current time


    def setInterestRates(self,r, dt,t):

        if(self.__dT < dt): # the rate model dt is greater than the coupon payment
            raise AttributeError('Euler discretized rate is larger interval. Bond price will be off by large amount.')
        nskip = (int)(self.__dT/dt)
        nsims = r.shape[0]
        ind_t = (int)(t/dt)
        ind_T = (int)(t/self.__dT)
        self.__r = np.zeros([nsims,self.__nperiods])
        for i in range(self.__nperiods):
            if (i+1)< ind_T :
                self.__r[:,i] = 0
            else:
                self.__r[:,i] = np.mean(r[:,ind_t:nskip*(i+1)],1)

    def BondPrice(self,t=0):
        BP = 0
        ind_T = (int)(t/self.__dT)
        fra_T = t/self.__dT - ind_T
        if ind_T > self.__nperiods:
            raise AttributeError('t is greater than maturity T')
        for i in range(ind_T,self.__nperiods):
            if i == ind_T:
                BP = BP + self.__C[i]*np.exp(-self.__r*(i+1)*(self.__dT*(1-fra_T)))
            else:
                BP = BP + self.__C[i]*np.exp(-self.__r*(i+1)*self.__dT)
        return np.mean(BP)

    '''

    def setInterestRates(self, r, dt):

        if (self.__dT < dt):  # the rate model dt is greater than the coupon payment
            raise AttributeError('Euler discretized rate is larger interval. Bond price will be off by large amount.')
        nskip = (int)(self.__dT / dt)
        nsims = r.shape[0]
        self.__r = np.zeros([nsims, self.__nperiods])
        ind_t0 = (int)((self.__dT-self.__t0)/dt)
        #ind_t0 = (int)( self.__t0 / dt)
        #self.__r[:, 0] = np.mean(r[:,0:ind_t0], 1)

        for i in range(self.__nperiods):
            self.__r[:, i] = np.mean(r[:, 0:(ind_t0+nskip * i)], 1)

    def BondPrice(self):
        BP = 0
        for i in range(self.__nperiods):
            BP = BP + self.__C[i] * np.exp(-self.__r[:,i] * (i * self.__dT + (self.__dT - self.__t0)))
        return np.mean(BP)

    def CallPrice(self,K):
        BP = 0
        for i in range(self.__nperiods):
            BP = BP + self.__C[i] * np.exp(-self.__r[:,i] * (i * self.__dT + (self.__dT - self.__t0)))
        array_zeros = np.zeros(len(BP))
        return np.mean(np.maximum(BP - K, array_zeros))

    def PutPrice(self,K):
        BP = 0
        for i in range(self.__nperiods):
            BP = BP + self.__C[i] * np.exp(-self.__r[:,i] * (i * self.__dT + (self.__dT - self.__t0)))
        array_zeros = np.zeros(len(BP))
        return np.mean(np.maximum(K-BP, array_zeros))


    def BondPrice_Semi_ClosedForm(self,A,B,r):
        BP = 0
        for i in range(len(A)):
            BP= BP+ A[i] * np.exp(-B[i] * r) * self.__C[i]
        return BP

    def CallPrice_Semi_ClosedForm(self,A,B,r,K):
        r_star = r[:,-1]
        BP = self.BondPrice_Semi_ClosedForm(A,B,r_star)
        array_zeros = np.zeros(len(BP))
        disc = np.exp(-np.mean(r,1)*self.__t0)
        C = np.maximum(BP - K, array_zeros)
        return np.mean(C*disc)

    def get_T_array(self):
        return   np.arange(self.__dT,self.__T+ self.__dT/2, self.__dT)


    def r_star(self,r, *data):
        A, B, K = data
        Pi = A*np.exp(-B*r)
        return np.sum(self.__C*Pi)- K

    def Vasicek_CallPrice_ClosedForm(self,vasicek,K1,r0):
        T_array = self.get_T_array()
        A,B = vasicek.getA_B(self.__t0,T_array)
        data = (A,B,K1)
        r = fsolve(self.r_star,r0, args=data)
        K = A*np.exp(-B*r)
        P = vasicek.BondPrice_ClosedForm(r0,0,T_array)
        P1 = vasicek.BondPrice_ClosedForm(r0,0,self.__t0)
        sigma_p = vasicek.getsigma_p(0,self.__t0,T_array)
        d1 = np.log(P/(K*P1)) + sigma_p/2
        d2 =d1 - sigma_p
        c=  P*norm.cdf(d1) - K*P1*norm.cdf(d2)
        return np.sum(c*self.__C)

    def CIR_CallPrice_ClosedForm(self, cir, K, r0):
        chisq1 = cir.getChiSq1(0,self.__t0,self.__T,K)
        chisq2 = cir.getChiSq2(0,self.__t0, self.__T,K)
        #P_S = cir.BondPrice_ClosedForm(r0, 0, self.__T)*self.__C[-1]
        #P_T = cir.BondPrice_ClosedForm(r0, 0, self.__t0)*self.__C[-1]
        P_S = cir.BondPrice_ClosedForm(r0, 0, self.__T)*self.__C[-1]
        P_T = cir.BondPrice_ClosedForm(r0, 0, self.__t0)*self.__C[-1]

        c = P_S *chisq1 - K * P_T *chisq2
        return c

    def G2pp_CallPrice_ClosedForm(self, g2pp, K):
        P_S = g2pp.BondPrice_ClosedForm(self.__T)*self.__C[-1]
        P_T = g2pp.BondPrice_ClosedForm(self.__t0)*self.__C[-1]
        BigSigma = g2pp.BigSigma(self.__t0, self.__T)
        d1 = np.log(K*P_T/P_S)/ BigSigma - BigSigma/2
        d2 = np.log(K*P_T/P_S)/ BigSigma + BigSigma/2
        c = -P_S*norm.cdf(d1) + P_T*norm.cdf(d2)*K
        return c