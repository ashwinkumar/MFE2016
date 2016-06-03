import numpy as np
import math
import InterestRateModel as IRM
from scipy.optimize import fsolve



class NumerixPrepayment:
    def __init__(self, PV0, WAC, rt,T):
        self.__PV0 = PV0
        self.__r = WAC/12
        self.__N = T*12

        self.__rt  = rt
        self.__ndiv = rt.shape[1]
        self.__nsims = rt.shape[0]
        self.__dT = 1/360

        self.__PV = np.zeros([self.__nsims,self.__N+1])
        self.__PV[:,0] = PV0
        self.__c = np.zeros([self.__nsims, self.__N+1])
        self.__IP = np.zeros([self.__nsims, self.__N+1])
        self.__TPP = np.zeros([self.__nsims, self.__N+1])

    def setInterestRateModel(self,IR):
        self.__IR = IR

    def SP(self,t):  ## t in months
        disc= 1/np.power(1+self.__r,self.__N - (t -1))
        return self.__PV[:,t-1]*self.__r*(1/(1- disc) -1)

    def PP(self,t):
        return  (self.__PV[:,t-1] - self.SP(t))*(1 - np.power(1- self.CPR(t),1/12))

    def IP(self,t):
        return self.__PV[:,t-1]*self.__r

    def setPV(self,t):
        self.__PV[:,t] = self.__PV[:,t-1] - (self.SP(t) + self.PP(t))

    '''
    def r_10yr(self,t): ## t in months

        #ind = (int)(((t-1))/self.__N)*self.__ndiv  # t is in years. t-1 corresponds to one month before
        ind = (t-1)*30
        step = 10*360
        if (ind+step)> self.__rt.shape[1]:
            raise AttributeError('interest rate model cannot be extended for 10 more years from given t')
        return -np.sum(self.__rt[:,ind:(ind+step)],1)/3600
    '''

    def r_10yr(self,t):## t in months
        #ind = (int)(((t - 1)) / self.__N) * self.__ndiv  # t is in years. t-1 corresponds to one month before
        ind = (t-1)*30
        return -np.log(self.__IR.BondPrice_ClosedForm(self.__rt[:,ind],t/12,t/12+10))/10
        #return -np.log(self.__IR.BondPrice_ClosedForm(self.__rt[:, t-1], t / 12, t / 12 + 10)) / 10

    def Ref_Incentive(self,t): ## t in months
        R = 12*self.__r
        return 0.28 + 0.14 * np.arctan(-8.57 + 430 * (R - self.r_10yr(t)))

    def Burnout(self,t): ## t in months
        return 0.3 + 0.7*self.__PV[:,t-1]/self.__PV0

    def Seasoning(self,t): ## t in months
        return min(1,t/30)

    def Seasonality(self,t): ## t in months
        SY = [0.98, 0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.10, 1.18, 1.22,1.23]  # SY[0] - corresponds to Dec
        return SY[t % 12]

    def CPR(self,t):
        return self.Ref_Incentive(t) * self.Burnout(t)*self.Seasoning(t)*self.Seasonality(t)

    def setc(self,t):
        self.__TPP[:,t] = self.SP(t) + self.PP(t)
        self.__IP[:,t] = self.IP(t)
        #self.__c[:,t] =self.SP(t) + self.PP(t) + self.IP(t)
        self.__c[:, t] = self.__TPP[:,t] + self.__IP[:,t]

    def generateC(self):
        for i in range(self.__N):
            self.setPV(i+1)
            self.setc(i+1)

    def setdt(self,x):
        '''
        rtx = self.__rt[:,0:(self.__N+1)] + x
        #self.__d= np.exp(-np.true_divide(rtx.cumsum(1), np.arange(1, rtx.shape[1] + 1)))
        self.__d = np.exp(-rtx.cumsum(1)*self.__dT)
        '''
        self.__d = np.zeros([self.__nsims, self.__N+1])
        rtx = self.__rt +x
        nskip = 30
        for i in range(self.__N):
            #self.__d[:, i+1] = np.exp(-np.mean(self.__rt[:, 0:(nskip * (i+1))], 1)*(i+1)/nskip)
            self.__d[:, i+1] = np.exp(-np.sum(rtx[:, 0:(nskip * (i+1))], 1)/360)


    def MBSPrice(self):
        self.setdt(0)
        self.generateC()
        P= np.sum(self.__c*self.__d,1)
        #print(np.mean(P))
        return np.mean(P)

    def MBS_IOPrice(self):
        self.setdt(0)
        self.generateC()
        P = np.sum(self.__IP * self.__d, 1)
        # print(np.mean(P))
        return np.mean(P)

    def MBS_POPrice(self):
        self.setdt(0)
        self.generateC()
        P = np.sum(self.__TPP * self.__d, 1)
        # print(np.mean(P))
        return np.mean(P)



    def solveSpread(self,x,MBSPrice):
        '''
        dx = np.zeros([self.__nsims,self.__N + 1])
        nskip = (int)(365 / 12)
        rtx = self.__rt + x
        for j in range(self.__N):
            dx[:,j + 1] = np.exp(-np.mean(rtx[:,0:(nskip * (j + 1))]) * (j + 1) / 12)
        '''
        self.setdt(x)
        P = np.mean(np.sum(self.__c * self.__d,1))
        return P - MBSPrice

    def OASSpread(self,MBSPrice):
        self.generateC()
        x = fsolve(self.solveSpread,-0.013,MBSPrice)
        return x

    def getDurationConvexity(self, ybps, MBSPrice):
        x = self.OASSpread(MBSPrice)
        P1 = self.solveSpread(x+ybps/10000,MBSPrice) + MBSPrice
        P2 = self.solveSpread(x-ybps/10000,MBSPrice) + MBSPrice
        duration = 10000 * (P2 - P1) / (2 * ybps * MBSPrice)
        convexity = math.pow(10, 8) * (P2 + P1 - 2 * MBSPrice) / (2 * ybps * ybps * MBSPrice)
        return duration,convexity

        '''
        nskip = (int)(365 / 12)
        dxup = np.zeros([self.__nsims,self.__N + 1])
        dxdwn = np.zeros([self.__nsims,self.__N + 1])
        rtup = self.__rt + x + ybps/10000
        rtdwn = self.__rt + x - ybps/10000
        for j in range(self.__N):
            #dxup[j + 1] = np.exp(-np.mean(rtup[0:(nskip * (j + 1))]) * (j + 1) / nskip)
            #dxdwn[j + 1] = np.exp(-np.mean(rtdwn[0:(nskip * (j + 1))]) * (j + 1) / nskip)
            dxup[:,j + 1] = np.exp(-np.mean(rtup[:,0:(nskip * (j + 1))]) * (j + 1) / 12)
            dxdwn[:,j + 1] = np.exp(-np.mean(rtdwn[:,0:(nskip * (j + 1))]) * (j + 1) / 12)

        P1 = np.sum(self.__c * dxup,1)
        P2 = np.sum(self.__c * dxdwn,1)
        duration=10000*(P2 - P1)/(2*ybps*MBSPrice)
        convexity=math.pow(10,8)*(P2 + P1- 2*MBSPrice)/(2*ybps*ybps*MBSPrice)
        '''

class PSA:
    def __init__(self, PV0, WAC, rt,T):
        self.__PV0 = PV0
        self.__r = WAC/12
        self.__N = T*12

        self.__rt  = rt
        self.__ndiv = rt.shape[1]
        self.__nsims = rt.shape[0]
        self.__dT = 1/360

        self.__PV = np.zeros(self.__N+1)
        self.__PV[0] = PV0
        self.__c = np.zeros(self.__N+1)

    def SP(self,t):  ## t in months
        disc= 1/np.power(1+self.__r,self.__N - (t -1))
        return self.__PV[t-1]*self.__r*(1/(1- disc) -1)

    def PP(self,t):
        return  (self.__PV[t-1] - self.SP(t))*(1 - np.power(1- self.CPR(t),1/12))

    def IP(self,t):
        return self.__PV[t-1]*self.__r

    def setPV(self,t):
        self.__PV[t] = self.__PV[t-1] - (self.SP(t) + self.PP(t))


    def CPR(self,t):
        if t<=30:
            return 0.002*t
        else:
            return 0.06

    def setc(self,t):
        self.__c[t] =self.SP(t) + self.PP(t) + self.IP(t)

    def generateC(self):
        for i in range(self.__N):
            self.setPV(i+1)
            self.setc(i+1)

    def setdt(self):
        self.__d = np.zeros([self.__nsims, self.__N+1])
        nskip = 30
        for i in range(self.__N):
            #self.__d[:, i+1] = np.exp(-np.mean(self.__rt[:, 0:(nskip * (i+1))], 1)*(i+1)/nskip)
            self.__d[:, i+1] = np.exp(-np.sum(self.__rt[:, 0:(nskip * (i+1))], 1)/360)

    def MBSPrice(self):
        self.setdt()
        self.generateC()
        P= np.sum(self.__c*self.__d,1)
        #print(np.mean(P))
        return np.mean(P)
