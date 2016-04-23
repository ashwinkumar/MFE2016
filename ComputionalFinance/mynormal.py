__author__= 'ashwin'
import myrandgen as rndgen
import numpy as np

class NormalDistribution:
    CONST_MIN = 0.00001
    list_methods = ["Box-Muller", "Polar-Marsaglia"]
    def __init__(self, method = "Polar-Marsaglia" , seed = 14):
        self.method = method
        self.__lgm = rndgen.LGMRandomGenerator(seed)

    @property
    def method(self):
        return self.__method

    @method.setter
    def method(self, method):
        if method in self.list_methods:
            self.__method = method
        else:
            raise AttributeError('unknown method for calculating normal distribution')

    def generatebyBoxMuller(self, n):
        num = 0
        while num <n:
            u1 = max(self.CONST_MIN, self.__lgm.getNextRdmNumber())
            u2 = max(self.CONST_MIN, self.__lgm.getNextRdmNumber())
            try:
                n1 = np.sqrt(-2 * np.log(u1)) * np.cos(2. * np.pi * u2)
                n2 = np.sqrt(-2 * np.log(u1)) * np.sin(2. * np.pi * u2)
            except ZeroDivisionError:
                print(u1, u2)
            yield n1
            yield n2
            num += 1

    def generatebyPolarMarsaglia(self, n):
        num = 0
        tot = 0
        while num <n:
            v1 = 2*self.__lgm.getNextRdmNumber() - 1
            v2 = 2*self.__lgm.getNextRdmNumber() - 1
            W = v1*v1 + v2*v2
            lnW = np.log(W)
            if W<=1 and W> 0 :
                yield v1* np.sqrt(-2 * lnW / W)
                yield v2* np.sqrt(-2 * lnW / W)
                num += 2
            tot += 2
        perc = num/tot
        #print(' Percentage of uniform numbers used %f' % perc)


    def generateNormalDistribution(self,n):
        if self.__method == self.list_methods[0]:
            self.normaldist = self.generatebyBoxMuller(n)
        elif self.__method == self.list_methods[1]:
            self.normaldist = self.generatebyPolarMarsaglia(n)
        else:
            raise AttributeError('unknown method for calculating normal distribution')

    def getRdmNumbers(self):
        return self.normaldist
