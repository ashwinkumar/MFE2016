__author__= 'ashwin'
import myrandgen as rndgen
import numpy as np

class NormalDistribution:
    list_methods = ["Box-Muller", "Polar-Marsaglia"]
    def __init__(self, method):
        self.method = method
        self.__lgm = rndgen.LGMRandomGenerator(14)

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
        '''self.__lgm.generateRdmNumber(n)
        U = self.__lgm.getRdmNumber().reshape(2, n / 2)
        myln = np.sqrt(-np.log(U[0]))
        mycoses = np.cos(2. * np.pi * U[1])
        mysines = np.sin(2. * np.pi * U[1])
        z1 = np.multiply(myln, mycoses) * np.sqrt(2)
        z2 = np.multiply(myln, mysines) * np.sqrt(2)
        return  np.concatenate([z1, z2])'''
        num = 0
        while num <n:
            u1 = self.__lgm.getNextRdmNumber();
            u2 = self.__lgm.getNextRdmNumber();
            yield np.sqrt(-2 * np.log(u1)) * np.cos(2. * np.pi * u2)
            yield np.sqrt(-2 * np.log(u1)) * np.sin(2. * np.pi * u2)
            num += 1

    def generatebyPolarMarsaglia(self, n):
        num = 0
        tot = 0
        while num <n:
            v1 = 2*self.__lgm.getNextRdmNumber() - 1
            v2 = 2*self.__lgm.getNextRdmNumber() - 1
            W = v1*v1 + v2*v2
            lnW = np.log(W)
            if W<=1:
                yield v1* np.sqrt(-2 * lnW / W)
                yield v2* np.sqrt(-2 * lnW / W)
                num += 1
            tot += 1
        perc = num/tot
        print(' Percentage of uniform numbers used %f' % perc)


    def generateNormalDistribution(self,n):
        if self.__method == self.list_methods[0]:
            self.normaldist = self.generatebyBoxMuller(n)
        elif self.__method == self.list_methods[1]:
            self.normaldist = self.generatebyPolarMarsaglia(n)
        else:
            raise AttributeError('unknown method for calculating normal distribution')


    def getRdmNumbers(self):
        return self.normaldist
