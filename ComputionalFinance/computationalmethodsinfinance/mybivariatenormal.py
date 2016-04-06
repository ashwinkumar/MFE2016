import mynormal as nl
import math

class BiVariateNormalDistribution:
    def __init__(self, p):
        self.__p = p
        self.__normal = nl.NormalDistribution("Polar-Marsaglia")

    def generateBiVariateRdmNum(self, n, normaldst):
        iterator = iter(normaldst)
        num =0
        while num < n:
            n1 = iterator.__next__()
            n2 = iterator.__next__()
            x1, y1 = n1 , self.__p*n1 + math.sqrt(1 - self.__p*self.__p)* n2
            yield x1
        num += 1


    def generateDistribution(self, n):
        self.__normal.generateNormalDistribution(n)
        normaldst = self.__normal.getRdmNumbers()
        self.generateBiVariateRdmNum(n,normaldst)


    def getRdmNumbers(self):
        return self.bivariate
