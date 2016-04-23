import mynormal as nl
import math
import statistics as st

class BiVariateNormalDistribution:
    def __init__(self, p=-0.7):
        self.__p = p
        self.__normal = nl.NormalDistribution("Polar-Marsaglia")

    def generateBiVariateRdmNum(self, n, normaldst):
            done = object()
            n1 = next(normaldst, done)
            n2 = next(normaldst, done)
            while (n1 is not done) and (n2 is not done):
                y1 =  self.__p * n1 + math.sqrt(1 - self.__p * self.__p) * n2
                #y2 = self.__p * n2 + math.sqrt(1 - self.__p * self.__p) * n1
                yield n1,y1
                n1 = next(normaldst, done)
                n2 = next(normaldst, done)

    def setrho(self, p):
        self.__p = p
    def generateDistribution(self, n):
        self.__normal.generateNormalDistribution(n)
        normaldst = self.__normal.getRdmNumbers()
        self.bivariate = self.generateBiVariateRdmNum(n, normaldst)

    def getRdmNumbers(self):
        return self.bivariate

    def getCorr(self):
        Z = list(self.bivariate)
        X1 = [item[0] for item in Z]
        Y1 = [item[1] for item in Z]
        mu_X1 = st.mean(X1)
        mu_Y1 = st.mean(Y1)
        if( len(X1) != len(Y1)):
            raise Exception
        l = len(X1)
        nr = 0
        dr1 = 0
        dr2 = 0
        for i in range(0,l):
            nr = nr + (X1[i] - mu_X1)*( Y1[i] - mu_Y1)
            dr1 = dr1 + (X1[i] - mu_X1)* (X1[i] - mu_X1)
            dr2 = dr2 + (Y1[i] - mu_Y1) * (Y1[i] - mu_Y1)

        nr = nr/ (l-1)
        dr1 = math.sqrt(dr1/ (l-1))
        dr2 = math.sqrt(dr2/ (l-1))
        dr = dr1*dr2
        if dr == 0:
            raise Exception
        else:
            return nr/dr


