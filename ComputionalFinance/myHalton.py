import math

class Halton:
    def __init__(self, m):
        self.m = m
        self.__k = 1
        self.__pow = [1]
        self.__curpow = 0

    @property
    def m(self):
        return self.__m

    @m.setter
    def m(self, m):
        if m <= 0 or not isinstance(m, int):
            raise AttributeError('m is not positve integer')
        else:
            self.__m = m

    def getPow(self,i):
        if i< self.__curpow:
            return self.__pow[i]
        else:
            while i >self.__curpow:
                tmp = self.__pow[self.__curpow]*self.__m
                self.__pow.append(tmp)
                self.__curpow = self.__curpow+1

            return self.__pow[i]


    def calculateNextRdmNumber(self,k):
        sum = 0
        i=1
        while k >0:
            sum = sum + (k % self.__m)/(self.getPow(i))
            k = (int) (k /self.__m)
            i= i+1
        return sum

    def generateNextRdmNumber(self):
        rdm = self.calculateNextRdmNumber(self.__k)
        self.__k = self.__k + 1
        return rdm

    def resetSequences(self):
        self.__k = 1
        self.__curpow = 0

    def generateRdmNumberSequences(self,n, reg = True):
        num = 0
        if reg == True:
            self.resetSequences()
        while num <n:
            tmp =  self.generateNextRdmNumber()
            yield tmp
            num = num+1



