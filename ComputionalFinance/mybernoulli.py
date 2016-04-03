__author__ = 'ashwin'
import statistics as stats
import numpy as np
import myrandgen as rndgen


class BernoulliDistribution:
    def __init__(self, p, x):
        self.p = p
        self.x = x

    @property
    def p(self):
        return self.__p

    @p.setter
    def p(self, p):
        if p.size == 0 :
            raise AttributeError('empty values of p in Bernoulli')
        else:
            self.__p = p

    @property
    def x(self):
        return self.__x

    @x.setter
    def x(self, x):
        if x.size == 0:
            raise AttributeError('empty values of X in Bernoulli')
        elif x.size != self.__p.size :
            raise AttributeError('size of X should equal to p in Bernoulli')
        else:
            self.__x = x

    def generateBernoulliDistribution(self, n):
        lgm = rndgen.LGMRandomGenerator(14)
        lgm.generateRdmNumber(n)
        l1= self.__p.size
        U = lgm.getRdmNumber()
        # try and catch
        if U.size != n:
            raise Exception

        cummprob = np.cumsum(self.__p)
        cummprob = np.insert(cummprob, 0, 0)
        self.bernoulli = np.empty(n, dtype= self.__x.dtype)
        for index in range(0, n):
            for j in range(0, l1):
                if cummprob[j] <= U[index] and U[index] <= cummprob[j+1] :
                    self.bernoulli[index] = self.__x[j]
                    break


    def getBernoulDistNumbers(self):
        return self.bernoulli


    def meanRdmNumber(self):
        """Return the sample arithmetic mean of data."""
        return stats.mean(self.bernoulli)

    def stdDeviationRdmNumber(self):
        """Calculates the population standard deviation."""
        return stats.stdev(self.bernoulli)
