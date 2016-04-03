__author__ = 'ashwin'
import statistics as stats
import numpy as np
import mybernoulli as bernoulli

class BinomialDistribution:
    def __init__(self, p, n):
        self.p = p
        self.n = n

    @property
    def p(self):
        return self.__p

    @p.setter
    def p(self, p):
        if p < 0 or p > 1:
            raise AttributeError('p values in (0,1) in Bernoulli')
        else:
            self.__p = p

    @property
    def n(self):
        return self.__n

    @p.setter
    def n(self, n):
        if n < 0:
            raise ValueError('n greater than 0')
        elif isinstance(n, int):
            self.__n = n
        else:
            raise TypeError('n positive integer')

    def generateBinomialDistribution(self, N):
        p_array = np.array((self.__p, 1- self.__p))
        X = np.array((1, 0))
        bdist = bernoulli.BernoulliDistribution(p_array, X)
        bdist.generateBernoulliDistribution(self.__n*N)
        bernouli = bdist.getBernoulDistNumbers()
        bernouli= bernouli.reshape(N, self.__n)
        self.binomdist = np.sum(bernouli, axis=1)


    def getBinomialDistribution(self):
        return self.binomdist

    def getCumulativeDistribution(self, xi):
        count= 0
        for x in self.binomdist:
            if x <xi:
                count = count+1
        return count/self.binomdist.size

    def meanRdmNumber(self):
        """Return the sample arithmetic mean of data."""
        return stats.mean(self.binomdist)

    def stdDeviationRdmNumber(self):
        """Calculates the population standard deviation."""
        return stats.stdev(self.binomdist)
