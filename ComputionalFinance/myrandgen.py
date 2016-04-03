__author__ = 'ashwin'
import statistics as stats
from array import *
import numpy as np

class RandomGenerator:
    def __init__(self, x0, a, b, m):
        self.x0 = x0
        self.a = a
        self.b = b
        self.m = m
        self.rdmarray = array('Q', [self.__x0])
        #self.normalrdmarray = array('d', [0])
        #self.rdmarray = np.array(x0, dtype=np.uint64)

    @property
    def x0(self):
        return self.__x0

    @x0.setter
    def x0(self, x0):
        if isinstance(x0, int):
            if x0 < 0:
                raise ValueError
            else:
                self.__x0 = x0
        else:
            raise TypeError

    @property
    def a(self):
        return self.__a

    @a.setter
    def a(self, a):
        if isinstance(a, int):
            if a < 0:
                raise ValueError
            else:
                self.__a = a
        else:
            raise TypeError

    @property
    def b(self):
        return self.__b

    @b.setter
    def b(self, b):
        if isinstance(b, int):
            if b < 0:
                raise ValueError
            else:
                self.__b = b
        else:
            raise TypeError

    @property
    def m(self):
        return self.__m

    @m.setter
    def m(self, m):
        if isinstance(m, int):
            if m < 0:
                raise ValueError
            else:
                self.__m = m
        else:
            raise TypeError

    def generateRdmNumber(self, n):
        flag = False
        for index in range(1, n):
            x_temp = self.__a * self.rdmarray[index - 1] + self.__b
            if not flag:
                flag = x_temp > self.__m
            if not flag:
                self.rdmarray.append(x_temp)
            else:
                self.rdmarray.append(x_temp % self.__m)
        self.normalrdmarray = np.frombuffer(self.rdmarray, dtype=np.uint64)
        self.normalrdmarray /= self.__m

    def getRdmNumber(self):
        return self.normalrdmarray

    def meanRdmNumber(self):
        """Return the sample arithmetic mean of data."""
        return stats.mean(self.normalrdmarray)

    def stdDeviationRdmNumber(self):
        """Calculates the population standard deviation."""
        return stats.stdev(self.normalrdmarray)


class LGMRandomGenerator(RandomGenerator):
    (CONST_2POW31, CONST_7POW5) = (2147483648, 16807)

    def __init__(self, x0):
        RandomGenerator.__init__(self, x0, self.CONST_7POW5, 0, 31)

    def specialMod(self, x):
        if x < self.CONST_2POW31:
            return x
        else:
            return self.specialMod(x >> 31 + (x & (self.CONST_2POW31 -1)))

    def generateRdmNumber(self, n):
        flag = False
        for index in range(1, n):
            x_temp = self.a * self.rdmarray[index - 1] + self.b
            if not flag:
                flag = x_temp > self.CONST_2POW31
            if not flag:
                self.rdmarray.append(x_temp)
            else:
                red_xtemp = ( x_temp >> 31) + (x_temp & (self.CONST_2POW31-1))
                y = self.specialMod(red_xtemp)
                self.rdmarray.append(y)

        #self.normalrdmarray[0]= self.rdmarray[0] / (self.CONST_2POW31 - 1)
        #for index in range(1, n-1):
        #    self.normalrdmarray.append(self.rdmarray[index] / (self.CONST_2POW31 - 1))
        self.normalrdmarray = np.frombuffer(self.rdmarray, dtype=np.uint64)
        self.normalrdmarray = self.normalrdmarray / (self.CONST_2POW31-1)

