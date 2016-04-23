__author__ = 'ashwin'
import sys
import statistics as stats
import math
import numpy as np
MAX_INT64 = sys.maxsize
def getCumulativeDistribution(xi, X):
    count = 0
    for x in X:
        if x < xi:
            count = count + 1
    return count / X.size

def getSurvivor(xi, X):
    return 1- getCumulativeDistribution(xi, X)


def meanfn(x):
    """Return the sample arithmetic mean of data."""
    return stats.mean(x)


def stdDeviationfn(x):
    """Calculates the population standard deviation."""
    return stats.stdev(x)

def cov(a,b):
    if len(a) != len(b):
        return

    a_mean = np.mean(a)
    b_mean = np.mean(b)

    sum = 0

    for i in range(0, len(a)):
        sum += ((a[i] - a_mean) * (b[i] - b_mean))

    return sum / (len(a) - 1)


def mycuberoot(x):
    return np.piecewise(x, [x >= 0, x < 0], [lambda x: 1.*x**(1/3),lambda x: -1.*((-x)**(1/3))])


'''

def meanfn(x):
    n = 0
    Sum = 0.0
    for v in x:
        Sum += v
        n += 1
    return Sum / n

def stdDeviationfn(x):
    mu = meanfn(x)
    n = 0
    Variance = 0.0
    for v in x:
        Variance += (v - mu)*(v -mu)
        n +=1
    return math.sqrt(Variance/(n-1))
'''