__author__ = 'ashwin'
import sys
import statistics as stats
import math
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