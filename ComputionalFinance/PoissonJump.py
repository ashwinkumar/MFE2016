import numpy as np
import math
class PoissonJump:
    def __init__(self, gamma,T,ndiv):
        if gamma == 0:
            self.__gamma = -1
        else:
            self.__gamma = 1/gamma
        self.__T = T
        self.__ndiv = ndiv
        self.__dt  =T/ndiv

    def getJt(self):
        jdt = np.zeros(self.__ndiv+1)
        if self.__gamma == -1:
            yield jdt
        else:
            jumptimes = list(self.getJumpTimes())
            for i in range(len(jumptimes)):
                jdt[jumptimes[i]]=1
            yield jdt

    def firstJumpTime(self):
        if self.__gamma == -1 :
            return self.__ndiv+1
        else:
            u = np.random.uniform(0, 1, 1)
            num = - self.__gamma * math.log(u)
            if num <= self.__T:
                return round(float(num) / self.__dt)
            else:
                return self.__ndiv+1



    def getJumpTimes(self):
        u = np.random.uniform(0, 1, 1)
        num = - self.__gamma * math.log(u)

        while num <= self.__T:
            yield round(float(num)/ self.__dt)
            u = np.random.uniform(0, 1, 1)
            num = num - self.__gamma * math.log(u)


