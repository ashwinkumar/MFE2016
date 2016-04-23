

class BinomialModel:
    def __init__(self, s0, u, d, r,t,n):
        self.s0= s0
        self.u = u
        self.d = d
        self.r = r
        self.t = t
        self.n = n

    @property
    def s0(self):
        return self.__s0

    @s0.setter
    def s0(self, s0):
        if s0 < 0 :
            raise AttributeError('S0 cannot be negative')
        self.__s0 = s0