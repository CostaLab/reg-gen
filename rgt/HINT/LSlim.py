import numpy as np


class LSlim:

    def __init__(self, length, order, distance, ess, q, t):
        self.length = length
        self.order = order
        self.distance = distance
        self.ess = ess,
        self.q = q
        self.type_ = t

        # parameters
        self.dependencyParameters = [None] * length
        self.componentMixtureParameters = [None] * length
        self.ancestorMixtureParameters = [None] * length

        self.localMixtureScore = [None] * (order + 1)

        a = (self.alphabets.getAlphabetLengthAt(0))
        l = 0
        while l < length:
            self.componentMixtureParameters[l] = [None] * (min(l, order) + 1)
            self.dependencyParameters[l] = [None] * (min(l, order) + 1)
            self.ancestorMixtureParameters[l] = [None] * (min(l, order) + 1)

            context = 1
            dist = 1
            o = 0
            while o < min(l, order) + 1:
                self.dependencyParameters[l][o] = [None] * context
                if o != 0:
                    dist = min(l - o + 1, distance)
                self.ancestorMixtureParameters[l][o] = [None] * dist
                context *= a
            o += 1
        l += 1
        self.init()

    def init(self):
        self.localMixtureScore = [None] * self.order + 1
        self.ancestorScore = [None] * self.order + 1
        self.e = [None] * self.order + 1