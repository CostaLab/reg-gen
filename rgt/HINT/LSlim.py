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
        init()

    def init(self):
        self.localMixtureScore = [None] * self.order + 1
        self.ancestorScore = [None] * self.order + 1
        self.e = [None] * self.order + 1
        try:
            self.logGamma = ArrayHandler.clone(self.componentMixtureParameters)
            self.componentMixturePotential = ArrayHandler.clone(self.componentMixtureParameters)
            self.componentMixtureLogNorm = [None] * self.length
            self.componentMixtureIndex = [None] * self.length
            self.ancestorMixturePotential = ArrayHandler.clone(self.ancestorMixtureParameters)
            self.dependencyPotential = ArrayHandler.clone(self.dependencyParameters)
            self.dependencyLogNorm = [None] * self.length
            self.dependencyIndex = [None] * self.length
            self.ancestorMixtureLogNorm = [None] * self.length
            self.ancestorMixtureIndex = [None] * self.length

            l = 0
            while l < self.length:
                self.dependencyLogNorm[l] = [None] * self.componentMixtureParameters[l].length
                self.dependencyIndex[l] = [None] * componentMixtureParameters[l].length