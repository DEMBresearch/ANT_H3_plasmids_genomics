class Edge:
    """This is the Edge class."""

    def __init__(self, ID, A, B, line, attNames):
        self.ID = 'd%i' % ID
        self.A = A
        self.B = B
        self.att = {attName: att for attName, att in zip(attNames, line)}

    # Retrieve node's parameters
    def getID(self, ):
        return self.ID

    def getA(self, ):
        return self.A

    def getB(self, ):
        return self.B

    def getPair(self, ):
        return (self.A, self.B)

    def getAtts(self, ):
        return self.att

    def writeGraphlm(self, keys, edgeAttNames):
        record = '<edge id="%s" source="%s" target="%s">\n' % (self.getID(), self.getA(), self.getB())
        record += ''.join(['<data key="d%i">%s</data>\n' %
                           (key, self.att[att]) for key, att in zip(keys, edgeAttNames)])
        record += '</edge>\n'
        return record
    # Change node's parameters
