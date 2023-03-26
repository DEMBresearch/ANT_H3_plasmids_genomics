from re import sub, compile

class Node:
    """This is the Node class."""

    def __init__(self, ID, name, line, attNames):
        self.ID = ID
        self.name = name
        self.att = {attName: att for attName, att in zip(attNames, line)}
        self.degree = 0
        self.indegree = 0
        self.outdegree = 0
        self.neighbours = set()
        self.X = 0.0
        self.Y = 0.0


    # Retrieve node's parameters
    def __str__(self):
        return self.ID
    
    def getID(self, ):
        return self.ID

    def getDegree(self, ):
        return self.degree

    def getIndegree(self, ):
        return self.indegree

    def getOutdegree(self, ):
        return self.outdegree
    
    def getNeighbours(self, ):
        return self.neighbours

    def getAtts(self, ):
        return self.att

    def getXY(self, ):
        return (self.X, self.Y)

    def getX(self, ):
        return self.X

    def getY(self, ):
        return self.Y

    def getAttVal(self, attribute):
        return self.att[attribute]
    
    def writeGraphlm(self, keys, nodeAttNames):
        self.att['x=float'] = str(self.getX())
        self.att['y=float'] = str(self.getY())
        self.att['degree=float'] = str(self.getDegree())
        record = '<node id="%s">\n' % self.getID()
        record += ''.join(['<data key="d%i">%s</data>\n' %
                           (key, sub(r'[<&]', 'and', self.att[att])) for key,att in zip(keys, nodeAttNames)])
        record += '</node>\n'
        
        return record

    def writeSummary(self, attNames):

        self.att['x=float'] = str(self.getX())
        self.att['y=float'] = str(self.getY())
        self.att['degree=float'] = str(self.getDegree())

        try:
            return '\t'.join([self.att[att] for att in attNames])
        except KeyError:
            for att in attNames:
                if att not in self.att:
                    print('Write summary error for ', self.name, self.att, att)
                    exit()


    # Change node's parameters
    def increaseDegree(self, ):
        self.degree += 1
        
    def addNeighbour(self, neighbour):
        self.neighbours.add(neighbour)
        
    def setDegree(self, degree): # unused
        self.degree = degree

    def setIndegree(self, indegree): # unused
        self.indegree = indegree

    def setOutdegree(self, outdegree): # unused
        self.outdegree = outdegree

    def setXY(self, coords):
        self.X = coords[0]
        self.Y = coords[1]