import numpy as np

def MbbDomain(Demand, Arg):
    BdBox = [0, 3, 0, 1]
    
    def DistFnc(P, BdBox):
        return dRectangle(P, BdBox[0], BdBox[1], BdBox[2], BdBox[3])

    def BndryCnds(Node, Element, BdBox):
        eps = 0.1 * np.sqrt((BdBox[1] - BdBox[0]) * (BdBox[3] - BdBox[2]) / Node.shape[0])
        LeftEdgeNodes = np.where(np.abs(Node[:, 0] - BdBox[0]) < eps)[0]
        LeftUpperNode = np.where((np.abs(Node[:, 0] - BdBox[0]) < eps) & (np.abs(Node[:, 1] - BdBox[3]) < eps))[0]
        RigthBottomNode = np.where((np.abs(Node[:, 0] - BdBox[1]) < eps) & (np.abs(Node[:, 1] - BdBox[2]) < eps))[0]
        FixedNodes = np.concatenate([LeftEdgeNodes, RigthBottomNode])
        Supp = np.zeros((len(FixedNodes), 3))
        Supp[:, 0] = FixedNodes
        Supp[:-1, 1] = 1
        Supp[-1, 2] = 1
        Load = np.array([[LeftUpperNode[0], 0, -0.5]])
        return [Supp, Load]

    def FixedPoints(BdBox):
        return BdBox

    switch_dict = {
        'Dist': lambda x: DistFnc(x, BdBox),
        'BC': lambda x, y, z: BndryCnds(x, y, z),
        'BdBox': BdBox,
        'PFix': lambda x: FixedPoints(x)
    }

    return switch_dict

# Rectangle distance function
def dRectangle(P, x1, x2, y1, y2):
    d = np.array([x1 - P[:, 0], P[:, 0] - x2, y1 - P[:, 1], P[:, 1] - y2]).T
    d = np.column_stack((d, np.max(d, axis=1)))
    return d