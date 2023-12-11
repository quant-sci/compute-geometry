import numpy as np
import matplotlib.pyplot as plt

def PolyMesher(Domain, NElem, MaxIter, P=None):
    if P is None:
        P = PolyMshr_RndPtSet(NElem, Domain)
    NElem = np.shape(P)[0]
    Tol = 5e-6
    It = 0
    Err = 1
    c = 1.5
    BdBox = Domain['BdBox']
    PFix = Domain['PFix'](BdBox)
    Area = (BdBox[1] - BdBox[0]) * (BdBox[3] - BdBox[2])
    Pc = P
    while (It <= MaxIter and Err > Tol):
        Alpha = c * np.sqrt(Area / NElem)
        P = Pc  # Lloyd's update
        R_P = PolyMshr_Rflct(P, NElem, Domain, Alpha)  # Generate the reflections
        P, R_P = PolyMshr_FixedPoints(P, R_P, PFix)  # Fixed Points 
        vor = Voronoi(np.vstack((P, R_P)))  # Construct Voronoi diagram
        Node = vor.vertices
        Element = [vor.regions[i] for i in vor.point_region[:NElem]]
        Pc, A = PolyMshr_CntrdPly(Element, Node, NElem)
        Area = np.sum(np.abs(A))
        Err = np.sqrt(np.sum((A**2) * np.sum((Pc - P)**2, axis=1))) * NElem / Area**1.5
        print(f'It: {It:3d}   Error: {Err:1.3e}')
        It = It + 1
        if NElem <= 2000:
            PolyMshr_PlotMsh(Node, Element, NElem)
    Node, Element = PolyMshr_ExtrNds(NElem, Node, Element)  # Extract node list
    Node, Element = PolyMshr_CllpsEdgs(Node, Element, 0.1)  # Remove small edges
    Node, Element = PolyMshr_RsqsNds(Node, Element)  # Reorder Nodes
    BC = Domain('BC', [Node, Element])
    Supp = BC[0]
    Load = BC[1]  # Recover BC arrays
    PolyMshr_PlotMsh(Node, Element, NElem, Supp, Load)  # Plot mesh and BCs
    return {'Node': Node, 'Element': Element, 'Supp': Supp, 'Load': Load}


def PolyMshr_RndPtSet(NElem, Domain):
    P = np.zeros((NElem, 2))
    BdBox = Domain['BdBox']
    Ctr = 0
    while Ctr < NElem:
        Y = np.zeros((NElem, 2))
        Y[:, 0] = ((BdBox[1] - BdBox[0]) * np.random.rand(NElem, 1) + BdBox[0]).reshape(-1)
        Y[:, 1] = ((BdBox[3] - BdBox[2]) * np.random.rand(NElem, 1) + BdBox[2]).reshape(-1)
        d = Domain['Dist'](Y)
        I = np.where(d[:, -1] < 0)[0]  # Index of seeds inside the domain
        NumAdded = min(NElem - Ctr, len(I))  # Number of seeds that can be added
        P[Ctr : Ctr + NumAdded, :] = Y[I[:NumAdded], :].reshape(-1, 2)
        Ctr = Ctr + NumAdded
    return P


def PolyMshr_FixedPoints(P, R_P, PFix):
    PP = np.concatenate([P, R_P], axis=0)

    for i in range(len(PFix)):
        B, I = np.argsort(np.sqrt((PP[:, 0] - PFix[i, 0])**2 + (PP[:, 1] - PFix[i, 1])**2))

        for j in range(1, 4):
            n = (PP[I[j], :] - PFix[i, :]) / np.linalg.norm(PP[I[j], :] - PFix[i, :])
            PP[I[j], :] = PP[I[j], :] - n * (B[j] - B[0])

    P = PP[:P.shape[0], :]
    R_P = PP[P.shape[0]:, :]

    return P, R_P


def PolyMshr_Rflct(P, NElem, Domain, Alpha):
    eps = 1e-8
    eta = 0.9
    d = Domain['Dist'](P)
    NBdrySegs = np.shape(d)[1] - 1
    n1 = (Domain['Dist'](P + np.tile([eps, 0], (NElem, 1))) - d) / eps
    n2 = (Domain['Dist'](P + np.tile([0, eps], (NElem, 1))) - d) / eps
    I = np.abs(d[:, :NBdrySegs]) < Alpha
    P1 = np.tile(P[:, 0], (1, NBdrySegs)).T
    P2 = np.tile(P[:, 1], (1, NBdrySegs)).T
    R_P = np.zeros((np.sum(I), 2))
    
    I = I.flatten()
    P1 = P1.flatten()
    P2 = P2.flatten()
    n1 = n1.flatten()
    n2 = n2.flatten()
    d = d.flatten()

    
    indices = np.where(I)
    R_P = np.zeros((len(indices[0]), 2))
    R_P[:, 0] = P1[indices] - 2 * n1[indices] * d[indices]
    R_P[:, 1] = P2[indices] - 2 * n2[indices] * d[indices]
    d_R_P = Domain['Dist'](R_P)
    J = np.logical_and(np.abs(d_R_P[:, -1]) >= eta * np.abs(d[indices]), d_R_P[:, -1] > 0)
    R_P = R_P[J, :]
    R_P = np.unique(R_P, axis=0)

    return R_P


def PolyMshr_CntrdPly(Element, Node, NElem):
    Pc = np.zeros((NElem, 2))
    A = np.zeros((NElem, 1))

    for el in range(NElem):
        vx = Node[Element[el], 0]
        vy = Node[Element[el], 1]
        nv = len(Element[el])
        vxS = np.roll(vx, shift=-1)
        vyS = np.roll(vy, shift=-1)
        temp = vx * vyS - vy * vxS
        A[el] = 0.5 * np.sum(temp)
        Pc[el, :] = 1 / (6 * A[el]) * [np.sum((vx + vxS) * temp), np.sum((vy + vyS) * temp)]

    return Pc, A


def PolyMshr_ExtrNds(NElem, Node0, Element0):
    map = np.unique(np.concatenate(Element0[:NElem]))
    cNode = np.arange(1, Node0.shape[0] + 1)
    cNode[np.setdiff1d(cNode, map) - 1] = np.max(map)
    Node, Element = PolyMshr_RbldLists(Node0, Element0[:NElem], cNode)
    return Node, Element


def PolyMshr_CllpsEdgs(Node0, Element0, Tol):
    while True:
        cEdge = []
        for el in range(len(Element0)):
            if len(Element0[el]) < 4:
                continue

            vx = Node0[Element0[el], 0]
            vy = Node0[Element0[el], 1]
            nv = len(vx)
            beta = np.arctan2(vy - np.sum(vy) / nv, vx - np.sum(vx) / nv)
            beta = np.mod(beta[np.roll(np.arange(nv), shift=1)] - beta, 2 * np.pi)
            betaIdeal = 2 * np.pi / len(Element0[el])
            Edge = np.column_stack((Element0[el], np.roll(Element0[el], shift=-1)))
            cEdge = np.concatenate([cEdge, Edge[beta < Tol * betaIdeal, :]])

        if len(cEdge) == 0:
            break

        cEdge = np.unique(np.sort(cEdge, axis=1), axis=0)
        cNode = np.arange(1, Node0.shape[0] + 1)

        for i in range(cEdge.shape[0]):
            cNode[cEdge[i, 1] - 1] = cNode[cEdge[i, 0] - 1]

        Node0, Element0 = PolyMshr_RbldLists(Node0, Element0, cNode)

    return Node0, Element0


def PolyMshr_RsqsNds(Node0, Element0):
    NNode0 = Node0.shape[0]
    NElem0 = len(Element0)
    ElemLnght = np.array([len(e) for e in Element0])
    nn = np.sum(ElemLnght**2)
    i = np.zeros(nn)
    j = np.zeros(nn)
    s = np.zeros(nn)
    index = 0

    for el in range(NElem0):
        eNode = Element0[el]
        ElemSet = np.arange(index + 1, index + ElemLnght[el]**2 + 1)
        i[ElemSet - 1] = np.tile(eNode, ElemLnght[el])
        j[ElemSet - 1] = np.tile(eNode, ElemLnght[el]).reshape(-1, ElemLnght[el]).T.flatten()
        s[ElemSet - 1] = 1
        index = index + ElemLnght[el]**2

    K = sparse.csr_matrix((s, (i.astype(int), j.astype(int))), shape=(NNode0, NNode0))
    p = symrcm(K)
    cNode = dict(zip(p[:NNode0], np.arange(1, NNode0 + 1)))
    Node, Element = PolyMshr_RbldLists(Node0, Element0, cNode)
    return Node, Element


def PolyMshr_RbldLists(Node0, Element0, cNode):
    Element = [np.array([cNode[node] for node in elem]) for elem in Element0]
    Node = Node0[list(cNode.keys()) - 1, :]
    return Node, Element


def PolyMshr_PlotMsh(Node, Element, NElem, Supp=None, Load=None):
    plt.clf()
    plt.axis('equal')
    plt.axis('off')
    Element = Element[:NElem]
    MaxNVer = max(map(len, Element))
    PadWNaN = lambda E: np.concatenate([E, np.full((MaxNVer - len(E),), np.nan)])
    ElemMat = np.vstack([PadWNaN(elem) for elem in Element])
    plt.fill(Node[:, 0], Node[:, 1], edgecolor='black', facecolor='none')
    plt.pause(1e-6)

    if Supp is not None:
        plt.plot(Node[Supp[:, 0] - 1, 0], Node[Supp[:, 0] - 1, 1], 'b>', markersize=8)

    if Load is not None:
        plt.plot(Node[Load[:, 0] - 1, 0], Node[Load[:, 0] - 1, 1], 'm^', markersize=8)

    plt.show()