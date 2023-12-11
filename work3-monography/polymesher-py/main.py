import sys
sys.path.append('./')

from polymesher import PolyMesher
from mbb import MbbDomain

domain = MbbDomain(1, 1)
poly = PolyMesher(domain, 100, 100)
print(poly)