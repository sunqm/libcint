import numpy
from pyscf import lib
from pyscf.pbc import gto, df
import scipy.special
cell = gto.M(
atom = '''H 0 0 0;
H 0 -1 1
H 0  1 1
H 1 0 -1''',
basis = {'H': [[0, (2., 1)],
               [0, (.625, 1)],
               [1, ( .8, 1)],
               #[2, (1.2, 1)]
              ]},
dimension = 0,
a = numpy.eye(3),
gs = [30]*3
)

zeta = .4
def weighted_stg(kpt=numpy.zeros(3), exx=False, gs=None):
    if gs is None:
        gs = self.gs
    Gv, Gvbase, kws = cell.get_Gv_weights(gs)
    G2 = numpy.einsum('gx,gx->g', Gv, Gv)
    coulG = 8*numpy.pi*zeta / (G2 + zeta**2)**2
    coulG *= kws
    return coulG

def weighted_yp(kpt=numpy.zeros(3), exx=False, gs=None):
    if gs is None:
        gs = self.gs
    Gv, Gvbase, kws = cell.get_Gv_weights(gs)
    G2 = numpy.einsum('gx,gx->g', Gv, Gv)
    coulG = 4*numpy.pi / (G2 + zeta**2)
    coulG *= kws
    return coulG

def get_eri(mydf):
    eriR = 0
    kptijkl = numpy.zeros((4,3))
    q = numpy.zeros(3)
    coulG = mydf.weighted_coulG(q, False, mydf.gs)
    for pqkR, pqkI, p0, p1 \
            in mydf.pw_loop(mydf.gs, kptijkl[:2], q, aosym='s2'):
        vG = coulG[p0:p1]
        pqkRv = pqkR * vG
        pqkIv = pqkI * vG
        eriR += lib.ddot(pqkRv, pqkR.T)
        eriR += lib.ddot(pqkIv, pqkI.T)
        pqkR = pqkI = None
    return eriR

mydf = df.AFTDF(cell)

eri0 = get_eri(mydf)
eri1 = cell.intor('cint2e_sph', aosym='s4')
print abs(eri0-eri1).max(), abs(eri0).max()

mydf.weighted_coulG = weighted_stg
eri0 = get_eri(mydf)
cell.set_f12_zeta(zeta)
eri1 = cell.intor('cint2e_stg_sph', aosym='s4')
print abs(eri0-eri1).max(), abs(eri0).max()

mydf.weighted_coulG = weighted_yp
eri0 = get_eri(mydf)
cell.set_f12_zeta(zeta)
eri1 = cell.intor('cint2e_yp_sph', aosym='s4')
print abs(eri0-eri1).max(), abs(eri0).max()


#def weighted_gtg_coul(kpt=numpy.zeros(3), exx=False, gs=None):
#    if gs is None:
#        gs = self.gs
#    Gv, Gvbase, kws = cell.get_Gv_weights(gs)
#    G2 = numpy.einsum('gx,gx->g', Gv, Gv)
#    G = numpy.sqrt(G2)
#    coulG = 4*numpy.pi / (G*numpy.sqrt(zeta)) * scipy.special.dawsn(G/(2*numpy.sqrt(zeta)))
#    coulG *= kws
#    return coulG
#
#mydf.weighted_coulG = weighted_gtg_coul
#eri0 = get_eri(mydf)
