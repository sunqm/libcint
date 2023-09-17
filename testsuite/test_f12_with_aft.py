import numpy
from pyscf import lib
from pyscf.pbc import gto, df
import scipy.special
cell = gto.M(
atom = '''H 0 0 0;
H 0 -1 1
H 0  1 1
H 1 0 -1''',
basis = {'H': [[0, (2., .8), (1.5, .3)],
               [0, (.625, 1)],
               [1, ( .8, 1)],
               [2, (1.2, 1)]
              ]},
dimension = 0,
a = numpy.eye(3),
mesh = [60]*3
)
nao = cell.nao

zeta = .4
def weighted_stg(kpt=numpy.zeros(3), exx=False, mesh=None):
    if mesh is None:
        mesh = self.mesh
    Gv, Gvbase, kws = cell.get_Gv_weights(mesh)
    G2 = numpy.einsum('gx,gx->g', Gv, Gv)
    coulG = 8*numpy.pi*zeta / (G2 + zeta**2)**2
    coulG *= kws
    return coulG

def weighted_yp(kpt=numpy.zeros(3), exx=False, mesh=None):
    if mesh is None:
        mesh = self.mesh
    Gv, Gvbase, kws = cell.get_Gv_weights(mesh)
    G2 = numpy.einsum('gx,gx->g', Gv, Gv)
    coulG = 4*numpy.pi / (G2 + zeta**2)
    coulG *= kws
    return coulG

def get_eri(mydf):
    eriR = 0
    kptijkl = numpy.zeros((4,3))
    q = numpy.zeros(3)
    coulG = mydf.weighted_coulG(q, False, mydf.mesh)
    for pqkR, pqkI, p0, p1 \
            in mydf.pw_loop(mydf.mesh, kptijkl[:2], q, aosym='s2'):
        vG = coulG[p0:p1]
        pqkRv = pqkR * vG
        pqkIv = pqkI * vG
        # rho(G) v(G) rho(-G) = rho(G) v(G) [rho(G)]^*
        eriR += lib.ddot(pqkRv, pqkR.T)
        eriR += lib.ddot(pqkIv, pqkI.T)
        pqkR = pqkI = None
    return eriR

def get_eri_ip1(mydf):
    eriR = 0
    kptijkl = numpy.zeros((4,3))
    q = numpy.zeros(3)
    coulG = mydf.weighted_coulG(q, False, mydf.mesh)
    Gv, Gvbase, kws = mydf.cell.get_Gv_weights(mydf.mesh)
    for pqkR, pqkI, p0, p1 \
            in mydf.pw_loop(mydf.mesh, kptijkl[:2], q, aosym='s1'):
        vG = coulG[p0:p1]
        G_vG = Gv[p0:p1].T * vG  # = -i\nabla * f12(G)
        pqkRv = (pqkR * G_vG[:,None,:]).reshape(-1,p1-p0)
        pqkIv = (pqkI * G_vG[:,None,:]).reshape(-1,p1-p0)
        # Imaginary part of rho(G) [-i\nabla*f12(G)] [rho(G)]^*
        eriR += lib.ddot(pqkIv, pqkR.reshape(-1,p1-p0).T)
        eriR -= lib.ddot(pqkRv, pqkI.reshape(-1,p1-p0).T)
        pqkR = pqkI = None
    return eriR.reshape(3,nao,nao,nao,nao)

def get_eri_laplacian(mydf):
    eriR = 0
    kptijkl = numpy.zeros((4,3))
    q = numpy.zeros(3)
    coulG = mydf.weighted_coulG(q, False, mydf.mesh)
    Gv, Gvbase, kws = mydf.cell.get_Gv_weights(mydf.mesh)
    G2 = numpy.einsum('gx,gx->g', Gv, Gv)
    for pqkR, pqkI, p0, p1 \
            in mydf.pw_loop(mydf.mesh, kptijkl[:2], q, aosym='s1'):
        vG = coulG[p0:p1]
        G2_vG = -G2[p0:p1] * vG  # = \nabla^2 * f12(G)
        pqkRv = (pqkR * G2_vG)
        pqkIv = (pqkI * G2_vG)
        eriR += lib.ddot(pqkRv, pqkR.T)
        eriR += lib.ddot(pqkIv, pqkI.T)
        pqkR = pqkI = None
    return eriR.reshape(nao,nao,nao,nao)

def get_eri_3c2e(mydf):
    from pyscf.gto.ft_ao import ft_ao
    eriR = 0
    kptijkl = numpy.zeros((4,3))
    q = numpy.zeros(3)
    coulG = mydf.weighted_coulG(q, False, mydf.mesh)
    Gv, Gvbase, kws = cell.get_Gv_weights(mydf.mesh)
    ao = ft_ao(cell, Gv)

    for pqkR, pqkI, p0, p1 \
            in mydf.pw_loop(mydf.mesh, kptijkl[:2], q, aosym='s2'):
        vG = coulG[p0:p1]
        pqkRv = pqkR * vG
        pqkIv = pqkI * vG
        # rho(G) v(G) rho(-G) = rho(G) v(G) [rho(G)]^*
        eriR += lib.ddot(pqkRv, ao[p0:p1].real)
        eriR += lib.ddot(pqkIv, ao[p0:p1].imag)
        pqkR = pqkI = None
    return eriR

def get_eri_2c2e(mydf):
    from pyscf.gto.ft_ao import ft_ao
    eriR = 0
    kptijkl = numpy.zeros((4,3))
    q = numpy.zeros(3)
    coulG = mydf.weighted_coulG(q, False, mydf.mesh)
    Gv, Gvbase, kws = cell.get_Gv_weights(mydf.mesh)
    ao = ft_ao(cell, Gv)
    return lib.dot(ao.T*coulG, ao.conj())


mydf = df.AFTDF(cell)

# Regular ERIs
eri0 = get_eri(mydf)
eri1 = cell.intor('int2e_sph', aosym='s4')
print('int2e', abs(eri0-eri1).max(), abs(eri0).max())

eri0 = get_eri_ip1(mydf)
eri1 = cell.intor('int2e_ip1_sph', aosym='s1', comp=3)
eri1 = -eri1 - eri1.transpose(0,2,1,3,4)
print('int2e_ip1', abs(eri0-eri1).max(), abs(eri0).max())

eri0 = get_eri_laplacian(mydf)
eri1_ipip1 = cell.intor('int2e_ipip1_sph', aosym='s1', comp=9).reshape(3,3,nao,nao,nao,nao)
eri1_ipvip1 = cell.intor('int2e_ipvip1_sph', aosym='s1', comp=9).reshape(3,3,nao,nao,nao,nao)
eri1 = eri1_ipip1 + eri1_ipip1.transpose(0,1,3,2,4,5)
eri1+= eri1_ipvip1 * 2
print('laplacian', abs(eri0-eri1.trace()).max(), abs(eri0).max())

eri1_ip1ip2 = cell.intor('int2e_ip1ip2_sph', aosym='s1', comp=9).reshape(3,3,nao,nao,nao,nao)
eri1 =(eri1_ip1ip2 + eri1_ip1ip2.transpose(0,1,3,2,4,5)
       + eri1_ip1ip2.transpose(0,1,2,3,5,4) + eri1_ip1ip2.transpose(0,1,3,2,5,4))
eri1 *= -1  # (-) sign associated to particle 2 in |r1-r2|
print('laplacian', abs(eri0-eri1.trace()).max(), abs(eri0).max())


# STG
mydf.weighted_coulG = weighted_stg
cell.set_f12_zeta(zeta)
eri0 = get_eri(mydf)
eri1 = cell.intor('int2e_stg_sph', aosym='s4')
print('int2e_stg', abs(eri0-eri1).max(), abs(eri0).max())

eri0 = get_eri_ip1(mydf)
eri1 = cell.intor('int2e_stg_ip1_sph', aosym='s1', comp=3)
eri1 = -eri1 - eri1.transpose(0,2,1,3,4)
print('int2e_stg_ip1', abs(eri0-eri1).max(), abs(eri0).max())

eri0 = get_eri_laplacian(mydf)
eri1_ipip1 = cell.intor('int2e_stg_ipip1_sph', aosym='s1', comp=9).reshape(3,3,nao,nao,nao,nao)
eri1_ipvip1 = cell.intor('int2e_stg_ipvip1_sph', aosym='s1', comp=9).reshape(3,3,nao,nao,nao,nao)
eri1 = eri1_ipip1 + eri1_ipip1.transpose(0,1,3,2,4,5)
eri1+= eri1_ipvip1 * 2
print('laplacian', abs(eri0-eri1.trace()).max(), abs(eri0).max())

eri1_ip1ip2 = cell.intor('int2e_stg_ip1ip2_sph', aosym='s1', comp=9).reshape(3,3,nao,nao,nao,nao)
eri1 =(eri1_ip1ip2 + eri1_ip1ip2.transpose(0,1,3,2,4,5)
       + eri1_ip1ip2.transpose(0,1,2,3,5,4) + eri1_ip1ip2.transpose(0,1,3,2,5,4))
eri1 *= -1  # (-) sign associated to particle 2 in |r1-r2|
print('laplacian', abs(eri0-eri1.trace()).max(), abs(eri0).max())


# Yukawa potential
mydf.weighted_coulG = weighted_yp
cell.set_f12_zeta(zeta)
eri0 = get_eri(mydf)
eri1 = cell.intor('int2e_yp_sph', aosym='s4')
print('int2e_yp', abs(eri0-eri1).max(), abs(eri0).max())

eri0 = get_eri_ip1(mydf)
eri1 = cell.intor('int2e_yp_ip1_sph', aosym='s1', comp=3)
eri1 = -eri1 - eri1.transpose(0,2,1,3,4)
print('int2e_yp_ip1', abs(eri0-eri1).max(), abs(eri0).max())

eri0 = get_eri_laplacian(mydf)
eri1_ipip1 = cell.intor('int2e_yp_ipip1_sph', aosym='s1', comp=9).reshape(3,3,nao,nao,nao,nao)
eri1_ipvip1 = cell.intor('int2e_yp_ipvip1_sph', aosym='s1', comp=9).reshape(3,3,nao,nao,nao,nao)
eri1 = eri1_ipip1 + eri1_ipip1.transpose(0,1,3,2,4,5)
eri1+= eri1_ipvip1 * 2
print('laplacian', abs(eri0-eri1.trace()).max(), abs(eri0).max())

eri1_ip1ip2 = cell.intor('int2e_yp_ip1ip2_sph', aosym='s1', comp=9).reshape(3,3,nao,nao,nao,nao)
eri1 =(eri1_ip1ip2 + eri1_ip1ip2.transpose(0,1,3,2,4,5)
       + eri1_ip1ip2.transpose(0,1,2,3,5,4) + eri1_ip1ip2.transpose(0,1,3,2,5,4))
eri1 *= -1  # (-) sign associated to particle 2 in |r1-r2|
print('laplacian', abs(eri0-eri1.trace()).max(), abs(eri0).max())


omega = 0.2
def weighted_rsh(kpt=numpy.zeros(3), exx=False, mesh=None):
    if mesh is None:
        mesh = self.mesh
    Gv, Gvbase, kws = cell.get_Gv_weights(mesh)
    G2 = numpy.einsum('gx,gx->g', Gv, Gv)
    coulG = 4*numpy.pi / G2 * numpy.exp(-G2/(4*omega**2))
    coulG *= kws
    return coulG

mydf.weighted_coulG = weighted_rsh
eri0 = get_eri(mydf)
with cell.with_range_coulomb(omega):
    eri1 = cell.intor('int2e_sph', aosym='s4', comp=1)
print('int2e', abs(eri0-eri1).max(), abs(eri0).max())
