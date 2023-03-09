import numpy as np
from pyscf import gto

for omega in [.3, .5, .9]:
    for l in range(6):
        for x in np.arange(1., 10., .5):
            mol = gto.M(atom=f'H 0 0 0; H 0 0 {x}',
                        basis = [[l, (0.25, 1.)]], cart=True)
            intor = "int2c2e"
            with mol.with_short_range_coulomb(omega):
                sr = mol.intor_by_shell(intor, (0,1))
            with mol.with_long_range_coulomb(omega):
                lr = mol.intor_by_shell(intor, (0,1))
            full = mol.intor_by_shell(intor, (0,1))
            ref = full - lr
            err = abs(ref - sr).max()
            if err > 1e-10 and err/abs(ref).max() > 1e-5:
                print('Failed', omega, l, x, err, abs(sr).max(), abs(ref).max())

for omega in [.3, .5, .9]:
    for l in range(6):
        for x in np.arange(1., 10., .5):
            mol = gto.M(atom=f'H 0 0 0; H 0 0 {x}; He .3 .1 .8',
                        basis = [[l, (0.25, 1.)]], cart=True)
            intor = "int3c2e"
            with mol.with_short_range_coulomb(omega):
                sr = mol.intor_by_shell(intor, (0,1,2))
            with mol.with_long_range_coulomb(omega):
                lr = mol.intor_by_shell(intor, (0,1,2))
            full = mol.intor_by_shell(intor, (0,1,2))
            ref = full - lr
            err = abs(ref - sr).max()
            if err > 1e-10 and err/abs(ref).max() > 1e-5:
                print('Failed', omega, l, x, err, abs(sr).max(), abs(ref).max())

for omega in [.3, .5, .9]:
    for l in range(6):
        for x in np.arange(1., 10., .5):
            mol = gto.M(atom=f'H 0 0 0; H 0 0 {x}; H .3 .1 .8; H .5 .5 .5',
                        basis = [[l, (0.25, 1.)]], cart=True)
            intor = "int2e"
            with mol.with_short_range_coulomb(omega):
                sr = mol.intor_by_shell(intor, (0,1,2,3))
            with mol.with_long_range_coulomb(omega):
                lr = mol.intor_by_shell(intor, (0,1,2,3))
            full = mol.intor_by_shell(intor, (0,1,2,3))
            ref = full - lr
            err = abs(ref - sr).max()
            if err > 1e-10 and err/abs(ref).max() > 1e-5:
                print('Failed', omega, l, x, err, abs(sr).max(), abs(ref).max())
