# Cmake flag -DWITH_RANGE_COULOMB=1 should be set

import time
import numpy
import pyscf
mol = pyscf.M(
atom = '''H 0 0 0;
H 0 -1 1
H 0  1 1
H 1 0 -1''',
basis = {'H': [[0, (2., .8), (1.5, .3)],
               [0, (.625, 1)],
               [1, ( .8, 1)],
               [2, (1.2, 1)]
              ]},
)

eri = mol.intor('int2e')
with mol.with_range_coulomb(.3):
    eri_lr = mol.intor('int2e')
with mol.with_range_coulomb(-.3):
    eri_sr = mol.intor('int2e')

print(abs(eri_lr + eri_sr - eri).max())
