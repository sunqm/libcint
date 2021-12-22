from pyscf import gto

# Disable basis normalization otherwise the basis contraction coefficients are
# scaled to meet basis normalization
gto.mole.NORMALIZE_GTO = False

mol = gto.M(
    atom='''H 0 0 -0.8; H  0 0 0.8''',
    unit='Bohr',
    basis='''
H S
6.   0.7  0.4
2.   0.6  0.3
0.8  0.5  0.2
H P
0.9  1.
''',
    cart=True)

print('contraction coefficients of shell 0:')
print(mol.bas_ctr_coeff(0))
print('contraction coefficients of shell 1:')
print(mol.bas_ctr_coeff(1))

print('integrals int1e_ipnuc_cart for shells (0, 1)')
shells = (0, 1)
print(mol.intor_by_shell('int1e_ipnuc_cart', shells))

print('integrals int2e_cart for shells (0, 1, 2, 2)')
shells = (0, 1, 2, 2)
print(mol.intor_by_shell('int2e_cart', shells))
