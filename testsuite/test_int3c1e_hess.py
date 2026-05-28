"""
Finite-difference validation of int3c1e second-derivative integrals.

Sign convention: libcint's nabla = d/dr (electron coord derivative).
Nuclear coordinate derivative: d/dR = -d/dr.
Therefore: FD of ip1 w.r.t. nuclear coord = -analytical_second_derivative.

Component ordering: for (nabla_A nabla_B), component index = alpha_A * 3 + beta_B
where alpha_A is the direction of the first nabla and beta_B is the second.
"""

import os
import numpy as np
from pyscf import gto


def _make_mol():
    """Create a test molecule with 3 atoms (distinct centers for bra/ket)."""
    return gto.M(
        atom='''
            H  0.0  0.0  0.0
            He 1.5  0.3  0.2
            Li 0.5  1.8  0.4
        ''',
        basis={'H': 'sto-3g', 'He': '6-31g', 'Li': 'sto-3g'},
        unit='Bohr',
    )


def _make_auxmol():
    """Aux mol with s, p, d, f shells on a single atom (separate from mol atoms)."""
    basis_str = {
        'He': gto.basis.parse('''
He    S
     1.5    1.0
He    P
     0.8    1.0
He    D
     1.2    1.0
He    F
     0.9    1.0
        '''),
    }
    return gto.M(
        atom='He 2.0 0.5 1.0',
        basis=basis_str,
        unit='Bohr',
    )


def _compute_3c1e(mol, auxmol, intor, comp=1):
    """Compute a 3c1e integral using mol+auxmol as supermol."""
    supermol = mol + auxmol
    nbas_mol = mol.nbas
    nbas_aux = auxmol.nbas
    shls_slice = (0, nbas_mol, 0, nbas_mol, nbas_mol, nbas_mol + nbas_aux)
    return supermol.intor(intor, comp=comp, shls_slice=shls_slice)


def test_int3c1e_ipip1():
    """
    Verify int3c1e_ipip1 via FD of int3c1e_ip1 w.r.t. bra atom.

    For shell triple (i on atom A, j NOT on A, k NOT on A):
      d/dR_A,β [ip1[α, i, j, k]] = -ipip1[α*3+β, i, j, k]
    """
    mol = _make_mol()
    auxmol = _make_auxmol()
    delta = 1e-5

    analytic = _compute_3c1e(mol, auxmol, 'int3c1e_ipip1', comp=9)
    nao = mol.nao
    naux = auxmol.nao
    assert analytic.shape == (9, nao, nao, naux)

    max_err = 0.0
    aoslice = mol.aoslice_by_atom()

    for iatm in range(mol.natm):
        i0, i1 = int(aoslice[iatm, 2]), int(aoslice[iatm, 3])
        # j indices NOT on iatm
        j_mask = np.ones(nao, dtype=bool)
        j_mask[i0:i1] = False

        for beta in range(3):
            coords_p = mol.atom_coords(unit='Bohr').copy()
            coords_m = mol.atom_coords(unit='Bohr').copy()
            coords_p[iatm, beta] += delta
            coords_m[iatm, beta] -= delta
            mol_p = mol.set_geom_(coords_p, unit='Bohr', inplace=False)
            mol_m = mol.set_geom_(coords_m, unit='Bohr', inplace=False)
            ip1_p = _compute_3c1e(mol_p, auxmol, 'int3c1e_ip1', comp=3)
            ip1_m = _compute_3c1e(mol_m, auxmol, 'int3c1e_ip1', comp=3)
            dip1 = (ip1_p - ip1_m) / (2 * delta)

            for alpha in range(3):
                comp_idx = alpha * 3 + beta
                # FD = -ipip1
                ana_block = analytic[comp_idx, i0:i1][:, j_mask, :]
                fd_block = -dip1[alpha, i0:i1][:, j_mask, :]
                err = np.max(np.abs(ana_block - fd_block))
                max_err = max(max_err, err)

    print(f"int3c1e_ipip1: max FD error = {max_err:.2e}")
    assert max_err < 1e-7, f"int3c1e_ipip1 FD error too large: {max_err}"


def test_int3c1e_ip1ip2():
    """
    Verify int3c1e_ip1ip2 by displacing aux atom.

    Since aux atoms are distinct from mol atoms:
      d/dR_C,β [ip1[α, i, j, k]] = -ip1ip2[α*3+β, i, j, k]  for k on C
    """
    mol = _make_mol()
    auxmol = _make_auxmol()
    delta = 1e-5

    analytic = _compute_3c1e(mol, auxmol, 'int3c1e_ip1ip2', comp=9)
    nao = mol.nao
    naux = auxmol.nao
    assert analytic.shape == (9, nao, nao, naux)

    max_err = 0.0
    aoslice_aux = auxmol.aoslice_by_atom()

    for katm in range(auxmol.natm):
        k0, k1 = int(aoslice_aux[katm, 2]), int(aoslice_aux[katm, 3])
        for beta in range(3):
            coords_p = auxmol.atom_coords(unit='Bohr').copy()
            coords_m = auxmol.atom_coords(unit='Bohr').copy()
            coords_p[katm, beta] += delta
            coords_m[katm, beta] -= delta
            auxmol_p = auxmol.set_geom_(coords_p, unit='Bohr', inplace=False)
            auxmol_m = auxmol.set_geom_(coords_m, unit='Bohr', inplace=False)
            ip1_p = _compute_3c1e(mol, auxmol_p, 'int3c1e_ip1', comp=3)
            ip1_m = _compute_3c1e(mol, auxmol_m, 'int3c1e_ip1', comp=3)
            dip1 = (ip1_p - ip1_m) / (2 * delta)

            for alpha in range(3):
                comp_idx = alpha * 3 + beta
                ana_block = analytic[comp_idx, :, :, k0:k1]
                fd_block = -dip1[alpha, :, :, k0:k1]
                err = np.max(np.abs(ana_block - fd_block))
                max_err = max(max_err, err)

    print(f"int3c1e_ip1ip2: max FD error = {max_err:.2e}")
    assert max_err < 1e-7, f"int3c1e_ip1ip2 FD error too large: {max_err}"


def test_int3c1e_ipvip1():
    """
    Verify int3c1e_ipvip1 by displacing ket atom.

    For shell triple (i NOT on B, j on B, k NOT on B):
      d/dR_B,β [ip1[α, i, j, k]] = -ipvip1[α*3+β, i, j, k]
    """
    mol = _make_mol()
    auxmol = _make_auxmol()
    delta = 1e-5

    analytic = _compute_3c1e(mol, auxmol, 'int3c1e_ipvip1', comp=9)
    nao = mol.nao
    naux = auxmol.nao
    assert analytic.shape == (9, nao, nao, naux)

    max_err = 0.0
    aoslice = mol.aoslice_by_atom()

    for jatm in range(mol.natm):
        j0, j1 = int(aoslice[jatm, 2]), int(aoslice[jatm, 3])
        # i indices NOT on jatm
        i_mask = np.ones(nao, dtype=bool)
        i_mask[j0:j1] = False

        for beta in range(3):
            coords_p = mol.atom_coords(unit='Bohr').copy()
            coords_m = mol.atom_coords(unit='Bohr').copy()
            coords_p[jatm, beta] += delta
            coords_m[jatm, beta] -= delta
            mol_p = mol.set_geom_(coords_p, unit='Bohr', inplace=False)
            mol_m = mol.set_geom_(coords_m, unit='Bohr', inplace=False)
            ip1_p = _compute_3c1e(mol_p, auxmol, 'int3c1e_ip1', comp=3)
            ip1_m = _compute_3c1e(mol_m, auxmol, 'int3c1e_ip1', comp=3)
            dip1 = (ip1_p - ip1_m) / (2 * delta)

            for alpha in range(3):
                comp_idx = alpha * 3 + beta
                ana_block = analytic[comp_idx][i_mask][:, j0:j1, :]
                fd_block = -dip1[alpha][i_mask][:, j0:j1, :]
                err = np.max(np.abs(ana_block - fd_block))
                max_err = max(max_err, err)

    print(f"int3c1e_ipvip1: max FD error = {max_err:.2e}")
    assert max_err < 1e-7, f"int3c1e_ipvip1 FD error too large: {max_err}"


def test_int3c1e_ipip2():
    """
    Verify int3c1e_ipip2 (d²/dR₃²) via double FD of int3c1e.

    Since aux atoms are separate from mol, no contamination.
    Sign: ipip2 = (nabla_k)^2 S. FD gives d²S/dR_C^2.
    d²S/dR_C,α dR_C,β = (-d/dr_k,α)(-d/dr_k,β) S = ipip2[α*3+β]
    So FD should directly match ipip2.
    """
    mol = _make_mol()
    auxmol = _make_auxmol()
    delta = 1e-4

    analytic = _compute_3c1e(mol, auxmol, 'int3c1e_ipip2', comp=9)
    nao = mol.nao
    naux = auxmol.nao
    assert analytic.shape == (9, nao, nao, naux)

    max_err = 0.0
    aoslice_aux = auxmol.aoslice_by_atom()

    for katm in range(auxmol.natm):
        k0, k1 = int(aoslice_aux[katm, 2]), int(aoslice_aux[katm, 3])
        for alpha in range(3):
            for beta in range(3):
                coords_pp = auxmol.atom_coords(unit='Bohr').copy()
                coords_pm = auxmol.atom_coords(unit='Bohr').copy()
                coords_mp = auxmol.atom_coords(unit='Bohr').copy()
                coords_mm = auxmol.atom_coords(unit='Bohr').copy()
                coords_pp[katm, alpha] += delta
                coords_pp[katm, beta] += delta
                coords_pm[katm, alpha] += delta
                coords_pm[katm, beta] -= delta
                coords_mp[katm, alpha] -= delta
                coords_mp[katm, beta] += delta
                coords_mm[katm, alpha] -= delta
                coords_mm[katm, beta] -= delta

                auxmol_pp = auxmol.set_geom_(coords_pp, unit='Bohr', inplace=False)
                auxmol_pm = auxmol.set_geom_(coords_pm, unit='Bohr', inplace=False)
                auxmol_mp = auxmol.set_geom_(coords_mp, unit='Bohr', inplace=False)
                auxmol_mm = auxmol.set_geom_(coords_mm, unit='Bohr', inplace=False)

                s_pp = _compute_3c1e(mol, auxmol_pp, 'int3c1e', comp=1)
                s_pm = _compute_3c1e(mol, auxmol_pm, 'int3c1e', comp=1)
                s_mp = _compute_3c1e(mol, auxmol_mp, 'int3c1e', comp=1)
                s_mm = _compute_3c1e(mol, auxmol_mm, 'int3c1e', comp=1)

                d2s = (s_pp - s_pm - s_mp + s_mm) / (4 * delta**2)

                comp_idx = alpha * 3 + beta
                err = np.max(np.abs(
                    analytic[comp_idx, :, :, k0:k1] - d2s[:, :, k0:k1]
                ))
                max_err = max(max_err, err)

    print(f"int3c1e_ipip2: max FD error = {max_err:.2e}")
    assert max_err < 1e-6, f"int3c1e_ipip2 FD error too large: {max_err}"


def test_translational_invariance():
    """
    Verify: ipip1 + ipvip1 + ip1ip2 = 0

    From d/dR₁ (ip1) + d/dR₂ (ip1) + d/dR₃ (ip1) = 0
    where each d/dR introduces a sign, but all three nabla operators
    are electron-coordinate derivatives, so the relation holds directly.
    """
    mol = _make_mol()
    auxmol = _make_auxmol()

    ipip1 = _compute_3c1e(mol, auxmol, 'int3c1e_ipip1', comp=9)
    ip1ip2 = _compute_3c1e(mol, auxmol, 'int3c1e_ip1ip2', comp=9)
    ipvip1 = _compute_3c1e(mol, auxmol, 'int3c1e_ipvip1', comp=9)

    total = ipip1 + ipvip1 + ip1ip2
    max_err = np.max(np.abs(total))
    print(f"Translational invariance (ipip1+ipvip1+ip1ip2=0): max error = {max_err:.2e}")
    assert max_err < 1e-12, \
        f"Translational invariance violated: {max_err}"


def test_shapes():
    """Verify all new integrals return correct shapes."""
    mol = _make_mol()
    auxmol = _make_auxmol()
    nao = mol.nao
    naux = auxmol.nao

    for name in ['int3c1e_ipip1', 'int3c1e_ip1ip2',
                 'int3c1e_ipvip1', 'int3c1e_ipip2']:
        result = _compute_3c1e(mol, auxmol, name, comp=9)
        expected = (9, nao, nao, naux)
        assert result.shape == expected, \
            f"{name}: shape {result.shape} != expected {expected}"
        print(f"{name}: shape OK {result.shape}")


def test_ipip1_symmetry():
    """Verify d²/dR_A,α dR_A,β == d²/dR_A,β dR_A,α for ipip1."""
    mol = _make_mol()
    auxmol = _make_auxmol()

    result = _compute_3c1e(mol, auxmol, 'int3c1e_ipip1', comp=9)
    # comp = alpha*3+beta
    for alpha in range(3):
        for beta in range(alpha + 1, 3):
            c_ab = result[alpha * 3 + beta]
            c_ba = result[beta * 3 + alpha]
            max_diff = np.max(np.abs(c_ab - c_ba))
            print(f"ipip1 symmetry ({alpha},{beta}): max diff = {max_diff:.2e}")
            assert max_diff < 1e-12, \
                f"ipip1 not symmetric for ({alpha},{beta}): {max_diff}"


def test_ipip2_symmetry():
    """Verify d²/dR_C,α dR_C,β == d²/dR_C,β dR_C,α for ipip2."""
    mol = _make_mol()
    auxmol = _make_auxmol()

    result = _compute_3c1e(mol, auxmol, 'int3c1e_ipip2', comp=9)
    for alpha in range(3):
        for beta in range(alpha + 1, 3):
            c_ab = result[alpha * 3 + beta]
            c_ba = result[beta * 3 + alpha]
            max_diff = np.max(np.abs(c_ab - c_ba))
            print(f"ipip2 symmetry ({alpha},{beta}): max diff = {max_diff:.2e}")
            assert max_diff < 1e-12, \
                f"ipip2 not symmetric for ({alpha},{beta}): {max_diff}"


if __name__ == '__main__':
    print("\n=== Shape tests ===")
    test_shapes()

    print("\n=== Symmetry tests ===")
    test_ipip1_symmetry()
    test_ipip2_symmetry()

    print("\n=== Translational invariance ===")
    test_translational_invariance()

    print("\n=== Finite-difference tests ===")
    test_int3c1e_ipip1()
    test_int3c1e_ip1ip2()
    test_int3c1e_ipvip1()
    test_int3c1e_ipip2()

    print("\n=== All tests passed ===")
