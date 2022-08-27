libcint
=======

version 5.1.6
2022-08-26


What is libcint
---------------

libcint is an open source library for analytical Gaussian integrals.
It provides C/Fortran API to evaluate one-electron / two-electron
integrals for Cartesian / real-spheric / spinor Gaussian type functions.


Features
--------

* Various GTO type:

  - Cartesian GTO:  s, p, 6d, 10f, 15g, 21h, 28i Gaussian type functions.
  - Real-spheric GTO:  s, p, 5d, 7f, 9g, 11h, 13i Gaussian type functions.
  - Spinor GTO:  J-adapted spinor Gaussian functions.

* One electron integrals.

  - Regular kinetic-like integrals.
  - Nuclear attraction-like integrals (Gaussian nuclear model are supported).

* Two electron integrals (value < 1e-15 are neglected) include

  - Coulomb repulsion
  - Gaunt interaction
  - Breit interaction
  - 2-center, 3-center and 4-center integrals
  - Long-range part and short-range part of range-separated Coulomb

* Common lisp script to generate C code for new integrals.
* Thread safe.
* Uniform API for all kind of integrals.
  - one electron integrals::

        not0 = fn1e_name(double *buf, int *atm, int natm, int *bas, int nbas, double *env);

  - two electron integrals::

        not0 = fn2e_name(double *buf, int *atm, int natm, int *bas, int nbas, double *env, NULL);

  - the return boolean (not0) gives the summary whether the integrals
    are completely 0.

* Minimal overhead of initialization

  - Pre-computation is not required.  Only basic info (see previous API)
    of basis function need to be initialized, within a plain integer or
    double precision array.  (For 2-e integral, there is an optional
    argument called optimizer which can be switched off by setting it to
    NULL.  Using optimizer should not affect the value of integral, but
    can increase the performance by ~10%.)

* Minimal dependence on external library.

  - BLAS is the only library needed.  Normally, the performance
    difference due to various BLAS implementations is less than 1%.

* Small memory usage.

  - Very few intermediate data are stored.  ~80% of the memory are
    allocated for holding the whole contracted Cartesion integrals,
    which typically should be less than 1 Mega bytes.


Getting libcint
---------------

The newest version is available on GitHub: ::

    git clone http://github.com/sunqm/libcint.git

It's very convenient to tryout Libcint with PySCF, which is a python
module for quantum chemistry program::

    http://github.com/sunqm/pyscf.git


Generating integrals
--------------------

If clisp was installed in the system, new integrals can be automatically
implemented.  You can add entries in ``script/auto_intor.cl`` and generate
code by::

    cd script
    clisp auto_intor.cl
    mv *.c ../src/autocode/

New entries should follow the format of those existed entries.
In one entry, you need to define the function name and the expression of
the integral.  The expression is consistent with Mulliken notation.
For one-electron integral, an entry can be::

    '("integral_name" spinor (number op-bra op-bra ... \| op-ket ...))

or::

    '("integral_name" spinor (number op-bra op-bra ... \| 1e-operator \| op-ket ...))


the entry of two-electron integral can be::

    '("integral_name" spinor (number op-bra-electron-1 ... \, op-ket-electron-1 ... \|
                                     op-bra-electron-2 ... \, op-ket-electron-2 ... ))

or::

    '("integral_name" spinor (number op-bra-electron-1 ... \, op-ket-electron-1 ... \|
                              r12 \| op-bra-electron-2 ... \, op-ket-electron-2 ... ))

* Parentheses must be paired.

* Line break is allowed.

* Note the _backslash_ in \| and \ is required.

* "integral_name" is the function name.  Valid name can be made up of
  letters, digits and underscore ("_").

* number can be an integer, a real number or a pure imaginary number. An
  imaginary number should be written as::

    #C(0 XXX)

* Supported operator-bra and operator-ket include

  - p     means    :math:`-i \nabla`
  - ip    means    :math:`\nabla`
  - r0    means    :math:`\vec{r} - (0,0,0)`
  - rc    means    :math:`\vec{r} - \vec{R}_(env[PTR_COMMON_ORIG])`
  - ri    means    :math:`\vec{r} - \vec{R}_i`
  - rj    means    :math:`\vec{r} - \vec{R}_j`
  - rk    means    :math:`\vec{r} - \vec{R}_k`
  - rl    means    :math:`\vec{r} - \vec{R}_l`
  - r              can be ri/rj/rk/rl; associate with the basis it operates
  - g     means    :math:`i/2 (\vec{R}_{bra} - \vec{R}_{ket}) \times \vec{r}`
  - sigma means    three pauli matrix
  - dot, cross     can be used to combine operator-bra or operator-ket

* Supported 1e-operator and 2e-operator include

  - rinv        means   :math:`1 / |\vec{r} - \vec{R}_(env[PTR_RINV_ORIG])|`
  - nuc         means   :math:`\sum_N Z_N / |\vec{r} - \vec{R}_N|`
  - nabla-rinv  means   :math:`\nabla (1 / |\vec{r} - \vec{R}_(env[PTR_RINV_ORIG])|)`
  - gaunt       means   :math:`\alpha_i \dot \alpha_j / |\vec{r}_i - \vec{r}_j|`
  - breit       means   :math:`-1/2\alpha_i \dot \alpha_j / |\vec{r}_i - \vec{r}_j| - 1/2 \alpha_i \dot r_{ij} \alpha_j \dot r_{ij} / |\vec{r}_i - \vec{r}_j|^3`

  Note sign - is not included in the gaunt integrals

Installation
------------

* Prerequisites

  - BLAS library
  - Python version 2.5 or higher (optional, for ``make test``)
  - Numpy (optional, for ``make test``)
  - clisp / SBCL (optional, for common lisp script)

* Build libcint::

    mkdir build; cd build
    cmake [-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>] ..
    make install

* Build libcint with examples and full or abridged tests (optional)::

    mkdir build; cd build
    cmake -DENABLE_EXAMPLE=1 -DENABLE_TEST=1 [-DQUICK_TEST=1] ..
    make
    make test ARGS=-V

* Build static library (optional)::

    mkdir build; cd build
    cmake -DBUILD_SHARED_LIBS=0 ..
    make install

* Compile with integer-8::

    mkdir build; cd build
    cmake -DI8=1 ..
    make install

* Long range part of range-separated Coulomb operator (optional)::

    mkdir build; cd build
    cmake -DWITH_RANGE_COULOMB ..
    make install


Available Integrals
-------------------

The available integrals can be found in the header file ``cint_funcs.h``. A simple
expression for each integral is also listed in the header file. The integral
function names and integral expressions correspond to the lisp symbol notations
in ``scripts/auto_intor.cl``

All integral functions have the same function signature: ::

    function_name(double *out, int *dims, int *shls, int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache);

Known problems
--------------

* Integral errors

  - Relative errors for regular ERIs are around 1e-12 and less.

  - Errors for short-range part of attenuated Coulomb interactions are generally
    larger than regular ERIs. Depending on the range-separation parameter,
    relative errors can reach 1e-10. However, comparing to computing integrals
    via "regular ERI - long-range ERI", errors are roughly one order of
    magnitude better.

  - Small integrals (< 1e-18 by default) are set to 0. If they are used in
    Schwarz inequality to estimate upper limit of an integral, the default
    integral cutoff might not be accurate enough. It can be adjusted by the
    parameter ``env[PTR_EXPCUTOFF]`` (since libcint 4.0). This parameter needs to be
    set to ``abs(log(cutoff_threshold))``.

* On 64-bit systems, ``make test`` stop with error: ::

    MKL FATAL ERROR: Cannot load libmkl_avx.so or libmkl_def.so.

  This problem is caused by the conflict between Python and MKL library.
  It can be fixed by adding ``-lmkl_avx`` or ``-lmkl_mc -lmkl_def`` to MKL link
  flags to replace the default BLAS link flags.  Be careful with the
  **order** of ``-lmkl_mc`` and ``-lmkl_def``.

* For basic ERIs, the code can handle highest angular momentum up to 7
  (present Rys-roots functions might be numerically unstable for
  nroots > 10 or l > 5).  But it has to be reduced to 5 or less for
  derivative or high order ERIs.  For every 4 derivative order,
  reduce 1 highest angular momentum for each shell.

* SIMD instructions can increase performance 5 ~ 50%.
  Please refer to **qcint** library (under GPL v3 license)::

        https://github.com/sunqm/qcint.git

* Tests and examples are not compiled by default. Compiling them by::

        cmake -DENABLE_EXAMPLE=1


How to cite
-----------

::

    @article{10.1002/jcc.23981,
      title = {Libcint: An efficient general integral library for Gaussian basis functions},
      author = {Sun, Qiming},
      journal = {Journal of Computational Chemistry},
      year = {2015},
      pages = {1664-1671},
      volume = {36},
      doi = {10.1002/jcc.23981},
      url = {http://dx.doi.org/10.1002/jcc.23981}
    }


Bug report
----------
Qiming Sun <osirpt.sun@gmail.com>

