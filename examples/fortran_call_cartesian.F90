integer function factorial(n)
integer :: n
integer :: i
factorial = 1
do i = 1, n
  factorial = factorial * i
end do
end function factorial

double precision function gto_norm(n, a)
! normalization factor of function r^n e^{-a r^2}
integer :: n
double precision :: a
integer,external :: factorial
gto_norm = 2**(2*n+3) * factorial(n+1) * (2*a)**(n+1.5) &
        / (factorial(2*n+2) * sqrt(3.14159265358979d0))
end function gto_norm

program cartesian
implicit none
integer,parameter :: CHARGE_OF  = 1
integer,parameter :: PTR_COORD  = 2
integer,parameter :: NUC_MOD_OF = 3
integer,parameter :: PTR_MASS   = 4
integer,parameter :: RADI_GRIDS = 5
integer,parameter :: ANG_GRIDS  = 6
integer,parameter :: ATM_SLOTS  = 6

integer,parameter :: ATOM_OF    = 1
integer,parameter :: ANG_OF     = 2
integer,parameter :: NPRIM_OF   = 3
integer,parameter :: NCTR_OF    = 4
integer,parameter :: KAPPA_OF   = 5
integer,parameter :: PTR_EXP    = 6
integer,parameter :: PTR_COEFF  = 7
integer,parameter :: GAUGE_OF   = 8
integer,parameter :: BAS_SLOTS  = 8

integer,parameter :: PTR_ENV_START = 20

integer :: natm = 2
integer :: nbas = 3
integer,allocatable :: atm(:,:)
integer,allocatable :: bas(:,:)
double precision,allocatable :: env(:)
double precision,external :: gto_norm

integer :: n, off
integer :: i, j, k, l
integer :: di, dj, dk, dl
integer :: shls(4)
double precision,allocatable :: buf1e(:,:,:), buf2e(:,:,:,:,:)
integer,external :: cgtos_cart
integer,external :: cint1e_ipovlp_cart, cint2e_ip1_cart
allocate (atm(ATM_SLOTS,natm))
allocate (bas(BAS_SLOTS,nbas))
allocate (env(10000))

off = PTR_ENV_START
do i = 1, natm
  atm(CHARGE_OF,i) = i
  atm(PTR_COORD,i) = off ! off is the offset in env
  env(off + 1) = i - 1
  env(off + 2) = i - 1
  env(off + 3) = i - 1
  off = off + 3
end do

n = 1

! basis #0 with kappa > 0  => p_1/2
bas(ATOM_OF  ,n)  = 0 ! the first atom, the index is 0-based
bas(ANG_OF   ,n)  = 1
bas(KAPPA_OF ,n)  = 1
bas(NPRIM_OF ,n)  = 1
bas(NCTR_OF  ,n)  = 1
bas(PTR_EXP  ,n)  = off ! offset of exponents in env
env(off + 1) = 1.d0
bas(PTR_COEFF,n) = off + 1 ! offset of contraction coefficeints
env(off + 2) = 1.d0 * gto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)))
off = off + 2
n = n + 1

! basis #1 with kappa = 0  => d_3/2, d_5/2,
! 2 primitive-GTO -> 2 contracted-GTO
bas(ATOM_OF  ,n)  = 0 ! the first atom
bas(ANG_OF   ,n)  = 2
bas(KAPPA_OF ,n)  = 0
bas(NPRIM_OF ,n)  = 2
bas(NCTR_OF  ,n)  = 2
bas(PTR_EXP  ,n)  = off
env(off + 0) = 3.d0
env(off + 1) = 5.d0
bas(PTR_COEFF,n) = off + 2
env(off + 2) = 1.d0 * gto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)))
env(off + 3) = 2.d0 * gto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+1))
env(off + 4) = 4.d0 * gto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)))
env(off + 5) = 8.d0 * gto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+1))
off = off + 6
n = n + 1

! basis #2 with kappa < 0  => f_5/2
bas(ATOM_OF  ,n)  = 1 ! the second atom
bas(ANG_OF   ,n)  = 3
bas(KAPPA_OF ,n)  = -4
bas(NPRIM_OF ,n)  = 1
bas(NCTR_OF  ,n)  = 1
bas(PTR_EXP  ,n)  = off
env(off + 0) = 1.d0
bas(PTR_COEFF,n) = off + 1
env(off + 1) = 1.d0 * gto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)))
off = off + 2
n = n + 1

!
! call one-electron spinor integrals
! the index of shell is 0-based
! the integral has 3 components
!
i = 0; shls(1) = i; di = cgtos_cart(i, bas)
j = 1; shls(2) = j; dj = cgtos_cart(j, bas)
allocate (buf1e(di,dj,3))
if (0 /= cint1e_ipovlp_cart(buf1e, shls, atm, natm, bas, nbas, env)) then
  print*, "This gradient integral is not 0.\n"
endif
deallocate (buf1e)

!
! call two-electron spinor integrals
! the index of shell is 0-based
! the integral has 3 components
!
i = 0; shls(1) = i; di = cgtos_cart(i, bas)
j = 1; shls(2) = j; dj = cgtos_cart(j, bas)
k = 2; shls(3) = k; dk = cgtos_cart(k, bas)
l = 2; shls(4) = l; dl = cgtos_cart(l, bas)
allocate (buf2e(di,dj,dk,dl,3))
if (0 /= cint2e_ip1_cart(buf2e, shls, atm, natm, bas, nbas, env)) then
  print*, "This gradient integral is not 0.\n"
endif
deallocate (buf2e)
deallocate (atm, bas, env)
end program cartesian
