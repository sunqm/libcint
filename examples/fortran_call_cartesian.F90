integer function factorial(n)
integer :: n
integer :: i
factorial = 1
do i = 1, n
  factorial = factorial * i
end do
end function factorial

! general contracted DZ basis [3s1p/2s1p] for H2
!     exponents    contract-coeff
! S   6.0          0.7               0.4
!     2.0          0.6               0.3
!     0.8          0.5               0.2
! P   0.9          1.

program cartesian
implicit none
integer,parameter :: CHARGE_OF  = 1
integer,parameter :: PTR_COORD  = 2
integer,parameter :: NUC_MOD_OF = 3
integer,parameter :: PTR_ZETA   = 4
integer,parameter :: ATM_SLOTS  = 6

integer,parameter :: ATOM_OF    = 1
integer,parameter :: ANG_OF     = 2
integer,parameter :: NPRIM_OF   = 3
integer,parameter :: NCTR_OF    = 4
integer,parameter :: KAPPA_OF   = 5
integer,parameter :: PTR_EXP    = 6
integer,parameter :: PTR_COEFF  = 7
integer,parameter :: BAS_SLOTS  = 8

integer,parameter :: PTR_ENV_START = 20

integer :: natm = 2
integer :: nbas = 4
integer,allocatable :: atm(:,:)
integer,allocatable :: bas(:,:)
double precision,allocatable :: env(:)
double precision,external :: CINTgto_norm

integer :: n, off
integer :: i, j, k, l
integer :: di, dj, dk, dl
integer :: shls(4)
double precision,allocatable :: buf1e(:,:,:), buf2e(:,:,:,:,:)
integer,external :: CINTcgto_cart
integer,external :: cint1e_ipovlp_cart, cint2e_ip1_cart
! On 32-bit machine, integer(4) :: opt
integer(8) :: opt
allocate (atm(ATM_SLOTS,natm))
allocate (bas(BAS_SLOTS,nbas))
allocate (env(10000))

off = PTR_ENV_START

i = 1
atm(CHARGE_OF,i) = 1
atm(PTR_COORD,i) = off ! note the 0-based index
env(off + 1) =  0 ! x (Bohr)
env(off + 2) =  0 ! y (Bohr)
env(off + 3) =-.8 ! z (Bohr)
i = i + 1
off = off + 3
atm(CHARGE_OF,i) = 1
atm(PTR_COORD,i) = off
env(off + 1) = 0
env(off + 2) = 0
env(off + 3) =.8 ! (Bohr)
i = i + 1
off = off + 3

n = 1
! basis #1, 
bas(ATOM_OF  ,n)  = 0 ! note that it's the first atom, the index is 0-based
bas(ANG_OF   ,n)  = 0
bas(NPRIM_OF ,n)  = 3
bas(NCTR_OF  ,n)  = 2
bas(PTR_EXP  ,n)  = off ! note the 0-based index
env(off + 1) = 6.
env(off + 2) = 2.
env(off + 3) = .8
off = off + 3
bas(PTR_COEFF,n) = off
env(off + 1) = .7 * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+1))
env(off + 2) = .6 * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+2))
env(off + 3) = .5 * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+3))
env(off + 4) = .4 * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+1))
env(off + 5) = .3 * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+2))
env(off + 6) = .2 * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+3))
off = off + 6
n = n + 1

! basis #2
bas(ATOM_OF  ,n)  = 0
bas(ANG_OF   ,n)  = 1
bas(NPRIM_OF ,n)  = 1
bas(NCTR_OF  ,n)  = 1
bas(PTR_EXP  ,n)  = off
env(off + 1) = .9
off = off + 1
bas(PTR_COEFF,n) = off
env(off + 1) = 1. * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+1))
off = off + 1
n = n + 1

! basis #3 == basis #1
bas(ATOM_OF  ,n) = 1 ! note that it's the second atom, the index is 0-based
bas(ANG_OF   ,n) = bas(ANG_OF   ,1)
bas(NPRIM_OF ,n) = bas(NPRIM_OF ,1)
bas(NCTR_OF  ,n) = bas(NCTR_OF  ,1)
bas(PTR_EXP  ,n) = bas(PTR_EXP  ,1)
bas(PTR_COEFF,n) = bas(PTR_COEFF,1)
n = n + 1

! basis #4 == basis #2
bas(ATOM_OF  ,n) = 1
bas(ANG_OF   ,n) = bas(ANG_OF   ,2)
bas(NPRIM_OF ,n) = bas(NPRIM_OF ,2)
bas(NCTR_OF  ,n) = bas(NCTR_OF  ,2)
bas(PTR_EXP  ,n) = bas(PTR_EXP  ,2)
bas(PTR_COEFF,n) = bas(PTR_COEFF,2)
n = n + 1


!
! call one-electron spinor integrals
! note the index of shell is 0-based
! the integral has 3 components
!
i = 0; shls(1) = i; di = CINTcgto_cart(i, bas)
j = 1; shls(2) = j; dj = CINTcgto_cart(j, bas)
allocate (buf1e(di,dj,3))
if (0 /= cint1e_ipovlp_cart(buf1e, shls, atm, natm, bas, nbas, env)) then
  print*, "This gradient integral is not 0.\n"
else
  print*, "This integral is 0.\n"
endif
deallocate (buf1e)

!
! call two-electron spinor integrals
! the index of shell is 0-based
! the integral has 3 components
!
i = 0; shls(1) = i; di = CINTcgto_cart(i, bas)
j = 1; shls(2) = j; dj = CINTcgto_cart(j, bas)
k = 2; shls(3) = k; dk = CINTcgto_cart(k, bas)
l = 2; shls(4) = l; dl = CINTcgto_cart(l, bas)
allocate (buf2e(di,dj,dk,dl,3))
if (0 /= cint2e_ip1_cart(buf2e, shls, atm, natm, bas, nbas, env, 0_8)) then
  print*, "This gradient integral is not 0.\n"
else
  print*, "This integral is 0.\n"
endif
deallocate (buf2e)

call cint2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env)
i = 0; shls(1) = i; di = CINTcgto_cart(i, bas)
j = 1; shls(2) = j; dj = CINTcgto_cart(j, bas)
k = 2; shls(3) = k; dk = CINTcgto_cart(k, bas)
l = 2; shls(4) = l; dl = CINTcgto_cart(l, bas)
allocate (buf2e(di,dj,dk,dl,3))
if (0 /= cint2e_ip1_cart(buf2e, shls, atm, natm, bas, nbas, env, opt)) then
  print*, "This gradient integral is not 0.\n"
else
  print*, "This integral is 0.\n"
endif
deallocate (buf2e)
call CINTdel_optimizer(opt)
deallocate (atm, bas, env)
end program cartesian
