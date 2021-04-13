
# general contracted DZ basis [3s1p/2s1p] for H₂
#     exponents    contract-coeff
# S   6.0          0.7               0.4
#     2.0          0.6               0.3
#     0.8          0.5               0.2
# P   0.9          1.

function gto_norm(n, a)
   # normalization factor of function rⁿ exp(-ar²)
    s = 2^(2n+3) * factorial(n+1) * (2a)^(n+1.5) / (factorial(2n+2) * √π)
    sqrt(s)
end

CHARGE_OF  = 1
PTR_COORD  = 2
NUC_MOD_OF = 3
PTR_ZETA   = 4
ATM_SLOTS  = 6

ATOM_OF    = 1
ANG_OF     = 2
NPRIM_OF   = 3
NCTR_OF    = 4
KAPPA_OF   = 5
PTR_EXP    = 6
PTR_COEFF  = 7
BAS_SLOTS  = 8

ptr_env = 20
env = fill(0.0, ptr_env)

i = 0
#                     CHARGE_OF,PTR_COORD
atm =           Int32[1         ptr_env  0 0 0 0]
ptr_env += 3
atm = vcat(atm, Int32[1         ptr_env  0 0 0 0])
ptr_env += 3
#             x  y   z (Bohr)
append!(env, [0, 0, -0.8])
append!(env, [0, 0,  0.8])

# basis for atom #0
# 3s -> 2s
append!(env, [6., 2., .8])
append!(env, [.7 * gto_norm(0, 6.),
              .6 * gto_norm(0, 2.),
              .5 * gto_norm(0, .8),
              .4 * gto_norm(0, 6.),
              .3 * gto_norm(0, 2.),
              .2 * gto_norm(0, .8)])
#                     ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, KAPPA_OF, PTR_EXP, PTR_COEFF
bas =           Int32[0        0       3         2         0        ptr_env  ptr_env+3  0]
ptr_env += 9
bas = vcat(bas, Int32[0        1       1         1        0         ptr_env  ptr_env+1  0])
append!(env, [.9])
append!(env, [1. * gto_norm(1, 0.9)])
ptr_env += 2

# basis functions for atom #1, they are the same to thoes of atom #0
bas = vcat(bas, bas)

shls = Int32[0, 1]

natm = length(atm)
nbas = length(bas)

@info "Input" atm bas #env

bas = permutedims(bas, (2,1))   # Column (Julia) -> Row (C) order
atm = permutedims(atm, (2,1))   # Column (Julia) -> Row (C) order

# Libcint API
const LIBCINT = joinpath(@__DIR__, "path/to/libcint.so/dylib/dll")

function CINTcgtos_spheric(bas_id, bas)
   @static if VERSION >= v"1.5.0" # clearer syntax
      @ccall LIBCINT.CINTcgtos_spheric(bas_id::Cint, bas::Ptr{Cint})::Cint
   else # works anyway
      ccall( (:CINTcgtos_spheric, LIBCINT), Cint, (Cint, Ptr{Cint},), bas_id, bas)
   end
end

function cint1e_ipnuc_sph(buf, shls, atm, natm, bas, nbas, env)
   @static if VERSION >= v"1.5.0" # clearer syntax
      @ccall LIBCINT.cint1e_ipnuc_sph(
                                      buf  :: Ptr{Cdouble},
                                      shls :: Ptr{Cint},
                                      atm  :: Ptr{Cint},
                                      natm :: Cint,
                                      bas  :: Ptr{Cint},
                                      nbas :: Cint,
                                      env  :: Ptr{Cdouble}
                                     )::Cvoid
   else # works anyway
      ccall( (:cint1e_ipnuc_sph, LIBCINT), Cvoid, (
                                      Ptr{Cdouble},
                                      Ptr{Cint},
                                      Ptr{Cint},
                                      Cint,
                                      Ptr{Cint},
                                      Cint,
                                      Ptr{Cdouble}, ),
                                      buf  ,
                                      shls ,
                                      atm  ,
                                      natm ,
                                      bas  ,
                                      nbas ,
                                      env  ,
                                     )
   end
end

di = CINTcgtos_spheric(0, bas) # C uses 0-indexing!
dj = CINTcgtos_spheric(1, bas) # C uses 0-indexing!
@assert di == (2bas[ANG_OF,1]+1) * bas[NCTR_OF,1]
@assert dj == (2bas[ANG_OF,2]+1) * bas[NCTR_OF,2]

buf = zeros(3, dj, di) # Column (Julia) -> Row (C) order

cint1e_ipnuc_sph(buf, shls, atm, natm, bas, nbas, env)
@info "∇⟨i|nuc|j⟩ = " permutedims(buf, (3,2,1))
