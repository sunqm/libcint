
AM_CPPFLAGS= -I$(top_srcdir)/src/include
AM_LDFLAGS = -L../src -lcint

noinst_PROGRAMS = cartesian_c spheric_c spinor_c \
		  cartesian_f spheric_f spinor_f \
		  time_c2h6 time_c60

time_c2h6_SOURCES = time_c2h6.c
time_c60_SOURCES = time_c2h6.c
cartesian_c_SOURCES = c_call_cartesian.c
spheric_c_SOURCES =   c_call_spheric.c
spinor_c_SOURCES =    c_call_spinor.c
cartesian_f_SOURCES = fortran_call_cartesian.F90
spheric_f_SOURCES =   fortran_call_spheric.F90
spinor_f_SOURCES =    fortran_call_spinor.F90

