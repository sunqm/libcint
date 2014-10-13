
AM_CFLAGS = -DNDEBUG -funroll-all-loops
AM_CPPFLAGS= -I$(top_srcdir)/src

lib_LTLIBRARIES = libcint.la

libcint_la_LDFLAGS = -release 2.0.6
libcint_la_LIBADD = $(BLAS_LIBS) $(LIBS)
libcint_la_SOURCES = c2f.c \
                     cart2sph.c \
                     cint1e.c  \
                     cint2e.c \
                     cint_bas.c \
                     fblas.c \
                     g1e.c \
                     g2e.c \
                     rys_roots.c \
                     misc.c \
                     optimizer.c \
                     autocode/auto_intor1.c \
                     autocode/auto_intor2.c \
                     autocode/auto_intor3.c \
                     autocode/grad1.c \
                     autocode/grad2.c \
                     autocode/grad3.c \
                     autocode/gaunt1.c

