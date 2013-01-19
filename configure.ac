#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([libcint], [1.0], [osirpt.sun@gmail.com])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_CONFIG_SRCDIR([src/cint_bas.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC
AC_PROG_FC

#AX_PYTHON
AX_BLAS

# new version autoconf provides
#LT_INIT([disable-static])
AM_PROG_LIBTOOL([disable-static])

AC_CHECK_HEADERS([stdlib.h math.h])

AC_FUNC_MALLOC
AC_CHECK_LIB([m], [sqrt])

AC_ARG_ENABLE([auto-intor],
              [AS_HELP_STRING([--enable-auto-intor],
                              [enable code-generate script [default=no]])],
              [case $enableval in
                 yes) enable_auto_intor=yes ;;
                  no) enable_auto_intor=no ;;
                   *) AC_MSG_ERROR([bad value $enableval for --enable-auto-intor]) ;;
               esac],
              [enable_auto_intor=no])
AM_CONDITIONAL([WITH_INTOR], [test "x$enable_auto_intor" = xyes])
if test "x$enable_auto_intor" = xyes ; then
    AC_CHECK_PROG(clisp_found, clisp, true,
                  AC_MSG_ERROR([clisp is not found.]))
fi

## output files
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 testsuite/Makefile
                 examples/Makefile])
AC_OUTPUT

