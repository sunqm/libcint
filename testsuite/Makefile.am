
libtestcint_la_CFLAGS = -g -O2

noinst_LTLIBRARIES = libtestcint.la
libtestcint_la_SOURCES = init_test_env.c
libtestcint_la_LIBADD = -L../src -lcint
libtestcint_la_LDFLAGS = -module -shared -rpath/nowhere -avoid-version -rpath $(abs_builddir)
libtestcint_la_CPPFLAGS = -I$(srcdir)/../src

TESTS_ENVIRONMENT = buildir=$(top_builddir) python
check_SCRIPTS = test_cint.py
TESTS = $(check_SCRIPTS)
