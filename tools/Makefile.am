CLEANFILES = *~
MAINTAINERCLEANFILES = \
        Makefile.in \
        stamp-*

AM_CPPFLAGS = -I$(top_srcdir) @POLOKA_CORE_CFLAGS@

AM_DEFAULT_SOURCE_EXT = .cc

bin_PROGRAMS = pka-imdiff pka-makesub

LDADD = $(top_builddir)/poloka/libpoloka-sub.la
