## Initialization 
AC_INIT([poloka-sub], 0.1.0, [poloka team])
AC_CONFIG_SRCDIR([poloka/kernelfit.cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR(autoconf)
AM_INIT_AUTOMAKE([foreign])

## Checks for programs to build poloka-core library
AC_PROG_CXX
LT_INIT([disable-static])

## Check for mandatory poloka-core
PKG_CHECK_MODULES([POLOKA_CORE],
		  [poloka-core],,
		  [AC_MSG_ERROR([Cound not find the poloka-core library])])

AC_CONFIG_FILES([Makefile
		 poloka-sub.pc
                 datacards/Makefile
                 doc/Doxyfile
                 doc/Makefile
                 poloka/Makefile
		 tools/Makefile
                 ])
AC_OUTPUT
