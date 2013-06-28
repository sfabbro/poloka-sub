CLEANFILES = *~
MAINTAINERCLEANFILES = \
	Makefile.in \
	stamp-*

lib_LTLIBRARIES = libpoloka-sub.la

src_includedir = $(includedir)/poloka
src_include_HEADERS = \
	convolution.h \
	convtemplate.h \
	detection.h \
	imagepair.h \
	imagesubtraction.h \
	kernelfit.h \
	kernelfitter.h \
	quali_box.h \
	scorescollection.h \
	sub.h \
	subtraction.h

libpoloka_sub_la_SOURCES = \
	$(src_include_HEADERS) \
	convolution.cc \
	detection.cc \
	imagepair.cc \
	imagesubtraction.cc \
	kernelfit.cc \
	kernelfitter.cc \
	quali_box.cc \
	scorescollection.cc \
	sub.cc \
	subtraction.cc

libpoloka_sub_la_CPPFLAGS = -I$(top_srcdir) @POLOKA_CORE_CFLAGS@ 

libpoloka_sub_la_LIBADD = @POLOKA_CORE_LIBS@
