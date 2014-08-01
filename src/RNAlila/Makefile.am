lib_LIBRARIES = libRNAlila.a

pkginclude_HEADERS = 	lila.h \
			moves.h

libRNAlila_a_SOURCES = 	lila.c \
			moves.c

AM_CFLAGS =  ${Glib2_CFLAGS} ${GSL_CFLAGS} ${ViennaRNA_CFLAGS} -g3 -O0
AM_CPPFLAGS = ${Glib2_LIBS} -I${includedir} -I.

EXTRA_DIST=VERSION.txt

# End of file
