lib_LIBRARIES = libRNAlila.a

pkginclude_HEADERS = 	lila.h \
			moves.h\
			topology.h

libRNAlila_a_SOURCES = 	io.c \
			lila.c \
			ds_utils.c \
			moves.c\
			topology.c\
			ranstrucS.c


AM_CFLAGS =  ${Glib2_CFLAGS} ${GSL_CFLAGS} ${ViennaRNA_CFLAGS} -g3 -O0
AM_CPPFLAGS = ${Glib2_LIBS} -I${includedir} -I.

#libRNAlila_a_CFLAGS = $(AM_CFLAGS) 

EXTRA_DIST=VERSION.txt

# End of file
