bin_PROGRAMS = RNAwalk

LDADD = $(top_srcdir)/src/RNAlila/libRNAlila.a ${ViennaRNA_LIBS} ${Glib2_LIBS}

AM_CFLAGS =  ${Glib2_CFLAGS} ${GSL_CFLAGS} ${ViennaRNA_CFLAGS}
AM_CPPFLAGS = -I$(top_srcdir)/src  ${Glib2_LIBS} -I${includedir} -I.

GENGETOPT_CMDL =  RNAwalk_cmdl.c RNAfold_cmdl.h

GENGETOPT_FILES =  RNAwalk.ggo

EXTRA_DIST =  ${GENGETOPT_FILES} ${GENGETOPT_CMDL}

# remove all generated man pages
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)


RNAwalk.c:    RNAwalk_cmdl.h

RNAwalk_cmdl.h  RNAwalk_cmdl.c: RNAwalk.ggo
	$(AM_V_GEN)gengetopt -i RNAwalk.ggo

RNAwalk_SOURCES=RNAwalk_cmdl.c RNAwalk.c

