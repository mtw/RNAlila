#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([RNAlila],
	m4_esyscmd_s([test -d .git && git describe --tags || cat VERSION]),
	[michael@wolfinger.eu])

# for automake
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE

#AC_CONFIG_SRCDIR([src/RNAlila/moves.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CANONICAL_HOST

# check for custom libs
PKG_CHECK_MODULES(ViennaRNA, RNAlib2 >= 2.2)

# check system programs
AC_PATH_PROG([MAKE],[make],[make])
AC_PATH_PROG([AR],[ar],[ar])

# Checks for libraries.
AC_LANG(C)
AC_CHECK_LIB([m], [isnan])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl], [gsl_histogram_alloc])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h unistd.h assert.h errno.h math.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_CONST
AC_C_INLINE
AC_EXEEXT
AC_OBJEXT
AC_TYPE_SIZE_T
AC_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memset pow strdup strcpy floor])

AC_CONFIG_FILES([Makefile src/Makefile src/RNAlila/Makefile])
AC_CONFIG_FILES([RNAlila.pc])
AC_OUTPUT

# Inform user
AC_MSG_NOTICE(
[
You can run make now.
])
