### configure.ac
### Process this file with autoconf to produce a configure script
###
### Encore
### Nick Davis (c) 2011-2012
### nick-davis@utulsa.edu
### http://insilico.utulsa.edu/encore

AC_PREREQ(2.61)
AC_INIT(encore, 1.0, nick-davis@utulsa.edu, http://insilico.utulsa.edu/encore)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE(encore, 1.0)
AC_CONFIG_SRCDIR([encore.cpp])
AC_CONFIG_HEADER([config.h])
AC_LANG_CPLUSPLUS

if test -z $CXXFLAGS; then
	CXXFLAGS='-O3'
fi

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL

# Use libtool for all builds
LT_INIT

AC_CANONICAL_BUILD
AC_CANONICAL_HOST
# Checks OS
case "$host_os" in
  darwin*)
        # Mac OS X
        OS_CXX="-DUNIX"
        SYS_MAC=1
        ;;
  *solaris*)
  		# Solaris
  		OS_CXX="-fast -xtarget=ultraT2 -xdepend=no"
  		SYS_SOLARIS=1
		# FIXME: Replace `main' with a function in `-lgcc':
		AC_CHECK_LIB([gcc], [main])
		# FIXME: Replace `main' with a function in `-lsocket':
		AC_CHECK_LIB([socket], [main])
		# FIXME: Replace `main' with a function in `-lnsl':
		AC_CHECK_LIB([nsl], [main])
        ;;

  *mingw32* | *cygwin*)
  		# Windows
        OS_CXX="-DWIN"
        SYS_WIN=1
        ;;
    *)
        #Default Case, assume UNIX
        OS_CXX="-DUNIX"
        SYS_UNIX=1
        ;;
esac

# OS conditionals
AM_CONDITIONAL([SYS_MAC], [test "x$SYS_MAC" = x1])
AM_CONDITIONAL([SYS_SOLARIS], [test "x$SYS_SOLARIS" = x1])
AM_CONDITIONAL([SYS_WIN], [test "x$SYS_WIN" = x1])
AM_CONDITIONAL([SYS_UNIX], [test "x$SYS_UNIX" = x1])

# Set OS C++ flags in Makefile.am
AC_SUBST(OS_CXX)

# Checks for libraries.
AC_CHECK_LIB([z], [gzopen],
	[],
	AC_MSG_ERROR(test for libz failed)
)

AC_CHECK_LIB([armadillo], [main],
	[],
	AC_MSG_ERROR(test for armadillo failed)
)


AC_CHECK_LIB([blas], [dgemm_],
	[],
	AC_MSG_ERROR(test for blas failed)
)


AC_CHECK_LIB([lapack], [dgesdd_],
	[],
	AC_MSG_ERROR(test for lapack failed)
)


AC_CHECK_LIB([plink], [main],
	[],
	AC_MSG_ERROR(test for plink failed)
)


AC_ARG_WITH([openmp],
	[AS_HELP_STRING([--with-openmp],
	[enable support for OpenMP parallelism])],
	[],
	[with_openmp=yes])
AS_IF([test "x$with_openmp" = xyes],
	[AC_CHECK_LIB([gomp], [GOMP_parallel_loop_dynamic_start],
		[],
		[if test "x$with_openmp" != xcheck; then
			AC_MSG_ERROR(
				[--with-openmp was given, but test for openmp failed])
		fi
		], [])])

# boost.m4
BOOST_REQUIRE([1.34])
BOOST_FIND_HEADER(boost/program_options.hpp, AC_MSG_ERROR([could not find required headers of Boost c++. Install boost c++ libraries.],))
BOOST_PROGRAM_OPTIONS
BOOST_SYSTEM
BOOST_FILESYSTEM

# Random Jungle
AC_CHECK_LIB([rjungle], [main])
# libr, used in Random Jungle
AC_CHECK_LIB([lr], [main])

# GNU Scientific Library (GSL)
# gsl.m4
AX_PATH_GSL(1.8,, AC_MSG_ERROR([could not find required version of GSL]))

# EC
AC_CHECK_LIB([ec], [main])
			 
AC_ARG_ENABLE([force32bit],
	[AS_HELP_STRING([--enable-force32bit],
	[force 32-bit code for 64-bit architectures])],
	[case "${enableval}" in
	  yes) force32bit=true ;;
	  no)  force32bit=false ;;
	  *) AC_MSG_ERROR([bad value ${enableval} for --enable-force32bit]) ;;
	esac],
	[force32bit=false])
AC_ARG_ENABLE([forcedynamic],
	[AS_HELP_STRING([--enable-forcedynamic],
	[force dynamic libraries])],
	[case "${enableval}" in
	  yes) forcedynamic=true ;;
	  no)  forcedynamic=false ;;
	  *) AC_MSG_ERROR([bad value ${enableval} for --enable-forcedynamic]) ;;
	esac],
	[forcedynamic=true])
# library and package conditionals
AM_CONDITIONAL([WITH_OPENMP], [test "x$with_openmp" = xyes])
AM_CONDITIONAL([FORCE_32BIT], [test "x$force32bit" = xtrue])
AM_CONDITIONAL([FORCE_DYNAMIC], [test "x$forcedynamic" = xtrue])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([arpa/inet.h fcntl.h netdb.h netinet/in.h stddef.h stdlib.h string.h sys/socket.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_C_VOLATILE

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_HEADER_MAJOR
AC_FUNC_STAT
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor memset pow socket sqrt strdup])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
