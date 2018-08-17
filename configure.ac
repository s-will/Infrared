AC_PREREQ([2.65])
AC_INIT([InfraRed], [0.3], [], [infrared], [http://www.tbi.univie.ac.at/~will] )

dnl minimal required version of the boost library 
BOOST_REQUIRED_VERSION=1.50.0

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADERS([src/config.h:config.h.in])

AC_CONFIG_MACRO_DIR([m4])

dnl include libtool
LT_PREREQ([2.2])
LT_INIT([shared])

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_CC

dnl automake initialisation (mandatory) and check for minimal automake API version
AM_INIT_AUTOMAKE([1.11])

dnl use the C++ compiler for the following checks
AC_LANG([C++])

dnl ----------------------------------------
dnl generate doxygen documentation
dnl

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN($PACKAGE_NAME, ["Doc/doxygen.cfg"], ["Doc"])
dnl ----------------------------------------


dnl --------------------
dnl Debug option
dnl
AC_MSG_CHECKING([whether to build with debug information])
debug_it=no
AC_ARG_ENABLE([debug],
    AC_HELP_STRING(
        [--enable-debug],
        [enable debug data generation (def=no)]
    ),
    debug_it="$enableval"
)
AC_MSG_RESULT([$debug_it])
if test "$debug_it" = "yes"; then
    AC_DEFINE([DEBUG],1,[Turn on debugging])
    AM_CXXFLAGS="$AM_CXXFLAGS -g"
else
    AC_DEFINE([NDEBUG],1,[Turn off debugging])
fi

AM_CXXFLAGS="$AM_CXXFLAGS -Wall -pedantic"


dnl ensure we are using c14 C++ standard
AX_CXX_COMPILE_STDCXX([14], [noext], [mandatory])

dnl ########################################
dnl HTD tree decomposition library
dnl ########################################
AC_ARG_WITH([htd],
    AC_HELP_STRING(
        [--with-htd=PREFIX],
        [Installation directory of the HTD tree decomposition library]
      )
)

if test "$HTD_PREFIX" = "" ; then
    AC_MSG_CHECKING([where to find the htd tree decomposition library])
    if test "$with_htd" = "" ; then
      if test -d "/usr/local/include/htd" ; then
        HTD_PREFIX="/usr/local"
      else
        HTD_PREFIX="/usr"
      fi
    else
      HTD_PREFIX="$with_htd"
    fi
    AC_MSG_RESULT([$HTD_PREFIX])
fi

HTD_CPPFLAGS="-I$HTD_PREFIX/include"
HTD_LDFLAGS="-L$HTD_PREFIX/lib"
HTD_LIBS="-lhtd"

AC_SUBST([HTD_CPPFLAGS])
AC_SUBST([HTD_LDFLAGS])
AC_SUBST([HTD_LIBS])

dnl ##############################################################################
dnl BOOST CHECK
dnl ##############################################################################

PYTHON_VERSION=3.6

AX_BOOST_BASE([$BOOST_REQUIRED_VERSION], [FOUND_BOOST=1;], [FOUND_BOOST=0;])
AX_BOOST_PYTHON
dnl AX_BOOST_GRAPH

dnl ###########  CHECKS  ############################################

dnl Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

dnl Checks for header files.
AC_HEADER_STDC


dnl #########################################################################
dnl check boost test results
dnl #########################################################################

dnl FOUND_BOOST is only defined if want_boost is "yes"
AS_IF([test $want_boost = "no" || test $FOUND_BOOST != 1], [
        AC_MSG_NOTICE([])
        AC_MSG_NOTICE([The Boost Library was not found!])
        AC_MSG_NOTICE([ -> If installed in a non-standard path, please use '--with-boost=PREFIX'.])
        AC_MSG_NOTICE([])
        DEPENDENCYNOTFOUND=1;
], [
        AM_CXXFLAGS="$BOOST_CPPFLAGS $AM_CXXFLAGS"
        AM_LDFLAGS="$BOOST_LDFLAGS $AM_LDFLAGS"
        LIBS="$LIBS"
])

AM_CXXFLAGS="$PYTHON_CPPFLAGS $AM_CXXFLAGS"
AM_LDFLAGS="$PYTHON_LDFLAGS $AM_LDFLAGS"
LIBS="$LIBS -l$BOOST_PYTHON_LIB"

dnl distribute additional compiler and linker flags
dnl --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])
AM_CONDITIONAL(enable_pkg_config, [test "x$HAVE_PKG_CONFIG" = "xyes"])


dnl files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/Infrared/Makefile])
AC_CONFIG_FILES([src/Redprint/Makefile])
AC_CONFIG_FILES([src/Test/Makefile])

dnl generate the final Makefile etc.
AC_OUTPUT
