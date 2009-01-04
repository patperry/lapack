
AC_PREREQ(2.61)
AC_INIT([Haskell lapack], [0.1], [patperry@stanford.edu])
AC_CONFIG_SRCDIR([lapack.cabal])

AC_CONFIG_HEADERS([cbits/config.h])

AC_PROG_CC
AC_CONFIG_MACRO_DIR([m4])

AC_ARG_WITH(compiler, 
    [AC_HELP_STRING([--with-compiler=<cmd>], [use Haskell compiler <cmd>])])

AX_LAPACK
if test "$ax_lapack_ok" = yes; then
  BUILD_PACKAGE_BOOL=True
else 
  BUILD_PACKAGE_BOOL=False
fi
AC_SUBST([BUILD_PACKAGE_BOOL])

AC_CONFIG_FILES([lapack.buildinfo])
AC_OUTPUT
