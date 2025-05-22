
# 
# Copyright (C) 2025  PML
#
# Created from https://github.com/fplll/fplll/blob/master/configure.ac
#
#   This file is part of pml. pml is free software: you
#   can redistribute it and/or modify it under the terms of the GNU Lesser
#   General Public License as published by the Free Software Foundation,
#   either version 2.1 of the License, or (at your option) any later version.
#
#   pml is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with pml. If not, see <http://www.gnu.org/licenses/>. 


# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([pml],[0.0.1])
AC_CONFIG_SRCDIR([flint-extras/nmod_extra/src/nmod_extra.c])
AM_MAINTAINER_MODE([enable])
AC_CANONICAL_TARGET


AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_HEADERS([config.h])

#GV do it cleaner for local include, compilation in lib (temporary) 
DEFAULT_CPPFLAGS="-O3 -Wall -I/usr/local/include -I../include" 

CPPFLAGS=${CPPFLAGS:-$DEFAULT_CPPFLAGS}

DEFAULT_LDFLAGS="-L/usr/local/lib"

LDFLAGS=${LDFLAGS:-$DEFAULT_LDFLAGS}

#Defaults 


LT_INIT


AC_PROG_INSTALL

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AX_CXX_COMPILE_STDCXX([11],[noext],[mandatory]) # GV outdated ? 


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_C_VOLATILE

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow rint sqrt strtol])


CHECK_ATOMIC()

#################################################
AS_BOX([          PML config               ])
#################################################


# Check gmp
#----------

AC_ARG_WITH(gmp,
      AS_HELP_STRING([--with-gmp=@<:@=DIR@:>@], [GMP install directory]), [
         CPPFLAGS="$CPPFLAGS -I$withval/include"
         LDFLAGS="$LDFLAGS -L$withval/lib"
      ])

AC_CHECK_LIB(gmp, __gmp_get_memory_functions, , [AC_MSG_ERROR(
   [GMP version >= 4.2.0 needed, see http://gmplib.org])])


# Check flint
#------------

AC_ARG_WITH(flint,
      AS_HELP_STRING([--with-flint=@<:@=DIR@:>@], [Flint install directory]), [
         CPPFLAGS="$CPPFLAGS -I$withval/include"
         LDFLAGS="$LDFLAGS -L$withval/lib"
      ])

AC_CHECK_LIB(flint, __flint_get_memory_functions, , [AC_MSG_ERROR(
   [ERROR FLINT])])

LIBS="$LIBS -lm -lpthread"


#################################################
AS_BOX([          mapml config               ])
#################################################


# Check maple
#------------

AC_MSG_CHECKING([use of maple])

AC_ARG_WITH(maplec,
      AS_HELP_STRING([--with-maplec=@<:@=DIR@:>@], [Maple include]), [
         MAPLE_INCDIR="$withval/include"])


AS_IF([test "x$with_maplec" != "x"], [
      AS_ECHO("x$with_maplec")
      MAPLE_CPPFLAGS="$CPPFLAGS -I${MAPLE_INCDIR}"
      AC_SUBST(MAPLE_CPPFLAGS)
      AC_CHECK_HEADER([$MAPLE_INCDIR/maplec.h],
                [AC_DEFINE([WITH_MAPLEC], [1],
                   [with maplec.h])],
                [AC_MSG_ERROR([Sorry, maplec.h location missing])])
      AC_ARG_WITH(maplelib,
            AS_HELP_STRING([--with-maplelib=@<:@=DIR@:>@], [Maple lib]), [
            MAPLE_LIBDIR="$withval"
            ])
      MAPLE_LDFLAGS="$LDFLAGS -L${MAPLE_LIBDIR}"
      AC_SUBST(MAPLE_LDFLAGS)
      MAPLE_LIBS="$LIBS -lmaplec"
      AC_SUBST(MAPLE_LIBS)
],[AC_MSG_RESULT(no)])

AM_CONDITIONAL(WITH_MAPLE, test "x$WITH_MAPLEC" = "xyes")


#GV Does not seem to work, howto do ? 

#AC_CHECK_LIB(maplec, __EvalMapleProc, , [AC_MSG_ERROR(
   #[ERROR MAPLELIB])])



################   PML FILES   #################################


AS_ECHO(["##--------------------------------------##"])

AC_CONFIG_LINKS([flint-extras/include/nmod_extra.h:flint-extras/nmod_extra/nmod_extra.h
flint-extras/include/fmpz_extra.h:flint-extras/fmpz_extra/fmpz_extra.h
flint-extras/include/fmpz_mat_extra.h:flint-extras/fmpz_mat_extra/fmpz_mat_extra.h
flint-extras/include/nmod_vec_extra.h:flint-extras/nmod_vec_extra/nmod_vec_extra.h
flint-extras/include/nmod_mat_extra.h:flint-extras/nmod_mat_extra/nmod_mat_extra.h
flint-extras/include/nmod_poly_extra.h:flint-extras/nmod_poly_extra/nmod_poly_extra.h
flint-extras/include/nmod_poly_fft.h:flint-extras/nmod_poly_extra/nmod_poly_fft.h
flint-extras/include/nmod_mat_poly.h:flint-extras/nmod_mat_poly_extra/nmod_mat_poly.h
flint-extras/include/sagemath_extra.h:flint-extras/sagemath_extra/sagemath_extra.h
flint-extras/include/nmod_poly_mat_approximant.h:flint-extras/nmod_poly_mat_extra/nmod_poly_mat_approximant.h
flint-extras/include/nmod_poly_mat_io.h:flint-extras/nmod_poly_mat_extra/nmod_poly_mat_io.h
flint-extras/include/nmod_poly_mat_arith.h:flint-extras/nmod_poly_mat_extra/nmod_poly_mat_arith.h 
flint-extras/include/nmod_poly_mat_multiply.h:flint-extras/nmod_poly_mat_extra/nmod_poly_mat_multiply.h
flint-extras/include/nmod_poly_mat_extra.h:flint-extras/nmod_poly_mat_extra/nmod_poly_mat_extra.h
flint-extras/include/nmod_poly_mat_utils.h:flint-extras/nmod_poly_mat_extra/nmod_poly_mat_utils.h
flint-extras/include/nmod_poly_mat_forms.h:flint-extras/nmod_poly_mat_extra/nmod_poly_mat_forms.h
])

AC_CONFIG_FILES([
Makefile
flint-extras/Makefile
flint-extras/include/Makefile
flint-extras/lib/Makefile
flint-extras/tests/Makefile
mapml/Makefile
mapml/src/Makefile
mapml/maple/Makefile
])

AC_OUTPUT

