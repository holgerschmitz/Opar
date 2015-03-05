# configure.ac
#
# Created on: 30 Jul 2012
# Author: Holger Schmitz
# Email: holger@notjustphysics.com
#
# Copyright 2012 Holger Schmitz
#
# This file is part of OPar.
#
# OPar is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OPar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OPar.  If not, see <http://www.gnu.org/licenses/>.
 
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([OPar],[0.9.1])

AC_CONFIG_SRCDIR([src/opar.cpp])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([src/config.hpp])

AM_INIT_AUTOMAKE
VERSION="0.9.1"
PACKAGE="opar"

AC_LANG(C++)

AC_SUBST(INTI_CFLAGS)
AC_SUBST(INTI_LIBS)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

LT_INIT
AC_SUBST([LIBTOOL_DEPS])

# AC_ARG_PROGRAM

# Checks for libraries.

AX_BOOST_BASE([1.48.0],,[AC_MSG_ERROR([OPar needs Boost version 1.48.0 to compile!])])

CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"

# LDFLAGS="$LDFLAGS -lschnek"

# Search for MPI support
# allows the  --with-mpi option (yes, no, auto)
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
           [compile with MPI (parallelization) support. If none is found,
           MPI is not used. Default: auto])
       ],,[with_mpi=auto])
       
# Search for MPI support
# allows the  --with-mpi option (yes, no, auto)
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
           [compile with MPI (parallelization) support. If none is found,
           MPI is not used. Default: auto])
       ],,[with_mpi=auto])
       
if test x"$with_mpi" != xno; then
  AX_MPI([
    using_mpi=yes
    CXX=$MPICXX
    CC=$MPICXX
    CPP=$MPICXX
    AC_DEFINE([HAVE_MPI],[1],[Define this macro if Opar is compiled with MPI support])
  ],[
    using_mpi=no
    if test x"$with_mpi" = xyes; then
      AC_MSG_ERROR([MPI support requested but no MPI library found!])
    fi
  ])
else
  using_mpi=no
fi

AX_HDF5([
  have_hdf5="yes"
  EXTRA_CPPFLAGS="$EXTRA_CPPFLAGS $HDF5_CFLAGS"
  EXTRA_LDFLAGS="$EXTRA_LDFLAGS $HDF5_LDFLAGS"
  EXTRA_LIBS="$EXTRA_LIBS $HDF5_LIBS"
],[
  have_hdf5="no"
])

#End :: Check HDF5 support

AX_SCHNEK([0.9.1],[
  EXTRA_CPPFLAGS="$EXTRA_CPPFLAGS $SCHNEK_CPPFLAGS"
  EXTRA_LDFLAGS="$EXTRA_LDFLAGS $SCHNEK_LDFLAGS"
  EXTRA_LIBS="$EXTRA_LIBS $SCHNEK_LIBS"
],[
  AC_MSG_ERROR([Schnek library needed by OPar but could not be found!])
])

CPPFLAGS="$CPPFLAGS $EXTRA_CPPFLAGS"
LDFLAGS="$LDFLAGS $EXTRA_LDFLAGS"
LIBS="$EXTRA_LIBS"

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_HEADER_STDC

AC_ARG_WITH([dimension],
            AC_HELP_STRING([--with-dimension=dim],
                           [The dimensionality of the code.
                            dim can be either 1,2 or 3.
                            Default is 1]))
            
if test x"$with_dimension" = x"1" ; then
  AC_DEFINE([ONE_DIMENSIONAL],[1],[Defined if the OPar is compiled in 1D configuration])
elif test x"$with_dimension" = x"2" ; then
  AC_DEFINE([TWO_DIMENSIONAL],[1],[Defined if the OPar is compiled in 2D configuration])
elif test x"$with_dimension" = x"3" ; then
  AC_DEFINE([THREE_DIMENSIONAL],[1],[Defined if the OPar is compiled in 3D configuration])
else
  AC_MSG_NOTICE([Invalid dimension! Defaulting to 1.])
  AC_DEFINE([ONE_DIMENSIONAL])
fi



AC_CONFIG_FILES([Makefile
                 src/Makefile
                 checks/Makefile])
AC_OUTPUT
