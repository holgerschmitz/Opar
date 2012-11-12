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

AC_INIT([OPar],[0.1])

AC_CONFIG_SRCDIR([src/opar.cpp])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([src/config.hpp])
AM_INIT_AUTOMAKE
# VERSION="0.1"
# PACKAGE="opar"

AC_LANG(C++)

AC_SUBST(INTI_CFLAGS)
AC_SUBST(INTI_LIBS)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET


# AC_ARG_PROGRAM

# Checks for libraries.

AX_BOOST_BASE([1.48.0],,[AC_MSG_ERROR([Schnek needs Boost version 1.48.0 to compile!])])

# LDFLAGS="$LDFLAGS -lschnek"

# Search for MPI support
# allows the  --with-mpi option (yes, no, auto)
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
           [compile with MPI (parallelization) support. If none is found,
           MPI is not used. Default: auto])
       ],,[with_mpi=auto])
       
if test x"$with_mpi" != xno; then
#  CXXFLAGS_save=$CXXFLAGS
#  LDFLAGS_save=$LDFLAGS
#  CXXFLAGS=""
#  LDFLAGS=""
  AX_MPI([
#    MPI_CXXFLAGS=`$MPICXX --showme:compile`
#    [echo Using MPI_CXXFLAGS=$MPI_CXXFLAGS]
#    AC_SUBST(MPI_CXXFLAGS)
#    MPI_LDFLAGS=`$MPICXX --showme:link`
#    [echo Using MPI_LDFLAGS=$MPI_LDFLAGS]
#    AC_SUBST(MPI_LDFLAGS)
#    CXXFLAGS="$CXXFLAGS_save $MPI_CXXFLAGS"
#    LDFLAGS="$LDFLAGS_save $MPI_LDFLAGS"
     CXX="mpic++"
    AC_DEFINE([HAVE_MPI],[1],[Define this macro if Schnek is compiled with MPI support])
  ],[
    if test x"$with_mpi" = xyes; then
      AC_MSG_ERROR([MPI support requested but no MPI library found!])
    fi
  ])
fi


# Search for HDF5 libraries
# To distinguish between serial and parallel, use the H5_HAVE_PARALLEL macro
# supplied by HDF5 itself
# Checks for libraries.
# start: Check for HDF5 support
AC_ARG_WITH(hdf5, [AC_HELP_STRING([--without-hdf5],[do not allow HDF5 output])], ok=$withval, ok=yes)
if test "$ok" = "yes"; then
  AC_CHECK_LIB(z,deflate, [],[AC_MSG_WARN([zlib is required for HDF5!])])

  save_LIBS_0="$LIBS" # save, to check later if we found any library
  AC_CHECK_LIB(hdf5, H5Pcreate, [
    AC_CHECK_HEADERS(hdf5.h, [LIBS="-lhdf5 $LIBS"
                     AC_DEFINE([SCHNEK_HAVE_HDF5],[1],[Define this macro if Schnek is compiled with HDF5 support])])])
  if test x"$save_LIBS_0" = x"$LIBS"; then
    AC_MSG_WARN([Couldn't find the HDF5 library!!  Switching to --without-hdf5.])
  fi

  if test "x$with_mpi" = "xyes"; then
    AC_CHECK_FUNCS(H5Pset_mpi H5Pset_fapl_mpio)
  fi
fi

AC_ARG_ENABLE([parallel-hdf5],
          [AS_HELP_STRING([--enable-parallel-hdf5], [use parallel file access for HDF5 (mpio) if available.])], 
          [],[])

if test x"$enable-parallel-hdf5" != "xno" ; then
  AC_DEFINE([SCHNEK_USE_HDF_PARALLEL],[1],[Defined if we use parallel file access for HDF5 (mpio) if available.])
fi

AC_LANG_PUSH([C++])
OLDCPPFLAGS=$CPPFLAGS
OLDLDFLAGS=$LDFLAGS
CPPFLAGS=$hdf5_CFLAGS
LDFLAGS=$hdf5_LIBS

CPPFLAGS=$OLDCPPFLAGS
LDFLAGS=$OLDLDFLAGS
AC_LANG_POP([C++])
#End :: Check HDF5 support





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
                 src/Makefile])
AC_OUTPUT
