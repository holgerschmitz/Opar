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
  AC_DEFINE([ONE_DIMENSIONAL])
elif test x"$with_dimension" = x"2" ; then
  AC_DEFINE([TWO_DIMENSIONAL])
elif test x"$with_dimension" = x"3" ; then
  AC_DEFINE([THREE_DIMENSIONAL])
else
  AC_MSG_NOTICE([Invalid dimension! Defaulting to 1.])
  AC_DEFINE([ONE_DIMENSIONAL])
fi


AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
