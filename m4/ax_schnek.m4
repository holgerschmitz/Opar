
AC_DEFUN([AX_SCHNEK],
[
AC_ARG_WITH([schnek],
  [AS_HELP_STRING([--with-schnek@<:@=ARG@:>@],
    [use Schnek library from a standard location (ARG=yes),
     from the specified location (ARG=<path>),
     or disable it (ARG=no)
     @<:@ARG=yes@:>@ ])],
    [
    if test "$withval" = "no"; then
        want_schnek="no"
    elif test "$withval" = "yes"; then
        want_schnek="yes"
        ac_schnek_path=""
    else
        want_schnek="yes"
        ac_schnek_path="$withval"
    fi
    ],
    [want_schnek="yes"])


AC_ARG_WITH([schnek-libdir],
        AS_HELP_STRING([--with-schnek-libdir=LIB_DIR],
        [Force given directory for Schnek libraries. Note that this will override library path detection, so use this parameter only if default library detection fails and you know exactly where your Schnek libraries are located.]),
        [
        if test -d "$withval"
        then
                ac_schnek_lib_path="$withval"
        else
                AC_MSG_ERROR(--with-schnek-libdir expected directory name)
        fi
        ],
        [ac_schnek_lib_path=""]
)

if test "x$want_schnek" = "xyes"; then
    schnek_lib_version_req=ifelse([$1], ,0.9.0,$1)
    schnek_lib_version_req_shorten=`expr $schnek_lib_version_req : '\([[0-9]]*\.[[0-9]]*\)'`
    schnek_lib_version_req_major=`expr $schnek_lib_version_req : '\([[0-9]]*\)'`
    schnek_lib_version_req_minor=`expr $schnek_lib_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
    schnek_lib_version_req_sub_minor=`expr $schnek_lib_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
    if test "x$schnek_lib_version_req_sub_minor" = "x" ; then
        schnek_lib_version_req_sub_minor="0"
        fi
    WANT_SCHNEK_VERSION=`expr $schnek_lib_version_req_major \* 100000 \+  $schnek_lib_version_req_minor \* 100 \+ $schnek_lib_version_req_sub_minor`
    AC_MSG_CHECKING(for Schnek version >= $schnek_lib_version_req)
    succeeded=no
    
    libsubdirs="lib"
    ax_arch=`uname -m`
    if test $ax_arch = x86_64 -o $ax_arch = ppc64 -o $ax_arch = s390x -o $ax_arch = sparc64; then
        libsubdirs="lib64 lib lib64"
    fi

    if test "$ac_schnek_path" != ""; then
        SCHNEK_CPPFLAGS="-I$ac_schnek_path/include"
        for ac_schnek_path_tmp in $libsubdirs; do
                if test -d "$ac_schnek_path"/"$ac_schnek_path_tmp" ; then
                        SCHNEK_LDFLAGS="-L$ac_schnek_path/$ac_schnek_path_tmp"
                        break
                fi
        done
    else
        if test "$SCHNEK_INCLUDEDIR" != ""; then
            SCHNEK_CPPFLAGS="-I$SCHNEK_INCLUDEDIR"
        else
            for ac_schnek_path_tmp in /usr /usr/local /opt /opt/local ; do
                if test -e "$ac_schnek_path_tmp/include/schnek/grid.hpp"; then
                    SCHNEK_CPPFLAGS="-I$ac_schnek_path_tmp/include"
                    break;
                fi
            done
        fi
        
        if test "$SCHNEK_LIBDIR" != ""; then
            SCHNEK_LDFLAGS="-L$HDF5_LIBDIR"
        fi
    fi

    dnl overwrite ld flags if we have required special directory with
    dnl --with-schnek-libdir parameter
    if test "$ac_schnek_lib_path" != ""; then
       SCHNEK_LDFLAGS="-L$ac_schnek_lib_path"
    fi

    SCHNEK_LIBS="-lschnek"
    
    CPPFLAGS_SAVED="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $SCHNEK_CPPFLAGS"
    export CPPFLAGS

    LDFLAGS_SAVED="$LDFLAGS"
    LDFLAGS="$LDFLAGS $SCHNEK_LDFLAGS"
    export LDFLAGS

    AC_REQUIRE([AC_PROG_CXX])
    
    AC_LANG_PUSH(C++)
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <schnek/schnek_config.hpp>
    ]], [[ 
      //schnek::Grid<double,2> grid();
      #if SCHNEK_VERSION >= $WANT_SCHNEK_VERSION
      // Everything is okay
      #else
      #  error Schnek version is too old
      #endif
    ]])],[
        AC_MSG_RESULT(yes)
    succeeded=yes
    found_system=yes
        ],[
        AC_MSG_RESULT(no)
        ])
    
    AC_LANG_POP([C++])


    if test "$succeeded" != "yes" ; then
        AC_MSG_NOTICE([[We could not detect the Schnek libraries (version $schnek_lib_version_req or higher).]])
        # execute ACTION-IF-NOT-FOUND (if present):
        ifelse([$3], , :, [$3])
    else
        AC_SUBST(SCHNEK_CPPFLAGS)
        AC_SUBST(SCHNEK_LIBS)
        AC_SUBST(SCHNEK_LDFLAGS)
        AC_DEFINE(HAVE_SCHNEK,,[define if the Schnek library is available])
        # execute ACTION-IF-FOUND (if present):
        ifelse([$2], , :, [$2])
    fi

    CPPFLAGS="$CPPFLAGS_SAVED"
    LDFLAGS="$LDFLAGS_SAVED"
fi

])
