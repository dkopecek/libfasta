AC_PREREQ([2.50])
AC_INIT([libfasta], [1.5], [dnk1618@gmail.com])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign subdir-objects])
LT_INIT

# If automake supports "silent rules", enable them by default
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# Checks for programs.
AC_PROG_CC
AM_PROG_LIBTOOL
AM_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL

# libtool versioning
# See http://sources.redhat.com/autobook/autobook/autobook_91.html#SEC91 for details

## increment if the interface has additions, changes, removals.
LT_CURRENT=1

## increment any time the source changes; set 0 to if you increment CURRENT
LT_REVISION=0

## increment if any interfaces have been added; set to 0
## if any interfaces have been changed or removed. removal has
## precedence over adding, so set to 0 if both happened.
LT_AGE=0

AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

CFLAGS="$CFLAGS -pipe -std=c99 -Wall -Wextra -Wundef -Wno-unknown-pragmas -Wformat"
CFLAGS="$CFLAGS -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE=1"

CFLAGS_DEBUGGING="-O0 -g -fno-inline-functions"

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h limits.h stdlib.h unistd.h ctype.h errno.h stdbool.h sys/stat.h assert.h stddef.h string.h sys/types.h stdio.h stdint.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_CHECK_FUNCS([malloc realloc atexit strchr strdup strerror lseek])

AC_ARG_ENABLE([debug],
     [AC_HELP_STRING([--enable-debug], [enable debugging flags (default=no)])],
     [case "${enableval}" in
       yes) debug=yes ;;
       no)  debug=no ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
     esac], [debug=no])

if test "$debug" = "yes"; then
   CFLAGS="$CFLAGS $CFLAGS_DEBUGGING"
else
   AC_DEFINE([NDEBUG], [1], [Compile without debugging code])
fi

AC_CONFIG_FILES([Makefile
		 src/Makefile
		 tests/Makefile
		 libfasta.pc])

AC_OUTPUT

echo "CFLAGS: $CFLAGS"
