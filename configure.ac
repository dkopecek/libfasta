AC_PREREQ([2.67])
AC_INIT([libfasta], [0.1], [xkopecek@fi.muni.cz])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign])

# If automake supports "silent rules", enable them by default
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# Checks for programs.
AC_PROG_CC
AM_PROG_LIBTOOL
AM_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL

CFLAGS="$CFLAGS -O0 -g -fno-inline-functions -pipe -std=c99 -Wall -Wundef -Wno-unknown-pragmas -Wno-unused"

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h stddef.h stdint.h stdlib.h string.h sys/file.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([atexit strchr strdup strerror])

AC_CONFIG_FILES([Makefile
		 src/Makefile])

AC_OUTPUT
