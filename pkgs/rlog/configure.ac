dnl Process this file with autoconf to produce a configure script.

AC_INIT
AC_CONFIG_SRCDIR([rlog/rlog.h]) dnl a source file from your sub dir
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

dnl This ksh/zsh feature conflicts with `cd blah ; pwd`
unset CDPATH

dnl Checking host/target/build systems, for make, install etc.
AC_CANONICAL_TARGET
dnl Perform program name transformation
AC_ARG_PROGRAM

AC_LANG(C++)

# Must be updated for public releases..
# CURRENT : REVISON : AGE
#
# +1 : ? : +1  =>  new interface that does not break old one
# +1 : 0 : 0   =>  new interface that breaks old one
#  ? : 0 : 0   =>  no new interfaces, but breaks apps
#  ? :+1 : ?   =>  internal changes, nothing breaks, and might work better..
#
LIBRLOG_VERSION=5:0:0

AC_SUBST(LIBRLOG_VERSION)
AC_SUBST(PACKAGE_VERSION)

AM_INIT_AUTOMAKE(rlog, 1.4) dnl searches for some needed programs

dnl
unset CDPATH
dnl make /usr/local the default for the installation
AC_PREFIX_DEFAULT(/usr/local)

if test "x$prefix" = "xNONE"; then
  prefix=$ac_default_prefix
  ac_configure_args="$ac_configure_args --prefix $prefix"
fi

dnl without this order in this file, automake will be confused!
dnl
AC_CONFIG_HEADERS([config.h])

dnl create only shared libtool-libraries
AC_ENABLE_SHARED(yes)

dnl set the following to yes, if you want to create static
dnl libtool-libraries, too.
AC_ENABLE_STATIC(no)

dnl checks for programs.
dnl first check for c/c++ compilers
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_LIBTOOL

dnl CXXFLAGS="$NOOPT_CXXFLAGS" dnl __kdevelop[noopt]__
dnl CFLAGS="$NOOPT_CFLAGS" dnl __kdevelop[noopt]__
dnl CXXFLAGS="$CXXFLAGS $USE_EXCEPTIONS" dnl __kdevelop[exc]__

ACX_PTHREAD

AC_ARG_ENABLE(printffp,
    AS_HELP_STRING([--disable-printffp],
	    [disables printf attribute on function pointers]),
    with_printffp=$enableval,with_printffp="test")

if test "x$with_printffp" = "xtest"; then

    # check if we can use printf attribute on a function..
    AC_MSG_CHECKING(if __printf__ attribute can apply to function pointers)
    AC_COMPILE_IFELSE( [[
	    void testfunc(const char *format, ...)
	    { }
	    void (*funcProto)(const char *, ...) 
		__attribute__ (( __format__ ( __printf__ , 1, 2 ))); 

	    int main()
	    {
		funcProto = testfunc;
		(*funcProto)("hello %s", "world");

		return 0;
	    }
	 ]],
	with_printffp="yes",
	with_printffp="no")

    AC_MSG_RESULT($with_printffp)
fi

if test "x$with_printffp" = "xyes"; then
    HAVE_PRINTF_FP_PROTOTYPE="1"
else
    HAVE_PRINTF_FP_PROTOTYPE="0"
fi
AC_SUBST(HAVE_PRINTF_FP_PROTOTYPE)

AC_ARG_ENABLE(vararg,
    AS_HELP_STRING([--disable-vararg],
	    [don't use vararg macros even if the compiler supports them]),
    enable_vararg=$enableval)
if test "x$enable_vararg" != "xno"; then
    # check if the compiler understands __VA_ARGS__ 
AC_MSG_CHECKING(if compiler has C99 variadac macro)
AC_COMPILE_IFELSE( [[
#include <stdio.h>

#define PRINT(FMT, ...) printf( FMT, __VA_ARGS__ )

int main()
{
    PRINT("hello %s", "world");
    return 0;
} ]], [C99_VARIADAC_MACROS="1" 
       AC_MSG_RESULT(yes)], 
    [C99_VARIADAC_MACROS="0"
     AC_MSG_RESULT(no)])


# check if the compiler understands pre-c99 variadac macros
AC_MSG_CHECKING(if compiler has pre-C99 variadac macro)
AC_COMPILE_IFELSE( [[
#include <stdio.h>

#define PRINT(FMT, ARGS...) printf( FMT, ##ARGS )

int main()
{
    PRINT("hello %s", "world");
    return 0;
} ]], [PREC99_VARIADAC_MACROS="1" 
       AC_MSG_RESULT(yes)], 
    [PREC99_VARIADAC_MACROS="0"
     AC_MSG_RESULT(no)])
else
dnl vararg support disabled
    C99_VARIADAC_MACROS=0
    PREC99_VARIADAC_MACROS=0
fi # enable_vararg

AC_SUBST(C99_VARIADAC_MACROS)
AC_SUBST(PREC99_VARIADAC_MACROS)


# check if the computer has hardware cycle counter
AC_ARG_ENABLE(rdtsc,
    AS_HELP_STRING([--enable-rdtsc],
	    [use rdtsc for benchmarking rlog in test programs]),
    with_rdtsc=$enableval,with_rdtsc="no")
if test "x$with_rdtsc" = "xyes"; then
    AC_DEFINE([USE_RDTSC], [1], [Define to use RDTSC calls for benchmark test])
else
    AC_DEFINE([USE_RDTSC], [0], [Define to use RDTSC calls for benchmark test])
fi

AC_ARG_ENABLE(valgrind,
    AS_HELP_STRING([--disable-valgrind],
		   [disables valgrind support code.]),
    with_valgrind=$enableval, with_valgrind="yes" )
USE_VALGRIND="0"
if test "x$with_valgrind" = "xyes"; then
    AC_CHECK_HEADER([valgrind/valgrind.h],
	AC_CHECK_DECLS([VALGRIND_PRINTF_BACKTRACE],
	    USE_VALGRIND="1",,[#include <valgrind/valgrind.h>]))
fi
AC_SUBST(USE_VALGRIND)


# allow documentation build to be disabled manually
AC_ARG_ENABLE(docs,
    AS_HELP_STRING([--disable-docs],[disable documentation build]),
    build_docs=$enableval, build_docs="yes")

# check for tools necessary to build documentation
AC_PATH_PROG(DOXYGEN, doxygen, [no])
AC_PATH_PROG(LATEX, latex, [no])
AC_PATH_PROG(PDFLATEX, pdflatex, [no])

AM_CONDITIONAL(BUILD_DOCS, test "x$build_docs" = "xyes" \
			     -a "x$DOXYGEN" != "xno" \
			     -a "x$LATEX" != "xno" \
			     -a "x$PDFLATEX" != "xno" )

AC_CHECK_FUNCS(localtime_r)

# check for files which may not be available if using old STL code
AC_CHECK_HEADERS(sstream)

LDFLAGS="$LDFLAGS $PTHREAD_LIBS $USER_LDFLAGS"
CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS $USER_INCLUDES"

AC_SUBST(AUTODIRS)
AC_CONFIG_FILES([Makefile \
	rlog/Makefile \
	docs/Makefile \
	config/Makefile \
	rlog/common.h \
	librlog.pc \
	rlog.spec \
	makedist2.sh])
AC_OUTPUT

