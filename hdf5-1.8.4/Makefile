# Makefile.in generated by automake 1.11 from Makefile.am.
# Makefile.  Generated from Makefile.in by configure.

# Copyright (C) 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
# 2003, 2004, 2005, 2006, 2007, 2008, 2009  Free Software Foundation,
# Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.



#
# Copyright by The HDF Group.
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the files COPYING and Copyright.html.  COPYING can be found at the root
# of the source code distribution tree; Copyright.html can be found at the
# root level of an installed copy of the electronic HDF5 document set and
# is linked from the top-level documents page.  It can also be found at
# http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have
# access to either file, you may request a copy from help@hdfgroup.org.
# 
# This makefile mostly just reinvokes make in the various subdirectories
# but does so in the correct order.  You can alternatively invoke make from
# each subdirectory manually.
#
# Top-level HDF5 Makefile(.in)

# pmake has issues if variables are undefined.  Solve this problem in
# top-level Makefile by defining .MAKEFLAGS target to -V before pmake can
# encounter any undefined variables.
# Automake resists putting anything but variable definitions first in
# a Makefile.in, so we'll put a placebo comment here and use sed in
# bin/reconfigure to turn it into the .MAKEFLAGS target.  Sigh.  -JL 2005
# Configure should set AM_MAKEFLAGS to -V to solve this problem in
# subdirectories.
# NOTE: This means that invoking pmake in a subdirectory will not work.
.MAKEFLAGS:
	$(MAKE) all


pkgdatadir = $(datadir)/hdf5
pkgincludedir = $(includedir)/hdf5
pkglibdir = $(libdir)/hdf5
pkglibexecdir = $(libexecdir)/hdf5
am__cd = CDPATH="$${ZSH_VERSION+.}$(PATH_SEPARATOR)" && cd
install_sh_DATA = $(install_sh) -c -m 644
install_sh_PROGRAM = $(install_sh) -c
install_sh_SCRIPT = $(install_sh) -c
INSTALL_HEADER = $(INSTALL_DATA)
transform = $(program_transform_name)
NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
build_triplet = x86_64-unknown-linux-gnu
host_triplet = x86_64-unknown-linux-gnu
DIST_COMMON = $(am__configure_deps) $(srcdir)/Makefile.am \
	$(srcdir)/Makefile.in $(top_srcdir)/config/commence.am \
	$(top_srcdir)/configure COPYING bin/COPYING bin/compile \
	bin/config.guess bin/config.sub bin/depcomp bin/install-sh \
	bin/ltmain.sh bin/missing bin/mkinstalldirs
subdir = .
ACLOCAL_M4 = $(top_srcdir)/aclocal.m4
am__aclocal_m4_deps = $(top_srcdir)/configure.in
am__configure_deps = $(am__aclocal_m4_deps) $(CONFIGURE_DEPENDENCIES) \
	$(ACLOCAL_M4)
am__CONFIG_DISTCLEAN_FILES = config.status config.cache config.log \
 configure.lineno config.status.lineno
mkinstalldirs = $(SHELL) $(top_srcdir)/bin/mkinstalldirs
CONFIG_HEADER = $(top_builddir)/src/H5config.h
CONFIG_CLEAN_FILES =
CONFIG_CLEAN_VPATH_FILES =
SOURCES =
DIST_SOURCES =
RECURSIVE_TARGETS = all-recursive check-recursive dvi-recursive \
	html-recursive info-recursive install-data-recursive \
	install-dvi-recursive install-exec-recursive \
	install-html-recursive install-info-recursive \
	install-pdf-recursive install-ps-recursive install-recursive \
	installcheck-recursive installdirs-recursive pdf-recursive \
	ps-recursive uninstall-recursive
RECURSIVE_CLEAN_TARGETS = mostlyclean-recursive clean-recursive	\
  distclean-recursive maintainer-clean-recursive
AM_RECURSIVE_TARGETS = $(RECURSIVE_TARGETS:-recursive=) \
	$(RECURSIVE_CLEAN_TARGETS:-recursive=) tags TAGS ctags CTAGS \
	distdir dist dist-all distcheck
ETAGS = etags
CTAGS = ctags
DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
distdir = $(PACKAGE)-$(VERSION)
top_distdir = $(distdir)
am__remove_distdir = \
  { test ! -d "$(distdir)" \
    || { find "$(distdir)" -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr "$(distdir)"; }; }
am__relativize = \
  dir0=`pwd`; \
  sed_first='s,^\([^/]*\)/.*$$,\1,'; \
  sed_rest='s,^[^/]*/*,,'; \
  sed_last='s,^.*/\([^/]*\)$$,\1,'; \
  sed_butlast='s,/*[^/]*$$,,'; \
  while test -n "$$dir1"; do \
    first=`echo "$$dir1" | sed -e "$$sed_first"`; \
    if test "$$first" != "."; then \
      if test "$$first" = ".."; then \
        dir2=`echo "$$dir0" | sed -e "$$sed_last"`/"$$dir2"; \
        dir0=`echo "$$dir0" | sed -e "$$sed_butlast"`; \
      else \
        first2=`echo "$$dir2" | sed -e "$$sed_first"`; \
        if test "$$first2" = "$$first"; then \
          dir2=`echo "$$dir2" | sed -e "$$sed_rest"`; \
        else \
          dir2="../$$dir2"; \
        fi; \
        dir0="$$dir0"/"$$first"; \
      fi; \
    fi; \
    dir1=`echo "$$dir1" | sed -e "$$sed_rest"`; \
  done; \
  reldir="$$dir2"
DIST_ARCHIVES = $(distdir).tar.gz
GZIP_ENV = --best
distuninstallcheck_listfiles = find . -type f -print
distcleancheck_listfiles = find . -type f -print
ACLOCAL = /home1/packages/automake/automake-1.9.6/bin/aclocal-1.9 -I /afs/ncsa/projects/hdf/packages/libtool_1.5.14/Linux_2.4/share/aclocal
ADD_PARALLEL_FILES = no
AMTAR = ${SHELL} /home/chaako/gitmerge/sceptic3D/hdf5-1.8.4/bin/missing --run tar

# H5_CFLAGS holds flags that should be used when building hdf5,
# but which should not be exported to h5cc for building other programs.
# AM_CFLAGS is an automake construct which should be used by Makefiles 
# instead of CFLAGS, as CFLAGS is reserved solely for the user to define.
AM_CFLAGS =  -std=c99 -pedantic -Wall -Wextra -Wundef -Wshadow -Wpointer-arith -Wbad-function-cast -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Winline -Wno-long-long -Wfloat-equal -Wmissing-format-attribute -Wmissing-noreturn -Wpacked -Wdisabled-optimization -Wformat=2 -Wunreachable-code -Wendif-labels -Wdeclaration-after-statement -Wold-style-definition -Winvalid-pch -Wvariadic-macros -Wnonnull -Winit-self -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wunused-macros -Wunsafe-loop-optimizations -Wc++-compat -Wvolatile-register-var -O -fomit-frame-pointer -finline-functions
AM_CPPFLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_BSD_SOURCE  -D_POSIX_SOURCE   -DNDEBUG -UH5_DEBUG_API
AM_CXXFLAGS =   
AM_FCFLAGS =    
AM_LDFLAGS = 
AM_MAKEFLAGS = 
AR = ar

# Set the paths for AFS installs of autotools for Linux machines
# Ideally, these tools should never be needed during the build.
AUTOCONF = /home1/packages/autoconf/autoconf-2.60/bin/autoconf
AUTOHEADER = /home1/packages/autoconf/autoconf-2.60/bin/autoheader
AUTOMAKE = /home1/packages/automake/automake-1.9.6/bin/automake-1.9
AWK = gawk
BYTESEX = little-endian
CC = gcc
CCDEPMODE = depmode=gcc3
CC_VERSION = /usr/bin/gcc (gcc-4.1.2)
CFLAGS = 
CLEARFILEBUF = yes
CODESTACK = no
CONFIG_DATE = Sat Sep 10 23:09:00 EDT 2011
CONFIG_MODE = production
CONFIG_USER = chaako@sceptic.psfc.mit.edu
CPP = gcc -E
CPPFLAGS = 
CXX = no
CXXCPP = g++ -E
CXXDEPMODE = depmode=gcc3
CXXFLAGS = 
CYGPATH_W = echo
DEBUG_PKG = 
DEFAULT_API_VERSION = v18
DEFS = -DHAVE_CONFIG_H
DEPDIR = .deps
DEPRECATED_SYMBOLS = yes
DIRECT_VFD = no
DSYMUTIL = 
DUMPBIN = 
DYNAMIC_DIRS = 
ECHO_C = 
ECHO_N = -n
ECHO_T = 
EGREP = /bin/grep -E
EXEEXT = 
EXTERNAL_FILTERS = deflate(zlib)

# Make sure that these variables are exported to the Makefiles
F9XMODEXT = mod
F9XMODFLAG = -I
F9XSUFFIXFLAG = 
FC = f95
FCFLAGS = 
FCFLAGS_f90 = 
FCLIBS =  -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lgfortranbegin -lgfortran -lm
FGREP = /bin/grep -F
FILTERS = shuffle,fletcher32,nbit,scaleoffset
FSEARCH_DIRS = 
GPFS = no
GREP = /bin/grep
H5_CFLAGS = -std=c99 -pedantic -Wall -Wextra -Wundef -Wshadow -Wpointer-arith -Wbad-function-cast -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Winline -Wno-long-long -Wfloat-equal -Wmissing-format-attribute -Wmissing-noreturn -Wpacked -Wdisabled-optimization -Wformat=2 -Wunreachable-code -Wendif-labels -Wdeclaration-after-statement -Wold-style-definition -Winvalid-pch -Wvariadic-macros -Wnonnull -Winit-self -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wunused-macros -Wunsafe-loop-optimizations -Wc++-compat -Wvolatile-register-var -O -fomit-frame-pointer -finline-functions
H5_CPPFLAGS = -D_POSIX_SOURCE   -DNDEBUG -UH5_DEBUG_API
H5_CXXFLAGS =  
H5_FCFLAGS =  
H5_LONE_COLON = yes
H5_VERSION = 1.8.4
HADDR_T = 
HAVE_DMALLOC = no
HDF5_HL = yes
HDF5_INTERFACES =  fortran
HDF_CXX = no
HDF_FORTRAN = yes
HID_T = 
HL = hl
HL_FOR = 
HSIZE_T = 
HSSIZE_T = 
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_PROGRAM = ${INSTALL}
INSTALL_SCRIPT = ${INSTALL}
INSTALL_STRIP_PROGRAM = $(install_sh) -c -s
INSTRUMENT = no
INSTRUMENT_LIBRARY = no
LD = /usr/bin/ld -m elf_x86_64
LDFLAGS = 
LIBOBJS = 
LIBS =  -lz -lm 
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LINUX_LFS = yes
LIPO = 
LL_PATH = 
LN_S = ln -s
LTLIBOBJS = 
LT_STATIC_EXEC = 
MAINT = #
MAKEINFO = ${SHELL} /home/chaako/gitmerge/sceptic3D/hdf5-1.8.4/bin/missing --run makeinfo
MKDIR_P = /bin/mkdir -p
MPE = no
MPI_GET_SIZE = 
NM = /usr/bin/nm -B
NMEDIT = 
OBJDUMP = objdump
OBJECT_NAMELEN_DEFAULT_F = 
OBJEXT = o
OTOOL = 
OTOOL64 = 
PACKAGE = hdf5
PACKAGE_BUGREPORT = help@hdfgroup.org
PACKAGE_NAME = HDF5
PACKAGE_STRING = HDF5 1.8.4
PACKAGE_TARNAME = hdf5
PACKAGE_URL = 
PACKAGE_VERSION = 1.8.4
PARALLEL = no
PATH_SEPARATOR = :
PERL = perl
PTHREAD = 
RANLIB = ranlib
ROOT = /home/chaako/gitmerge/sceptic3D/hdf5-1.8.4
RUNPARALLEL = 
RUNSERIAL = 
R_INTEGER = 
R_LARGE = 
SEARCH = VPATH=$(srcdir) $(top_builddir)/src $(top_srcdir)/src
SED = /bin/sed
SETX = set -x
SET_MAKE = 
SHELL = /bin/sh
SIZE_T = 
STATIC_EXEC = no
STATIC_SHARED = 
STRICT_FORMAT_CHECKS = no
STRIP = strip
TESTPARALLEL = 
THREADSAFE = no
TIME = time
TR = /usr/bin/tr
TRACE_API = no
UNAME_INFO = Linux sceptic.psfc.mit.edu 2.6.18-194.26.1.el5 #1 SMP Fri Oct 29 14:21:16 EDT 2010 x86_64 x86_64 x86_64 GNU/Linux
USE_FILTER_DEFLATE = yes
USE_FILTER_FLETCHER32 = yes
USE_FILTER_NBIT = yes
USE_FILTER_SCALEOFFSET = yes
USE_FILTER_SHUFFLE = yes
USE_FILTER_SZIP = no
USINGMEMCHECKER = no
VERSION = 1.8.4
WORDS_BIGENDIAN = no
abs_builddir = /home/chaako/gitmerge/sceptic3D/hdf5-1.8.4
abs_srcdir = /home/chaako/gitmerge/sceptic3D/hdf5-1.8.4
abs_top_builddir = /home/chaako/gitmerge/sceptic3D/hdf5-1.8.4
abs_top_srcdir = /home/chaako/gitmerge/sceptic3D/hdf5-1.8.4
ac_ct_CC = gcc
ac_ct_CXX = 
ac_ct_DUMPBIN = 
ac_ct_FC = 
am__include = include
am__leading_dot = .
am__quote = 
am__tar = ${AMTAR} chof - "$$tardir"
am__untar = ${AMTAR} xf -
bindir = ${exec_prefix}/bin
build = x86_64-unknown-linux-gnu
build_alias = 
build_cpu = x86_64
build_os = linux-gnu
build_vendor = unknown
builddir = .
datadir = ${datarootdir}
datarootdir = ${prefix}/share

# Install directories that automake doesn't know about
docdir = $(exec_prefix)/doc
dvidir = ${docdir}
enable_shared = yes
enable_static = yes
exec_prefix = ${prefix}
host = x86_64-unknown-linux-gnu
host_alias = 
host_cpu = x86_64
host_os = linux-gnu
host_vendor = unknown
htmldir = ${docdir}
includedir = ${prefix}/include
infodir = ${datarootdir}/info
install_sh = ${SHELL} /home/chaako/gitmerge/sceptic3D/hdf5-1.8.4/bin/install-sh
libdir = ${exec_prefix}/lib
libexecdir = ${exec_prefix}/libexec
localedir = ${datarootdir}/locale
localstatedir = ${prefix}/var
lt_ECHO = echo
mandir = ${datarootdir}/man
mkdir_p = /bin/mkdir -p
oldincludedir = /usr/include
pdfdir = ${docdir}
prefix = /home/chaako/gitmerge/sceptic3D/hdf5-1.8.4
program_transform_name = s,x,x,
psdir = ${docdir}
sbindir = ${exec_prefix}/sbin
sharedstatedir = ${prefix}/com
srcdir = .
sysconfdir = ${prefix}/etc
target_alias = 
top_build_prefix = 
top_builddir = .
top_srcdir = .

# Shell commands used in Makefiles
RM = rm -f
CP = cp

# Some machines need a command to run executables; this is that command
# so that our tests will run.
# We use RUNTESTS instead of RUNSERIAL directly because it may be that
# some tests need to be run with a different command.  Older versions
# of the makefiles used the command
# $(LIBTOOL) --mode=execute
# in some directories, for instance.
RUNTESTS = $(RUNSERIAL)

# Libraries to link to while building
LIBHDF5 = $(top_builddir)/src/libhdf5.la
LIBH5TEST = $(top_builddir)/test/libh5test.la
LIBH5F = $(top_builddir)/fortran/src/libhdf5_fortran.la
LIBH5FTEST = $(top_builddir)/fortran/test/libh5test_fortran.la
LIBH5CPP = $(top_builddir)/c++/src/libhdf5_cpp.la
LIBH5TOOLS = $(top_builddir)/tools/lib/libh5tools.la
LIBH5_HL = $(top_builddir)/hl/src/libhdf5_hl.la
LIBH5F_HL = $(top_builddir)/hl/fortran/src/libhdf5hl_fortran.la
LIBH5CPP_HL = $(top_builddir)/hl/c++/src/libhdf5_hl_cpp.la

# Scripts used to build examples
# If only shared libraries have been installed, have h5cc build examples with
# shared libraries instead of static libraries
H5CC = $(bindir)/h5cc
H5CC_PP = $(bindir)/h5pcc
H5FC = $(bindir)/h5fc
H5FC_PP = $(bindir)/h5pfc
H5CPP = $(bindir)/h5c++

# The trace script; this is used on source files from the C library to
# insert tracing macros.
TRACE = perl $(top_srcdir)/bin/trace

# .chkexe files are used to mark tests that have run successfully.
# .chklog files are output from those tests.
# *.clog are from the MPE option.

# Some files/directories generated during check that should be cleaned
CHECK_CLEANFILES = *.chkexe *.chklog *.clog *-tmp
TESTPARALLEL_DIR = 

# Define subdirectories to build.
# Add this directory to SUBDIRS so that examples get built after tools
# but before examples in extra interfaces (c++ and fortran).
# Since we're explicitly listing DIST_SUBDIRS, we also need to list
# directories that are only conditionally built (so that their Makefiles
# are cleaned as well).
# Note that `make clean' will not affect the examples or doc directories. 

# Conditionals.  These conditionals are defined during configure
# Define each variable to empty if it is not used to placate pmake
#TESTPARALLEL_DIR = testpar
CXX_DIR = 
#CXX_DIR = c++
#FORTRAN_DIR = 
FORTRAN_DIR = fortran
#HDF5_HL_DIR = 
HDF5_HL_DIR = hl
SUBDIRS = src test $(TESTPARALLEL_DIR) tools . $(CXX_DIR) $(FORTRAN_DIR) \
          $(HDF5_HL_DIR) perform

DIST_SUBDIRS = src test testpar tools . c++ fortran hl perform examples

# Some files generated during configure that should be cleaned
DISTCLEANFILES = config/stamp1 config/stamp2
all: all-recursive

.SUFFIXES:
am--refresh:
	@:
$(srcdir)/Makefile.in: # $(srcdir)/Makefile.am $(top_srcdir)/config/commence.am $(am__configure_deps)
	@for dep in $?; do \
	  case '$(am__configure_deps)' in \
	    *$$dep*) \
	      echo ' cd $(srcdir) && $(AUTOMAKE) --foreign'; \
	      $(am__cd) $(srcdir) && $(AUTOMAKE) --foreign \
		&& exit 0; \
	      exit 1;; \
	  esac; \
	done; \
	echo ' cd $(top_srcdir) && $(AUTOMAKE) --foreign Makefile'; \
	$(am__cd) $(top_srcdir) && \
	  $(AUTOMAKE) --foreign Makefile
.PRECIOUS: Makefile
Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@case '$?' in \
	  *config.status*) \
	    echo ' $(SHELL) ./config.status'; \
	    $(SHELL) ./config.status;; \
	  *) \
	    echo ' cd $(top_builddir) && $(SHELL) ./config.status $@ $(am__depfiles_maybe)'; \
	    cd $(top_builddir) && $(SHELL) ./config.status $@ $(am__depfiles_maybe);; \
	esac;

$(top_srcdir)/configure: # $(am__configure_deps)
	$(am__cd) $(srcdir) && $(AUTOCONF)
$(ACLOCAL_M4): # $(am__aclocal_m4_deps)
	$(am__cd) $(srcdir) && $(ACLOCAL) $(ACLOCAL_AMFLAGS)
$(am__aclocal_m4_deps):

mostlyclean-libtool:
	-rm -f *.lo

clean-libtool:
	-rm -rf .libs _libs

distclean-libtool:
	-rm -f libtool config.lt

# This directory's subdirectories are mostly independent; you can cd
# into them and run `make' without going through this Makefile.
# To change the values of `make' variables: instead of editing Makefiles,
# (1) if the variable is set in `config.status', edit `config.status'
#     (which will cause the Makefiles to be regenerated when you run `make');
# (2) otherwise, pass the desired values on the `make' command line.
$(RECURSIVE_TARGETS):
	@failcom='exit 1'; \
	for f in x $$MAKEFLAGS; do \
	  case $$f in \
	    *=* | --[!k]*);; \
	    *k*) failcom='fail=yes';; \
	  esac; \
	done; \
	dot_seen=no; \
	target=`echo $@ | sed s/-recursive//`; \
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  echo "Making $$target in $$subdir"; \
	  if test "$$subdir" = "."; then \
	    dot_seen=yes; \
	    local_target="$$target-am"; \
	  else \
	    local_target="$$target"; \
	  fi; \
	  ($(am__cd) $$subdir && $(MAKE) $(AM_MAKEFLAGS) $$local_target) \
	  || eval $$failcom; \
	done; \
	if test "$$dot_seen" = "no"; then \
	  $(MAKE) $(AM_MAKEFLAGS) "$$target-am" || exit 1; \
	fi; test -z "$$fail"

$(RECURSIVE_CLEAN_TARGETS):
	@failcom='exit 1'; \
	for f in x $$MAKEFLAGS; do \
	  case $$f in \
	    *=* | --[!k]*);; \
	    *k*) failcom='fail=yes';; \
	  esac; \
	done; \
	dot_seen=no; \
	case "$@" in \
	  distclean-* | maintainer-clean-*) list='$(DIST_SUBDIRS)' ;; \
	  *) list='$(SUBDIRS)' ;; \
	esac; \
	rev=''; for subdir in $$list; do \
	  if test "$$subdir" = "."; then :; else \
	    rev="$$subdir $$rev"; \
	  fi; \
	done; \
	rev="$$rev ."; \
	target=`echo $@ | sed s/-recursive//`; \
	for subdir in $$rev; do \
	  echo "Making $$target in $$subdir"; \
	  if test "$$subdir" = "."; then \
	    local_target="$$target-am"; \
	  else \
	    local_target="$$target"; \
	  fi; \
	  ($(am__cd) $$subdir && $(MAKE) $(AM_MAKEFLAGS) $$local_target) \
	  || eval $$failcom; \
	done && test -z "$$fail"
tags-recursive:
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  test "$$subdir" = . || ($(am__cd) $$subdir && $(MAKE) $(AM_MAKEFLAGS) tags); \
	done
ctags-recursive:
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  test "$$subdir" = . || ($(am__cd) $$subdir && $(MAKE) $(AM_MAKEFLAGS) ctags); \
	done

ID: $(HEADERS) $(SOURCES) $(LISP) $(TAGS_FILES)
	list='$(SOURCES) $(HEADERS) $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	mkid -fID $$unique
tags: TAGS

TAGS: tags-recursive $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	set x; \
	here=`pwd`; \
	if ($(ETAGS) --etags-include --version) >/dev/null 2>&1; then \
	  include_option=--etags-include; \
	  empty_fix=.; \
	else \
	  include_option=--include; \
	  empty_fix=; \
	fi; \
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  if test "$$subdir" = .; then :; else \
	    test ! -f $$subdir/TAGS || \
	      set "$$@" "$$include_option=$$here/$$subdir/TAGS"; \
	  fi; \
	done; \
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	shift; \
	if test -z "$(ETAGS_ARGS)$$*$$unique"; then :; else \
	  test -n "$$unique" || unique=$$empty_fix; \
	  if test $$# -gt 0; then \
	    $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	      "$$@" $$unique; \
	  else \
	    $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	      $$unique; \
	  fi; \
	fi
ctags: CTAGS
CTAGS: ctags-recursive $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	test -z "$(CTAGS_ARGS)$$unique" \
	  || $(CTAGS) $(CTAGSFLAGS) $(AM_CTAGSFLAGS) $(CTAGS_ARGS) \
	     $$unique

GTAGS:
	here=`$(am__cd) $(top_builddir) && pwd` \
	  && $(am__cd) $(top_srcdir) \
	  && gtags -i $(GTAGS_ARGS) "$$here"

distclean-tags:
	-rm -f TAGS ID GTAGS GRTAGS GSYMS GPATH tags

distdir: $(DISTFILES)
	$(am__remove_distdir)
	test -d "$(distdir)" || mkdir "$(distdir)"
	@srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  sed -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   sed '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  if test -f $$file || test -d $$file; then d=.; else d=$(srcdir); fi; \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | sed -e 's,/[^/]*$$,,'`; \
	    if test -d "$(distdir)/$$file"; then \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -fpR $(srcdir)/$$file "$(distdir)$$dir" || exit 1; \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    cp -fpR $$d/$$file "$(distdir)$$dir" || exit 1; \
	  else \
	    test -f "$(distdir)/$$file" \
	    || cp -p $$d/$$file "$(distdir)/$$file" \
	    || exit 1; \
	  fi; \
	done
	@list='$(DIST_SUBDIRS)'; for subdir in $$list; do \
	  if test "$$subdir" = .; then :; else \
	    test -d "$(distdir)/$$subdir" \
	    || $(MKDIR_P) "$(distdir)/$$subdir" \
	    || exit 1; \
	  fi; \
	done
	@list='$(DIST_SUBDIRS)'; for subdir in $$list; do \
	  if test "$$subdir" = .; then :; else \
	    dir1=$$subdir; dir2="$(distdir)/$$subdir"; \
	    $(am__relativize); \
	    new_distdir=$$reldir; \
	    dir1=$$subdir; dir2="$(top_distdir)"; \
	    $(am__relativize); \
	    new_top_distdir=$$reldir; \
	    echo " (cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) top_distdir="$$new_top_distdir" distdir="$$new_distdir" \\"; \
	    echo "     am__remove_distdir=: am__skip_length_check=: am__skip_mode_fix=: distdir)"; \
	    ($(am__cd) $$subdir && \
	      $(MAKE) $(AM_MAKEFLAGS) \
	        top_distdir="$$new_top_distdir" \
	        distdir="$$new_distdir" \
		am__remove_distdir=: \
		am__skip_length_check=: \
		am__skip_mode_fix=: \
	        distdir) \
	      || exit 1; \
	  fi; \
	done
	-test -n "$(am__skip_mode_fix)" \
	|| find "$(distdir)" -type d ! -perm -777 -exec chmod a+rwx {} \; -o \
	  ! -type d ! -perm -444 -links 1 -exec chmod a+r {} \; -o \
	  ! -type d ! -perm -400 -exec chmod a+r {} \; -o \
	  ! -type d ! -perm -444 -exec $(install_sh) -c -m a+r {} {} \; \
	|| chmod -R a+r "$(distdir)"
dist-gzip: distdir
	tardir=$(distdir) && $(am__tar) | GZIP=$(GZIP_ENV) gzip -c >$(distdir).tar.gz
	$(am__remove_distdir)

dist-bzip2: distdir
	tardir=$(distdir) && $(am__tar) | bzip2 -9 -c >$(distdir).tar.bz2
	$(am__remove_distdir)

dist-lzma: distdir
	tardir=$(distdir) && $(am__tar) | lzma -9 -c >$(distdir).tar.lzma
	$(am__remove_distdir)

dist-xz: distdir
	tardir=$(distdir) && $(am__tar) | xz -c >$(distdir).tar.xz
	$(am__remove_distdir)

dist-tarZ: distdir
	tardir=$(distdir) && $(am__tar) | compress -c >$(distdir).tar.Z
	$(am__remove_distdir)

dist-shar: distdir
	shar $(distdir) | GZIP=$(GZIP_ENV) gzip -c >$(distdir).shar.gz
	$(am__remove_distdir)

dist-zip: distdir
	-rm -f $(distdir).zip
	zip -rq $(distdir).zip $(distdir)
	$(am__remove_distdir)

dist dist-all: distdir
	tardir=$(distdir) && $(am__tar) | GZIP=$(GZIP_ENV) gzip -c >$(distdir).tar.gz
	$(am__remove_distdir)

# This target untars the dist file and tries a VPATH configuration.  Then
# it guarantees that the distribution is self-contained by making another
# tarfile.
distcheck: dist
	case '$(DIST_ARCHIVES)' in \
	*.tar.gz*) \
	  GZIP=$(GZIP_ENV) gunzip -c $(distdir).tar.gz | $(am__untar) ;;\
	*.tar.bz2*) \
	  bunzip2 -c $(distdir).tar.bz2 | $(am__untar) ;;\
	*.tar.lzma*) \
	  unlzma -c $(distdir).tar.lzma | $(am__untar) ;;\
	*.tar.xz*) \
	  xz -dc $(distdir).tar.xz | $(am__untar) ;;\
	*.tar.Z*) \
	  uncompress -c $(distdir).tar.Z | $(am__untar) ;;\
	*.shar.gz*) \
	  GZIP=$(GZIP_ENV) gunzip -c $(distdir).shar.gz | unshar ;;\
	*.zip*) \
	  unzip $(distdir).zip ;;\
	esac
	chmod -R a-w $(distdir); chmod a+w $(distdir)
	mkdir $(distdir)/_build
	mkdir $(distdir)/_inst
	chmod a-w $(distdir)
	test -d $(distdir)/_build || exit 0; \
	dc_install_base=`$(am__cd) $(distdir)/_inst && pwd | sed -e 's,^[^:\\/]:[\\/],/,'` \
	  && dc_destdir="$${TMPDIR-/tmp}/am-dc-$$$$/" \
	  && am__cwd=`pwd` \
	  && $(am__cd) $(distdir)/_build \
	  && ../configure --srcdir=.. --prefix="$$dc_install_base" \
	    $(DISTCHECK_CONFIGURE_FLAGS) \
	  && $(MAKE) $(AM_MAKEFLAGS) \
	  && $(MAKE) $(AM_MAKEFLAGS) dvi \
	  && $(MAKE) $(AM_MAKEFLAGS) check \
	  && $(MAKE) $(AM_MAKEFLAGS) install \
	  && $(MAKE) $(AM_MAKEFLAGS) installcheck \
	  && $(MAKE) $(AM_MAKEFLAGS) uninstall \
	  && $(MAKE) $(AM_MAKEFLAGS) distuninstallcheck_dir="$$dc_install_base" \
	        distuninstallcheck \
	  && chmod -R a-w "$$dc_install_base" \
	  && ({ \
	       (cd ../.. && umask 077 && mkdir "$$dc_destdir") \
	       && $(MAKE) $(AM_MAKEFLAGS) DESTDIR="$$dc_destdir" install \
	       && $(MAKE) $(AM_MAKEFLAGS) DESTDIR="$$dc_destdir" uninstall \
	       && $(MAKE) $(AM_MAKEFLAGS) DESTDIR="$$dc_destdir" \
	            distuninstallcheck_dir="$$dc_destdir" distuninstallcheck; \
	      } || { rm -rf "$$dc_destdir"; exit 1; }) \
	  && rm -rf "$$dc_destdir" \
	  && $(MAKE) $(AM_MAKEFLAGS) dist \
	  && rm -rf $(DIST_ARCHIVES) \
	  && $(MAKE) $(AM_MAKEFLAGS) distcleancheck \
	  && cd "$$am__cwd" \
	  || exit 1
	$(am__remove_distdir)
	@(echo "$(distdir) archives ready for distribution: "; \
	  list='$(DIST_ARCHIVES)'; for i in $$list; do echo $$i; done) | \
	  sed -e 1h -e 1s/./=/g -e 1p -e 1x -e '$$p' -e '$$x'
distuninstallcheck:
	@$(am__cd) '$(distuninstallcheck_dir)' \
	&& test `$(distuninstallcheck_listfiles) | wc -l` -le 1 \
	   || { echo "ERROR: files left after uninstall:" ; \
	        if test -n "$(DESTDIR)"; then \
	          echo "  (check DESTDIR support)"; \
	        fi ; \
	        $(distuninstallcheck_listfiles) ; \
	        exit 1; } >&2
distcleancheck: distclean
	@if test '$(srcdir)' = . ; then \
	  echo "ERROR: distcleancheck can only run from a VPATH build" ; \
	  exit 1 ; \
	fi
	@test `$(distcleancheck_listfiles) | wc -l` -eq 0 \
	  || { echo "ERROR: files left in build directory after distclean:" ; \
	       $(distcleancheck_listfiles) ; \
	       exit 1; } >&2
check-am: all-am
check: check-recursive
all-am: Makefile all-local
installdirs: installdirs-recursive
installdirs-am:
install: install-recursive
install-exec: install-exec-recursive
install-data: install-data-recursive
uninstall: uninstall-recursive

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am

installcheck: installcheck-recursive
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) INSTALL_PROGRAM="$(INSTALL_STRIP_PROGRAM)" \
	  install_sh_PROGRAM="$(INSTALL_STRIP_PROGRAM)" INSTALL_STRIP_FLAG=-s \
	  `test -z '$(STRIP)' || \
	    echo "INSTALL_PROGRAM_ENV=STRIPPROG='$(STRIP)'"` install
mostlyclean-generic:

clean-generic:

distclean-generic:
	-test -z "$(CONFIG_CLEAN_FILES)" || rm -f $(CONFIG_CLEAN_FILES)
	-test . = "$(srcdir)" || test -z "$(CONFIG_CLEAN_VPATH_FILES)" || rm -f $(CONFIG_CLEAN_VPATH_FILES)
	-test -z "$(DISTCLEANFILES)" || rm -f $(DISTCLEANFILES)

maintainer-clean-generic:
	@echo "This command is intended for maintainers to use"
	@echo "it deletes files that may require special tools to rebuild."
clean: clean-recursive

clean-am: clean-generic clean-libtool mostlyclean-am

distclean: distclean-recursive
	-rm -f $(am__CONFIG_DISTCLEAN_FILES)
	-rm -f Makefile
distclean-am: clean-am distclean-generic distclean-libtool \
	distclean-tags

dvi: dvi-recursive

dvi-am:

html: html-recursive

html-am:

info: info-recursive

info-am:

install-data-am:

install-dvi: install-dvi-recursive

install-dvi-am:

install-exec-am: install-exec-local

install-html: install-html-recursive

install-html-am:

install-info: install-info-recursive

install-info-am:

install-man:

install-pdf: install-pdf-recursive

install-pdf-am:

install-ps: install-ps-recursive

install-ps-am:

installcheck-am: installcheck-local

maintainer-clean: maintainer-clean-recursive
	-rm -f $(am__CONFIG_DISTCLEAN_FILES)
	-rm -rf $(top_srcdir)/autom4te.cache
	-rm -f Makefile
maintainer-clean-am: distclean-am maintainer-clean-generic

mostlyclean: mostlyclean-recursive

mostlyclean-am: mostlyclean-generic mostlyclean-libtool \
	mostlyclean-local

pdf: pdf-recursive

pdf-am:

ps: ps-recursive

ps-am:

uninstall-am: uninstall-local

.MAKE: $(RECURSIVE_CLEAN_TARGETS) $(RECURSIVE_TARGETS) ctags-recursive \
	install-am install-strip tags-recursive

.PHONY: $(RECURSIVE_CLEAN_TARGETS) $(RECURSIVE_TARGETS) CTAGS GTAGS \
	all all-am all-local am--refresh check check-am clean \
	clean-generic clean-libtool ctags ctags-recursive dist \
	dist-all dist-bzip2 dist-gzip dist-lzma dist-shar dist-tarZ \
	dist-xz dist-zip distcheck distclean distclean-generic \
	distclean-libtool distclean-tags distcleancheck distdir \
	distuninstallcheck dvi dvi-am html html-am info info-am \
	install install-am install-data install-data-am install-dvi \
	install-dvi-am install-exec install-exec-am install-exec-local \
	install-html install-html-am install-info install-info-am \
	install-man install-pdf install-pdf-am install-ps \
	install-ps-am install-strip installcheck installcheck-am \
	installcheck-local installdirs installdirs-am maintainer-clean \
	maintainer-clean-generic mostlyclean mostlyclean-generic \
	mostlyclean-libtool mostlyclean-local pdf pdf-am ps ps-am tags \
	tags-recursive uninstall uninstall-am uninstall-local


# List all build rules defined by HDF5 Makefiles as "PHONY" targets here.
# This tells the Makefiles that these targets are not files to be built but
# commands that should be executed even if a file with the same name already
# exists.
.PHONY: build-check-clean build-check-p build-check-s build-lib build-progs \
        build-tests check-clean check-install check-p check-s check-vfd \
        install-doc lib progs tests uninstall-doc _exec_check-s _test help

help:
	@$(top_srcdir)/bin/makehelp

# Define rules for lib, progs, check, and tests.
# These simply involve recursing into subdirectories.
test _test: check

lib progs check-p check-s:
	@set -x; for d in $(SUBDIRS); do                                \
	  if test $$d != .; then                                        \
	   (cd $$d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;            \
	  fi;                                                           \
	done

# Make all, tests, and (un)install also recurse into perform directory
all-local:
	@cd perform && $(MAKE) $(AM_MAKEFLAGS) all

tests:
	@set -x; for d in $(SUBDIRS) perform; do                        \
	  if test $$d != .; then                                        \
	   (cd $$d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;            \
	  fi;                                                           \
	done

install-exec-local:
	@cd perform && $(MAKE) $(AM_MAKEFLAGS) install
uninstall-local:
	@cd perform && $(MAKE) $(AM_MAKEFLAGS) uninstall

# Check-clean also recurses into examples directory
check-clean:
	@set -x; for d in $(SUBDIRS) examples perform; do               \
	  if test $$d != .; then                                        \
	   (cd $$d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;            \
	  fi;                                                           \
	done
	$(RM) -rf prefix-tmp destdir-tmp

# Some C++ compilers/linkers will create a directory named ii_files in
# the root directory, which should be cleaned.
mostlyclean-local:
	if test -d ii_files; then                                       \
	  $(RM) -rf ii_files;                                           \
	fi

# 'make install-all' also installs examples
install-all:
	@$(MAKE) $(AM_MAKEFLAGS) install
	@$(MAKE) $(AM_MAKEFLAGS) install-examples
uninstall-all:
	@$(MAKE) $(AM_MAKEFLAGS) uninstall
	@$(MAKE) $(AM_MAKEFLAGS) uninstall-examples

# Install examples in this directory and recursively
install-examples uninstall-examples:
	@set -x; for d in examples $(HDF5_INTERFACES); do               \
	    (cd $$d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;           \
	done

# Placeholder to remind users that docs are now in a separate repository.
install-doc:
	@echo "docs no longer live in this tree.  Use install-examples to install examples."

uninstall-doc:
	@echo "docs no longer live in this tree.  Use install-examples to install examples."

# `make check-install' or `make installcheck' checks that examples can 
#  be successfully built
installcheck-local:
	@(cd examples && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;

# check-install is just a synonym for installcheck
check-install: installcheck

# check-all-install tests all installation methods.
# Install via different mechanism and then compare against the default.
# Fine if only libXXX.a files are different since they may have been ranlib'ed.
check-all-install:
	@echo Installing to default location
	$(MAKE) install
	@echo Installing to different prefix location
	$(MAKE) prefix=${ROOT}/prefix-tmp install
	@echo Compare against the default installation.
	@echo Fine if only libXXX.a files are different.
	-diff -r prefix-tmp ${prefix}
	@echo Installing to different $$\DESTDIR location
	env DESTDIR=${ROOT}/destdir-tmp $(MAKE) install
	@echo Compare against the default installation.
	@echo Fine if only libXXX.a files are different.
	-diff -r destdir-tmp${prefix} ${prefix}

# Only source files in the src directory include tracing information,
# so 'make trace' only needs to recurse into that directory.
trace:
	@(cd src && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;

# Run tests with different Virtual File Drivers.
# Currently, only invoke check-vfd in the test directory.
check-vfd:
	@set -x; for d in src test; do                                 \
	  if test $$d != .; then                                       \
	    (cd $$d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;          \
	  fi;                                                          \
	done

# Automake wants to make config.status depend on configure.  This
# makes sense, but config.status can't always be regenerated
# properly, which can cause builds to fail.
# This is a problem for our Daily Tests, which need to be able to
# 'make distclean' reliably before running configure.
# The simple solution is to override the dependency Automake supplies
# for config.status so that it will never be regenerated.
$(top_builddir)/config.status:

# Don't include conclude.am in root Makefile; tests target needs to
# recurse into perform directory as well as reguar subdirs.

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
