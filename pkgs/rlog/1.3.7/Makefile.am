SUBDIRS = rlog docs config

EXTRA_DIST = AUTHORS COPYING INSTALL README rlog.lsm Doxyfile win32

# not a GNU package. You can remove this line, if
# have all needed files, that a GNU package needs
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = librlog.pc

documentation:
	cd docs && make

