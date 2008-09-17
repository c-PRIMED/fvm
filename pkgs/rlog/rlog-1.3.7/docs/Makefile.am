
SOURCES = index.dox channel.dox component.dox

docdir = $(datadir)/doc
pkgdocdir = $(docdir)/$(PACKAGE)
#pkghtmldir = $(pkgdocdir)/html


# have both html and html/index.html to force build of index.html if necessary
EXTRA_DIST = latex/refman.pdf html html/index.html $(SOURCES)

pkgdoc_DATA = \
    latex/refman.pdf

#pkghtml_DATA = \
#    html/*

if BUILD_DOCS
latex/refman.pdf: html/index.html
	-cd latex && make

html/index.html: ../Doxyfile $(SOURCES)
	-rm -rf html
	-cd .. && @DOXYGEN@ Doxyfile

clean-local:
	-rm -rf html
	-rm -rf latex
endif

all-local: html/index.html

install-data-local: all-local
	-echo installing to ${DESTDIR}${pkgdocdir}/html
	/bin/sh ../mkinstalldirs ${DESTDIR}${pkgdocdir}/html
	install html/* ${DESTDIR}${pkgdocdir}/html


