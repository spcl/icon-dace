# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
@DX_RULES@

.NOTPARALLEL:

.PHONY: man install-man do-install-man \
        html install-html do-install-html \
        pdf install-pdf do-install-pdf \
        ps install-ps do-install-ps

if DX_FLAG_man
STAMP_MAN = $(DX_DOCDIR)/$(PACKAGE).tag
STAMP_INSTALL_MAN = do-install-man
else !DX_FLAG_man
STAMP_MAN =
STAMP_INSTALL_MAN =
endif !DX_FLAG_man

if DX_FLAG_html
STAMP_HTML = $(DX_DOCDIR)/$(PACKAGE).tag
STAMP_INSTALL_HTML = do-install-html
else !DX_FLAG_html
STAMP_HTML =
STAMP_INSTALL_HTML =
endif !DX_FLAG_html

if DX_FLAG_pdf
STAMP_PDF = $(DX_DOCDIR)/$(PACKAGE).pdf
STAMP_INSTALL_PDF = do-install-pdf
else !DX_FLAG_pdf
STAMP_PDF =
STAMP_INSTALL_PDF =
endif !DX_FLAG_pdf

if DX_FLAG_ps
STAMP_PS = $(DX_DOCDIR)/$(PACKAGE).ps
STAMP_INSTALL_PS = do-install-ps
else !DX_FLAG_ps
STAMP_PS =
STAMP_INSTALL_PS =
endif !DX_FLAG_ps

man: $(STAMP_MAN)
install-man: $(STAMP_INSTALL_MAN)
do-install-man: man
	@$(NORMAL_INSTALL)
	@echo " $(MKDIR_P) '$(DESTDIR)$(mandir)/man3'"; \
	$(MKDIR_P) "$(DESTDIR)$(mandir)/man3" || exit 1; \
	echo " $(INSTALL_DATA) '$(DX_DOCDIR)/man/man3'/* '$(DESTDIR)$(mandir)/man3'"; \
	$(INSTALL_DATA) "$(DX_DOCDIR)/man/man3"/* "$(DESTDIR)$(mandir)/man3"

html: $(STAMP_HTML)
install-html: $(STAMP_INSTALL_HTML)
do-install-html: html
	@$(NORMAL_INSTALL)
	@echo " $(MKDIR_P) '$(DESTDIR)$(htmldir)'"; \
	$(MKDIR_P) "$(DESTDIR)$(htmldir)" || exit 1; \
	echo " cp -R '$(DX_DOCDIR)/html'/* '$(DESTDIR)$(htmldir)'"; \
	cp -R "$(DX_DOCDIR)/html"/* "$(DESTDIR)$(htmldir)"

pdf: $(STAMP_PDF)
install-pdf: $(STAMP_INSTALL_PDF)
do-install-pdf: pdf
	@$(NORMAL_INSTALL)
	@echo " $(MKDIR_P) '$(DESTDIR)$(pdfdir)'"; \
	$(MKDIR_P) "$(DESTDIR)$(pdfdir)" || exit 1; \
	echo " $(INSTALL_DATA) '$(DX_DOCDIR)/$(PACKAGE).pdf' '$(DESTDIR)$(pdfdir)'"; \
	$(INSTALL_DATA) "$(DX_DOCDIR)/$(PACKAGE).pdf" "$(DESTDIR)$(pdfdir)"

ps: $(STAMP_PS)
install-ps: $(STAMP_INSTALL_PS)
do-install-ps: ps
	@$(NORMAL_INSTALL)
	@echo " $(MKDIR_P) '$(DESTDIR)$(psdir)'"; \
	$(MKDIR_P) "$(DESTDIR)$(psdir)" || exit 1; \
	echo " $(INSTALL_DATA) '$(DX_DOCDIR)/$(PACKAGE).ps' '$(DESTDIR)$(psdir)'"; \
	$(INSTALL_DATA) "$(DX_DOCDIR)/$(PACKAGE).ps" "$(DESTDIR)$(psdir)"

doc: man html pdf ps
install-doc: install-man install-html install-pdf install-ps

EXTRA_DIST = libmtime.css libmtime.dox libmtime.tex

MOSTLYCLEANFILES = $(DX_CLEANFILES) -r $(DX_DOCDIR)

