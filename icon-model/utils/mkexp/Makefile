program := $(shell python -Ssc 'from package_info import name; print(name)')
version := $(shell python -Ssc 'from package_info import version; print(version)')
package = $(program)-$(version)

prefix = /usr/local

PYTHON = python
PIP = pip

all:
	@set -x;\
	if ./set_doc_release doc/$(program).fodt > doc/$(program).fodt~;\
	then mv doc/$(program).fodt~ doc/$(program).fodt;\
	else $(RM) doc/$(program).fodt~; fi

package:
	@echo $(package)

clean:
	$(RM) -r build dist MANIFEST *.pyc __pycache__ mkexp.egg-info

distclean: clean
	$(RM) doc/*.pdf

doc: doc/mkexp.pdf

check:
	PYTHON=$(PYTHON) $(PYTHON) test.py -v

# dist: doc
# 	$(PYTHON) setup.py sdist
dist:
	mkdir -p dist
	git archive --prefix $(package)/ -o dist/$(package).tar.gz HEAD
	ln -f doc/$(program).pdf dist/$(package).pdf

install: doc
	$(PIP) install --prefix=$(prefix) .

.PHONY: all clean distclean doc check dist install

%.pdf: %.fodt
	soffice --headless --convert-to pdf --outdir $(@D) $<
