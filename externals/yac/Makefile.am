SUBDIRS = src

if !EXTERNAL_MTIME
SUBDIRS += mtime
endif

if !EXTERNAL_LAPACK
SUBDIRS += clapack
endif

if ENABLE_PYTHON_BINDINGS
SUBDIRS += python
endif

SUBDIRS += contrib

if !ENABLE_LIB_ONLY
SUBDIRS += examples
endif

SUBDIRS += tests

EXTRA_DIST = config/checksrc config/checkdoc
