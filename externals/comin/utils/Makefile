FILES_IN  := $(shell ls *.F90.py *.h.py)
FILES_OUT := $(FILES_IN:.py=)

all: ${FILES_OUT}

%.F90: %.F90.py
	python3 $< > ../src/$@

%.h: %.h.py
	python3 $< > ../include/$@
