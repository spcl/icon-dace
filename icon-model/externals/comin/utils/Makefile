FILES_IN := \
  ../src/comin_descrdata_types.F90
FILES_OUT := \
  ../src/comin_descrdata_query_domain.inc \
  ../src/comin_descrdata_query_global.inc \
  ../include/comin_header_c_ext_descrdata_query_domain.h \
  ../include/comin_header_c_ext_descrdata_query_global.h \
  ../replay_tool/src/comin_descrdata_save.F90 \
  ../replay_tool/src/comin_descrdata_load.F90

all: ${FILES_OUT}

../src/%.F90: %.F90.py ${FILES_IN}
	python3 $< > $@

../src/%.inc: %.inc.py ${FILES_IN}
	python3 $< > $@

../replay_tool/src/%.F90: %.F90.py ${FILES_IN}
	python3 $< > $@

../include/%.h: %.h.py ${FILES_IN}
	python3 $< > $@

clean:
	rm -f ${FILES_OUT}
