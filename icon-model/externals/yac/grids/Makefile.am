# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

icon_grids = \
    icon_grid_0030_R02B03_G.nc \
    icon_grid_0036_R02B04_O.nc \
    icon_grid_0043_R02B04_G.nc
TEST_DATA_URL = http://icon-downloads.mpimet.mpg.de/grids/public/mpim

if HAVE_NETCDF

if ENABLE_EXAMPLES
all-local: $(icon_grids)
else
check-local: $(icon_grids)
endif

else

if ENABLE_PYTHON_BINDINGS
check-local: $(icon_grids)
endif

endif

$(icon_grids):
	@url="$(TEST_DATA_URL)/`echo $@ | cut -c 11-14`/$@"; \
	cmd="$(CURL) -sL '$$url' -o $@"; \
	if test -n '$(AM_V_GEN)'; then \
	  echo "  FETCH   " $@; \
	else \
	  echo "$$cmd"; \
	fi; \
	eval "$$cmd"

distclean-local:
	-rm -f $(icon_grids)

EXTRA_DIST = \
    GR30_lsm.nc
