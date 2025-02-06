# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

SOURCES = $(filter-out variables.yaml,$(wildcard *.yaml))
TARGETS = $(SOURCES:.yaml=.html)

all: $(TARGETS)

clean:
	$(RM) $(TARGETS)

%.html: %.yaml
	icon_variables render $< $@
