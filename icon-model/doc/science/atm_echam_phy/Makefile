# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

SHELL = /bin/bash

I_SOURCES=icon_atm_echam_phy_scidoc.tex
FILES_DYNAMICS=../dynamics/model_equations_ICONAM.tex \
               ../dynamics/terrain_following_ICONAM.tex \
               ../dynamics/vortex_bracket_ICONAM.tex \
               ../dynamics/specials_triangle_ICONAM.tex \
               ../dynamics/fig_primal_dual.pdf \
               ../dynamics/fig_main_boxes.pdf \
               ../dynamics/fig_edge_volume.pdf \
               ../dynamics/fig_vortex_boxes.pdf

icon_atm_echam_phy_scidoc: $(FILES_DYNAMICS) $(I_SOURCES)
	pdflatex icon_atm_echam_phy_scidoc
	makeindex icon_atm_echam_phy_scidoc
	bibtex icon_atm_echam_phy_scidoc
	pdflatex icon_atm_echam_phy_scidoc
	pdflatex icon_atm_echam_phy_scidoc
