#!/bin/ksh

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

FILES_DYNAMICS='model_equations_ICONAM.tex terrain_following_ICONAM.tex vortex_bracket_ICONAM.tex specials_triangle_ICONAM.tex fig_primal_dual.pdf fig_main_boxes.pdf fig_edge_volume.pdf fig_vortex_boxes.pdf'
for FILE in $FILES_DYNAMICS; do
ln -s ../dynamics/$FILE $FILE
done
make icon_atm_echam_phy_scidoc
