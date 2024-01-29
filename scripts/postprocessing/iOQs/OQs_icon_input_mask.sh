#! /bin/sh

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

FX=$POOL/$GRID/$GRID${LEV}_fx.nc
CDO=${CDO:-cdo -s}


make -j 4 -f - << EOF

$(basename $0): mask_3d.nc mask_antarc.nc mask_arc.nc mask_glo.nc mask_atl.nc mask_atl_arc.nc mask_indopacific.nc FORCE

mask_3d.nc: ${FX}
	${CDO} div -selvar,wet_c ${FX} -selvar,wet_c ${FX} mask_3d.nc

mask_antarc.nc: ${FX}
	${CDO} -mul -selvar,wet_c ${FX} -gtc,1 -setvrange,6,6 -selvar,regio_c ${FX}  mask_antarc.nc

mask_arc.nc: ${FX}
	${CDO} -mul -selvar,wet_c ${FX} -gtc,1 -setvrange,1,4 -selvar,regio_c ${FX} mask_arc.nc

mask_glo.nc: ${FX}
	${CDO} -mul -selvar,wet_c ${FX} -gtc,1 -setvrange,1,9 -selvar,regio_c ${FX} mask_glo.nc
	${CDO} -sellevidx,1  mask_glo.nc  mask_glo_surf.nc
	${CDO} zonmean -remapnn,r360x180 -setctomiss,0 -div mask_glo.nc mask_glo.nc mask_glo_zon.nc

mask_atl.nc: ${FX}
	${CDO} -mul -selvar,wet_c ${FX} -gtc,1 -setvrange,4,4 -selvar,regio_c ${FX} mask_atl.nc

mask_atl_arc.nc: ${FX}
	${CDO} -mul -selvar,wet_c ${FX} -gtc,1 -setvrange,1,5 -selvar,regio_c ${FX} mask_atl_arc.nc
	${CDO} zonmean -remapnn,r360x180 -setctomiss,0 -div mask_atl_arc.nc mask_atl_arc.nc mask_atl_arc_zon.nc

mask_indopacific.nc: ${FX}
	${CDO} -mul -selvar,wet_c ${FX} -gtc,1 -setvrange,7,8 -selvar,regio_c ${FX} mask_indopacific.nc
	${CDO} zonmean -remapnn,r360x180 -setctomiss,0  -div mask_indopacific.nc mask_indopacific.nc mask_indopacific_zon.nc

FORCE:

EOF
