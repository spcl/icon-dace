#!/bin/ksh
#   
# create_hddata_tarfile.com - Create HD data tar file for data publication linked from HD code repo.
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
# ****** Create data tar file that will be included in code publication  **************
# 
#        It should comprise HD parameter files for 0.5°, 5 Min. and Euro5Min
#
#   tar cf archive.tar foo bar --transform 's,^,bazdir/,'
#       --> files are archived as bazdir/foo and  bazdir/bar
#   you can run tar with multiple --transform options.
#
DDIR=/work/gg0302/g260122/HD/input
cd $DDIR
#
VS5=vs5_1         # e.g. vs5_0
VS1=vs1_11        # e.g. vs1_10_ext
#
# 0.5 degree
tar cvhf hd${VS5}_input_data.tar 05deg/hdpara_${VS1}.nc 05deg/masks_05.nc 05deg/hdstart_05.nc 
tar uvhf hd${VS5}_input_data.tar 5min/hdpara_${VS5}.nc 5min/masks_5min.nc 5min/hdstart_5min.nc
tar uvhf hd${VS5}_input_data.tar euro5min/hdpara_${VS5}_euro5min.nc euro5min/masks_euro5min.nc euro5min/hdstart_euro5min.nc
#
cat > readme_hd_data_${VS5}.txt << EOF
HD model parameter files

This dataset comprises global parameter data that are necessary to run the Hydrological Discharge (HD) model, which has been published on Zenodo. The HD model calculates the lateral transport of water over the land surface to simulate discharge into the oceans. The HD model parameter dataset comprises parameter fields at 0.5° global resolution and at 5 Min. resolution (global, Europe). Details for both resolutions are provided below.  

Authors: Stefan Hagemann, Tobias Stacke   
Copyright 2021: Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
Licensed under the Creative Commons Attribution-ShareAlike 4.0 International License (CC BY-SA 4.0; https://creativecommons.org/licenses/)

#
# Data description

HD model parameter file at 5 Min resolution: hdpara_${VS5}.nc 

River directions and digital elevation data were provided by Bernhard Lehner (pers. comm., 2014) and were derived from the HydroSHEDS (Lehner et al., 2006) database and from the Hydro1K dataset for areas north of 60°N (https://lta.cr.usgs.gov/HYDRO1K). 
For a number of rivers (most of them north of 60°N), flow directions and model orography were manually corrected based on available GIS data, such as from DIVA (https://www.diva-gis.org/gdata), CCM River and Catchment Database
(Vogt et al. 2007), SMHI (Swedish Meteorological and Hydrological Institute), NVE (Norges vassdrags- og energidirektorats), SYKE (Finnish Environment Institute).

This corrected dataset is referred to as HD${VS5} in the following.
The HD model parameters for overland flow, base flow and river flow are generated as described in Hagemann and Dümenil (1998) and Hagemann et al. (2020). However, different to the HD model vs. 4 described in Hagemann et al. (2020), Vs5 utilizes inland water fractions from the ESA CCI Water Bodies Map v4.0 (Lamarche et al. 2017) and wetland fractions from the Global Lakes and Wetlands Database (Lehner and Döll 2004) instead of the previously used lake and wetlands fractions.

The HD parameter dataset contains 14 variables which are shortly described in the following table. 

    FLAG   | Land sea mask | - 
    FDIR   | Flow direction | - | defined as written below
    ALF_K  | HD model parameter Overland flow k | d-1 
    ALF_N  | HD model parameter Overland flow n | - 
    ARF_K  | HD model parameter Riverflow k | d-1 
    ARF_N  | HD model parameter Riverflow n | - 
    AGF_K  | HD model parameter Baseflow flow k | d-1 
    AREA   | Grid cell area | m-2 | based on own computation
    FILNEW | River flow target indices for longitudes | - 
    FIBNEW | River flow target indices for latitudes | - 
    DISTANCE    | Distance between gridboxes in flow direction | m 
    RIVERLENGTH | Distance between gridbox and the river mouth (or final sink) | km 
    CAT_AREA    | Upstream catchment area of gridbox | km²
    CAT_ID      | Catchment ID of gridbox | -

    Flow directions in variable FDIR are defined as on the Num Pad of a PC keyboard:
    
                        7  8  9
                         \ | /
                          \|/
                        4--5--6
                          /|\
                         / | \
                        1  2  3
    
    Special directions:  5 = Sink point, i.e. no outflow
                         0 = River mouth point in the ocean      
                        -1 = Ocean point, but no river mouth



Forcing data masks file: masks_5min.nc

In the offline HD model version, this file is usually only used to obtain the grid information of the forcing data, i.e. of surface runoff and drainage (subsurface runoff). However, it contains four variables that are read in by the model, and that are actually used un coupled applications within the MPI-ESM. Even though these variables are not used in the HD model offline version, it was decided to keep them in order to allow future developments regarding the usage of these data and to keep some consistency with the HD model code implemented in MPI-ESM.
   ALAKE  | Lake fraction within a grid box     | Lamarche et al. 2017
   GLAC   | Glacier fraction within a grid box  | Hagemann 2002
   SLF    | Land fraction within a grid box     | Lamarche et al. 2017
   SLM    | Land Sea Mask                       | Lamarche et al. 2017
For simplicity, the data provided at the HD model resolution. Hence, these masks can be used when the forcing data are interpolated to the HD model resolution before they are read during the model run. 

This tar archive also include a subset of this global dataset for the European domain, hdpara_${VS5}_euro5min.nc.

In addition, a global 0.5° HD parameter file (hdpara_${VS1}.nc) is included which is an update of the version 1.10. Compared to vs. 1.10, ${VS1} now also (such as Vs 5.) utilizes inland water fractions from the ESA CCI Water Bodies Map v4.0 (Lamarche et al. 2017) and wetland fractions from the Global Lakes and Wetlands Database (Lehner and Döll 2004) instead of the previously used lake and wetlands fractions. In addition several flow directions and catchment borders were updated .

Note that vs 1.10 was consistent to the parameter files used in previous offline and coupled applications of the HD model at 0.5° resolution (see, e.g. studies cited in Sect. 2.1 of Hagemann et al., 2020). Compared to previous versions, only some flow directions (mainly over Europe) were updated for version 1.10. Except for DISTANCE and RIVERLENGTH, it comprises the same variables as for the 5 Min. version, but flow directions and parameters are generated as described in Hagemann and Dümenil (1998) and Hagemann and Dümenil Gates (2001). Here, the 0.5 degree mask file mask_05.nc comprises those masks that were utilized in the HD parameter generation. Only the land fraction is taken from Hagemann (2002) where the HD land sea mask indicates land.


References
Hagemann, S., L. Dümenil (1998) A parameterization of the lateral waterflow for the global scale. Clim. Dyn. 14 (1), 17-31

Hagemann, S., L. Dümenil Gates (2001) Validation of the hydrological cycle of ECMWF and NCEP reanalyses using the MPI hydrological discharge model, J. Geophys. Res. 106, 1503-1510

Hagemann, S., 2002: An improved land surface parameter dataset for global and regional climate models, MPI Report No. 336, Max Planck Institute for Meteorology, Hamburg, Germany 

Hagemann, S., T. Stacke and H. Ho-Hagemann (2020) High resolution discharge simulations over Europe and the Baltic Sea catchment. Front. Earth Sci., 8:12. doi: 10.3389/feart.2020.00012.

Lamarche, C., Santoro, M., Bontemps, S., d’Andrimont, R., Radoux, J., Giustarini, L., Brockmann, C., Wevers, J., Defourny, P. and Arino, O. (2017) Compilation and validation of SAR and optical data products for a complete and global map of inland/ocean water tailored to the climate modeling community. Remote Sensing, 9(1), p.36.

Lehner, B., P. Döll (2004) Development and validation of a global database of lakes, reservoirs and wetlands.
J. Hydrol., 296: 1-22, doi:10.1016/j.jhydrol.2004.03.028.

Vogt, J.V. et al. (2007): A pan-European River and Catchment Database. European Commission - JRC, Luxembourg, (EUR 22920 EN) 120 pp. 

EOF
#
#
cat > history_data.md << EOF2
#  History of major HD parameter file changes

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## HD Model Version 5.1

### Fixed
- Removal of two bugs in the parameter file for Vs. 5.0
  - Accidentally, there was no wetland impact on riverflow.
  - For the wetland fraction at 1 km derived from GLWD wetlands, GLWD type 10 (50-100%) was erroneously associated with a fraction of 0.875 instead of 0.75, and GLWD type 12 was erroneously omitted. The latter is now associated with 0.125. 

### Changed
- Corrected usage of GLWD wetlands
  - Removal of unspecified wetland areas (GLWD types 10-12) for the impact on riverflow that mainly only occur over North America. Now, GLWD wetland types 4-9 affect river flow.
  - Removal of GLWD rivers from ESA waterbodies. For the calculation of lake fraction from ESA water bodies, the fraction of GLWD rivers (type 3) was substracted from the ESA waterbody fraction.
- Scaling factors on 5 Min. model parameters for overland flow (2) and baseflow OLF (4) were implemented into the HD parameter file and taken out of the run script. 

## Reference: HD model Parameter Version: 5.0
EOF2
#
tar ufh hd${VS5}_input_data.tar readme_hd_data_${VS5}.txt LICENSE_HD_parameter.txt history_data.md
gzip hd${VS5}_input_data.tar
echo "hd${VS5}_input_data.tar.gz was generated in " $DDIR
#
