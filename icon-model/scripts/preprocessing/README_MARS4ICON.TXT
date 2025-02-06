----------------------------------------------------------------------------
Notes on MARS4ICON usage
----------------------------------------------------------------------------
Initial version: F.Prill, DWD (2011-12-16)
Shortened text: F. Prill, DWD (2014-05-28)

Document version:
  $LastChangedDate$
  $Rev: 7638 $ -- $Author$

----------------------------------------------------------------------------



1. GRIB API
-----------

Necessary:  grib_api Version 1.9.9
Check, eg., with "grib_info"; cf. settings in ".profile".

Important: Use grib_api *release* version
           This is necessary for setting "CI" missing values to 0 !
             (. ${GRIB_API}/bin/grib_api-setenv -r -V1.9.9)

Table of GRIB keys:
  http://www.ecmwf.int/publications/manuals/d/gribapi/param
see also
  firefox --no-remote /usr/local/pkg/grib_api/1.8.1/html/GRIB_PARAMETER_html/index.htm )

The IFS2ICON scripts handles GRIB1 input fields by their table codes
and GRIB2 input fields by their short names.

Get GRIB info of retrieved file:
  grib_ls ifs_oper_T255_2011010100.grb

Get GRIB version (GRIB1/GRIB2) of a given file:
  grib_get  -w count=1 -p editionNumber ifs_oper_T1279_2011080100.grb



2. MARS Retrieval
-----------------

Script for retrieval of IFS data (via MARS, ECMWF): See ICON svn repository,

~/trunk/icon-dev/scripts/preprocessing/mars4icon_smi

Run as:
  mars4icon_smi -d 2011010100 -s 0 -r 255 -l 1/to/91 2>&1 | tee log.mars4icon.txt

Note: Requires valid ECMWF certificate (command "eccert" on DWD cluster)

Note: At the end of this script, volumetric soil water layer fields (SWVL1-4) are
      renamed into "soil moisture index" (SMIL1); coded as GRIB indices 80-83. This
      is essential because values are converted without altering the grib index,
      according to the procedure described in

      http://www.ecmwf.int/products/changes/soil_hydrology_cy32r3/#compute_smi


