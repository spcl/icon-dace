; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;----------------------------------------------------------------------
; mtgrm_ascii.ncl
;----------------------------------------------------------------------
; Contains
;   - Transformation of a NetCDF meteogram file to
;     rudimentary ASCII output.
; Note: This is an EXAMPLE SCRIPT, much too slow to cope with
;       LARGE DATA SETS! For real applications specifiy, for example
;       the time step instead of printing all of them!
; Note: We do not distinguish between different types of vertical
;       axis/level indices!!
;
; Invocation:
;   ncl -n mtgrm_ascii.ncl  DataFileName='"MTGRM_PE000_patch001.nc"'
;
;----------------------------------------------------------------------
; 09/2011 : F. Prill, DWD
;----------------------------------------------------------------------

; open netCDF file
cdf_file = addfile(DataFileName,"r") 
nstations = dimsizes(cdf_file->station_lon)

colwidth = 15 ; table column width

; loop over meteogram stations:
do istation = 0,(nstations-1)
  
  print(str_get_nl())
  title_str = "STATION NAME: " + charactertostring(cdf_file->station_name(istation,:)) \
  +           "   (idx: "  + cdf_file->station_idx(istation) \
  +           ", blk: " + cdf_file->station_blk(istation) + ")"
  border_str = ""
  do i=1,strlen(title_str)
    border_str = border_str + "="
  end do
  print(tostring(border_str))
  print(tostring(title_str))
  print(tostring(border_str))
  print(str_get_nl())

  print("Initial date : " + charactertostring(cdf_file->date(1,:)))
  print("   HSURF     : " + cdf_file->station_hsurf(istation))
  print("   FR_LAND   : " + cdf_file->station_frland(istation))
  print("   LAT       : " + cdf_file->station_lat(istation))
  print("   LON       : " + cdf_file->station_lon(istation))
  print("   FC        : " + cdf_file->station_fc(istation))
  print("   SOIL TYPE : " + cdf_file->station_soiltype(istation))
  
  ; build title string with variable names and units
  ; atmospheric (3D) variables:
  nvars = dimsizes(cdf_file->var_name)
  column_str = "            LVL"
  unit_str   = "               "
  do ivar=0,(nvars(0)-1)
    scolumn = charactertostring(cdf_file->var_name(ivar,:))
    do i=1,(colwidth-strlen(scolumn))
      column_str = column_str + " "
    end do
    column_str = column_str + scolumn

    scolumn = "[" + charactertostring(cdf_file->var_unit(ivar,:)) + "]"
    do i=1,(colwidth-strlen(scolumn))
      unit_str = unit_str + " "
    end do
    unit_str = unit_str + scolumn
  end do

  ; surface (2D) variables:
  nsfcvars = dimsizes(cdf_file->sfcvar_name)
  sfccolumn_str = ""
  sfcunit_str   = ""
  do ivar=0,(nsfcvars(0)-1)
    scolumn = charactertostring(cdf_file->sfcvar_name(ivar,:))
    do i=1,(colwidth-strlen(scolumn))
      sfccolumn_str = sfccolumn_str + " "
    end do
    sfccolumn_str = sfccolumn_str + scolumn

    scolumn = "[" + charactertostring(cdf_file->sfcvar_unit(ivar,:)) + "]"
    do i=1,(colwidth-strlen(scolumn))
      sfcunit_str = sfcunit_str + " "
    end do
    sfcunit_str = sfcunit_str + scolumn
  end do

  ntimesteps = dimsizes(cdf_file->time_step)
  max_nlevs  = dimsizes(cdf_file->var_levels)
  do itime=0,(ntimesteps-1)
    print(str_get_nl())
    title_str = "DATE " + charactertostring(cdf_file->date(itime,:))  \
    +           " (time step: " + cdf_file->time_step(itime) + ")"  \
    +           str_get_nl()
    border_str = ""
    do i=1,(strlen(title_str)-1)
      border_str = border_str + "-"
    end do
    print(tostring(title_str+border_str))

    ; ------------------------------------------------
    ; Atmospheric (3D) variables:    

    ; build a matrix containing values of sampled variables
    print(str_get_nl() + "ATMOSPHERIC VARIABLES:" + str_get_nl())
    print(tostring(column_str))
    print(tostring(unit_str))

    values = new( (/ max_nlevs(1), nvars(0)+1 /), double)
    do ilev=0,(max_nlevs(1)-1)
      values(ilev, 0) = ilev
    end do
    do ivar=0,(nvars(0)-1)
      do ilev=0,(cdf_file->var_nlevs(ivar)-1)
        values(ilev, ivar+1) = cdf_file->values(itime, ilev, ivar, istation)
      end do ; ilev
    end do ; ivar
    ; write matrix
    opt        = True
    opt@tspace = 40
    write_matrix (values, "14f"+tostring(colwidth)+".3 ", opt)

    ; ------------------------------------------------
    ; Surface (2D) variables:    

    ; build a matrix containing values of sampled variables
    print(str_get_nl() + "SURFACE VARIABLES:" + str_get_nl())
    print(tostring(sfccolumn_str))
    print(tostring(sfcunit_str))

    sfcvalues = new( (/ 1, nsfcvars(0) /), double)
    do ivar=0,(nsfcvars(0)-1)
      sfcvalues(0,ivar) = cdf_file->sfcvalues(itime, ivar, istation)
    end do ; ivar
    ; write matrix
    opt        = True
    opt@tspace = 40
    write_matrix (sfcvalues, "14f"+tostring(colwidth)+".3 ", opt)

  end do ; itime

end do ; istation
