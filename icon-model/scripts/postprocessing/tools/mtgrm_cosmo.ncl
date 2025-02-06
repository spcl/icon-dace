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
; mtgrm_cosmo.ncl
;----------------------------------------------------------------------
; Contains
;   - Transformation of a NetCDF meteogram file to
;     ASCII output similar to COSMO long meteogram output.
;
; Invocation:
;   ncl -n mtgrm_cosmo.ncl  DataFileName='"METEOGRAM_patch001.nc"' itime=0
;
;----------------------------------------------------------------------
; 11/2011 : F. Prill, DWD
;----------------------------------------------------------------------

; Utility functions----------------------------------------------------
procedure transform_unit(varname, src_unit, src_val, dst_unit, dst_val)
begin
  dst_unit = src_unit
  dst_val  = 1.*src_val
  if (((varname .eq. "QV") .or. (varname .eq. "QC") .or. (varname .eq. "QI") .or. (varname .eq. "QR") .or. (varname .eq. "QS")) .and. (src_unit .eq. "kg kg-1")) then
    dst_unit = "mg kg-1"
    dst_val  = 1.e3 * src_val
    return
  end if
  if ((str_lower(varname) .eq. "qv_s") .and. (src_unit .eq. "kg/kg")) then
    dst_unit = "g/kg"
    dst_val  = 1.e3 * src_val
    return
  end if
  if (src_unit .eq. "Pa") then
    dst_unit = "hPa"
    dst_val  = 0.01 * src_val
    return
  end if
  if ((varname .eq. "W") .and. (src_unit .eq. "m/s")) then
    dst_unit = "cm/s"
    dst_val  = 1.e2 * src_val
    return
  end if
  if (src_unit .eq. "K") then
    dst_unit = "dgr C"
    dst_val  = src_val-273.15
    return
  end if
  if (src_unit .eq. "m H2O") then
    dst_unit = "mm H2O"
    dst_val  = 1000.*src_val
    return
  end if
  if (((varname .eq. "RAIN_GSP") .or. (varname .eq. "SNOW_GSP") .or. (varname .eq. "RAIN_CON") .or. (varname .eq. "SNOW_CON")) .and. (src_unit .eq. "kg/m2")) then
    dst_unit = "g/m2"
    dst_val  = 1000.*src_val
    return
  end if
end

procedure sfc_get_unit_val(cdf_file, ivar, itime, istation, str_u, str_v)
local var_name, var_unit
begin
  var_name = charactertostring(cdf_file->sfcvar_name(ivar,:))
  var_unit = charactertostring(cdf_file->sfcvar_unit(ivar,:))
  transform_unit(var_name, var_unit, cdf_file->sfcvalues(itime, ivar, istation), str_u, str_v)
end

function center (zstr : string, length)
local i, blength, result
begin
  result  = ""
  blength = length - strlen(str_strip(zstr))
  do i=1,floor(blength/2.)
    result = result + " "
  end do
  result = result + str_strip(zstr)
  do i=1,ceil(blength/2.)
    result = result + " "
  end do
  return(result)
end

function lalign (zstr : string, length)
local i, blength, result
begin
  blength = length - strlen(str_strip(zstr))
  result = str_strip(zstr)
  do i=1,blength
    result = result + " "
  end do
  return(result)
end

function ralign (zstr : string, length)
local i, blength, result
begin
  blength = length - strlen(str_strip(zstr))
  result = ""
  do i=1,blength
    result = result + " "
  end do
  result = result + str_strip(zstr)
  return(result)
end

; Constants -----------------------------------------------------------

colwidth = 13; table column width

;----------------------------------------------------------------------

; open netCDF file
cdf_file = addfile(DataFileName,"r")
nstations = dimsizes(cdf_file->station_lon)

nvars    = dimsizes(cdf_file->var_name)
nsfcvars = dimsizes(cdf_file->sfcvar_name)

; prepare some data structures

VAR_PHALF    = -1
VAR_W        = -1
VAR_TKVM     = -1
VAR_TKVH     = -1
VAR_T_SO     = -1
VAR_W_SO     = -1
VAR_W_SO_ICE = -1
VAR_QV_S     = -1
VAR_RUNOFF_S = -1
VAR_RUNOFF_G = -1
VAR_ROOTDP   = -1
VAR_LAI      = -1
VAR_PLCOV    = -1
VAR_TSNOW    = -1
VAR_WSNOW    = -1
VAR_TS       = -1
VAR_WI       = -1
VAR_TG       = -1
VAR_FRESHSNW = -1
VAR_RHO_SNOW = -1
VAR_H_SNOW   = -1

do ivar=0,(nvars(0)-1)
 zname = charactertostring(cdf_file->var_name(ivar,:))
 ; model levels
 if(zname .eq. "P") then
   VAR_PMAIN = ivar
 end if
 if(zname .eq. "T") then
   VAR_T = ivar
 end if
 if(zname .eq. "QV") then
   VAR_QV = ivar
 end if
 if(zname .eq. "QC") then
   VAR_QC = ivar
 end if
 if(zname .eq. "QI") then
   VAR_QI = ivar
 end if
 if(zname .eq. "QR") then
   VAR_QR = ivar
 end if
 if(zname .eq. "QS") then
   VAR_QS = ivar
 end if
 if(zname .eq. "REL_HUM") then
   VAR_REL_HUM = ivar
 end if
 if(zname .eq. "CLC") then
   VAR_CLC = ivar
 end if
 if(zname .eq. "U") then
   VAR_U = ivar
 end if
 if(zname .eq. "V") then
   VAR_V = ivar
 end if
 ; half levels
 if(str_lower(zname) .eq. "phalf") then
   VAR_PHALF = ivar
 end if
 if(str_lower(zname) .eq. "w") then
   VAR_W = ivar
 end if
 if(str_lower(zname) .eq. "tkvm") then
   VAR_TKVM = ivar
 end if
 if(str_lower(zname) .eq. "tkvh") then
   VAR_TKVH = ivar
 end if
 ; soil levels
 if(str_lower(zname) .eq. "t_so_1") then
   VAR_T_SO = ivar
 end if
 if(str_lower(zname) .eq. "w_so_1") then
   VAR_W_SO = ivar
 end if
 if(str_lower(zname) .eq. "w_so_ice_1") then
   VAR_W_SO_ICE = ivar
 end if
end do
do ivar=0,(nsfcvars(0)-1)
 zname = charactertostring(cdf_file->sfcvar_name(ivar,:))
 if(zname .eq. "P_SFC") then
   VAR_PS = ivar
 end if
 if(zname .eq. "TCM") then
   VAR_TCM = ivar
 end if
 if(zname .eq. "TCH") then
   VAR_TCH = ivar
 end if
 if(zname .eq. "Z0") then
   VAR_Z0 = ivar
 end if
 if(zname .eq. "SHFL") then
   VAR_SHFL = ivar
 end if
 if(zname .eq. "LHFL") then
   VAR_LHFL = ivar
 end if
 if(str_lower(zname) .eq. "qv_s") then
   VAR_QV_S = ivar
 end if
 if(str_lower(zname) .eq. "runoff_s_1") then
   VAR_RUNOFF_S = ivar
 end if
 if(str_lower(zname) .eq. "runoff_g_1") then
   VAR_RUNOFF_G = ivar
 end if
 if(str_lower(zname) .eq. "la_ind") then
   VAR_LAI = ivar
 end if
 if(str_lower(zname) .eq. "pl_cov") then
   VAR_PLCOV = ivar
 end if
 if(str_lower(zname) .eq. "ro_dept") then
   VAR_ROOTDP = ivar
 end if
 if(zname .eq. "VIO3") then
   VAR_VIO3 = ivar
 end if
 if(zname .eq. "HMO3") then
   VAR_HMO3 = ivar
 end if
 if(str_lower(zname) .eq. "t_snow_1") then
   VAR_TSNOW = ivar
 end if
 if(str_lower(zname) .eq. "w_snow_1") then
   VAR_WSNOW = ivar
 end if
 if(str_lower(zname) .eq. "w_snow_1") then
   VAR_WSNOW = ivar
 end if
 if(str_lower(zname) .eq. "t_s_1") then
   VAR_TS = ivar
 end if
 if(str_lower(zname) .eq. "w_i_1") then
   VAR_WI = ivar
 end if
 if(str_lower(zname) .eq. "t_g") then
   VAR_TG = ivar
 end if
 if(str_lower(zname) .eq. "freshsnw_1") then
   VAR_FRESHSNW = ivar
 end if
 if(str_lower(zname) .eq. "rho_snow_1") then
   VAR_RHO_SNOW = ivar
 end if
 if(str_lower(zname) .eq. "h_snow_1") then
   VAR_H_SNOW = ivar
 end if
 if(zname .eq. "T2M") then
   VAR_T2M = ivar
 end if
 if(zname .eq. "TD2M") then
   VAR_TD2M = ivar
 end if
 if(zname .eq. "U10M") then
   VAR_U10M = ivar
 end if
 if(zname .eq. "V10M") then
   VAR_V10M = ivar
 end if
 if(zname .eq. "SOBT") then
   VAR_SOBT = ivar
 end if
 if(zname .eq. "SOBS") then
   VAR_SOBS = ivar
 end if
 if(zname .eq. "THBT") then
   VAR_THBT = ivar
 end if
 if(zname .eq. "THBS") then
   VAR_THBS = ivar
 end if
 if(zname .eq. "PABS") then
   VAR_PABS = ivar
 end if
 if(zname .eq. "ALB") then
   VAR_ALB = ivar
 end if
 if(zname .eq. "RAIN_GSP") then
   VAR_RAIN_GSP = ivar
 end if
 if(zname .eq. "SNOW_GSP") then
   VAR_SNOW_GSP = ivar
 end if
 if(zname .eq. "RAIN_CON") then
   VAR_RAIN_CON = ivar
 end if
 if(zname .eq. "SNOW_CON") then
   VAR_SNOW_CON = ivar
 end if
end do

; definition of table columns -----------------------------------------

print_list       = (/ VAR_PMAIN, VAR_T, VAR_QV, VAR_QC, VAR_QI, VAR_QR, VAR_QS,VAR_REL_HUM,VAR_CLC, VAR_U, VAR_V /)
print_names      = (/ "Pmain",   "T",   "QV",   "QC",   "QI"  , "QR",   "QS"  ,"REL_HUM"  ,"CLC",   "U",   "V"   /)

print_list_half  = (/ VAR_PHALF, VAR_W, VAR_TKVM, VAR_TKVH   /)
print_names_half = (/ "Phalf",   "W",   "TKVM",   "TKVH"     /)


; loop over meteogram stations:
do istation = 0,(nstations-1)
; istation = 0

; print header --------------------------------------------------------

 str_u  = ""
 str_v  = cdf_file->sfcvalues(0, 0, 0)
 str_u2 = str_u
 str_v2 = str_v
 str_u3 = str_u
 str_v3 = str_v

  print(str_get_nl())
  title_str = "STATION NAME: " + charactertostring(cdf_file->station_name(istation,:)) \
  +           "   (idx: "  + cdf_file->station_idx(istation) \
  +             ", blk: " + cdf_file->station_blk(istation) + ")"
  border_str = ""
  do i=1,strlen(title_str)
    border_str = border_str + "="
  end do
  print(tostring(border_str))
  print(tostring(title_str))
  print(tostring(border_str +str_get_nl()))

  print("Initial date : " + charactertostring(cdf_file->date(1,:)))
  print("   HSURF   ( m ): " + sprintf("%7.3f",     cdf_file->station_hsurf(istation)   ))
  print("   FR_LAND ( % ): " + sprintf("%7.3f", 100*cdf_file->station_frland(istation)  ))
  print("   LAT   ( dgr ): " + sprintf("%7.3f",     cdf_file->station_lat(istation)     ))
  print("   LON   ( dgr ): " + sprintf("%7.3f",     cdf_file->station_lon(istation)     ))
  print("   FC    (1E4/s): " + sprintf("%7.3f",1.e4*cdf_file->station_fc(istation)      ))
  print("   SOIL TYPE    : " + sprinti("%7i",       cdf_file->station_soiltype(istation)))

  ; build title string with variable names and units
  ; atmospheric (3D) variables, main levels:
  column_str = "    "+center("K", colwidth+1)
  unit_str   = "    "+center("",  colwidth+1)
  do ivar=0,(dimsizes(print_list)-1)
    scolumn = center(charactertostring(cdf_file->var_name(print_list(ivar),:)),colwidth)
    column_str = column_str + scolumn

    var_name = charactertostring(cdf_file->var_name(print_list(ivar),:))
    var_unit = charactertostring(cdf_file->var_unit(print_list(ivar),:))
    transform_unit(var_name, var_unit, str_v, str_u, str_v)
    scolumn = center(str_u,colwidth)
    unit_str = unit_str + scolumn
  end do
  column_str = column_str + center("SPEED", colwidth+1)
  var_name = "SPEED"
  var_unit = charactertostring(cdf_file->var_unit(VAR_U,:))
  transform_unit(var_name, var_unit, str_v, str_u, str_v)
  unit_str   = unit_str   + center(str_u,   colwidth+1)
  column_str = column_str + center("HML", colwidth+1)
  unit_str   = unit_str   + center("m",   colwidth+1)

  ; atmospheric (3D) variables, half levels:
  column_str_half = "    "+center("K", colwidth+1)
  unit_str_half   = "    "+center("",  colwidth+1)
  do ivar=0,(dimsizes(print_list_half)-1)
    if (print_list_half(ivar) .ge. 0) then
      scolumn = center(charactertostring(cdf_file->var_name(print_list_half(ivar),:)),colwidth)
      column_str_half = column_str_half + scolumn

      var_name = charactertostring(cdf_file->var_name(print_list_half(ivar),:))
      var_unit = charactertostring(cdf_file->var_unit(print_list_half(ivar),:))
      transform_unit(var_name, var_unit, str_v, str_u, str_v)
      scolumn = center(str_u,colwidth)
      unit_str_half = unit_str_half + scolumn
    end if
  end do
  column_str_half = column_str_half + center("HHL", colwidth+1)
  unit_str_half   = unit_str_half   + center("m",   colwidth+1)

  ntimesteps = dimsizes(cdf_file->time_step)
  max_main_levs  = cdf_file->var_nlevs(VAR_PMAIN)
  max_half_levs  = cdf_file->var_nlevs(VAR_PHALF)
  max_soil_levs = 0
  if (VAR_T_SO .ge. 0) then
    max_soil_levs  = cdf_file->var_nlevs(VAR_T_SO)
  end if

; Enable the following line to print meteogram for all time slices:
;  do itime=0,(ntimesteps-1)

    title_str = "      Actual date:     " + charactertostring(cdf_file->date(itime,:))  \
    +           " (time step: " + cdf_file->time_step(itime) + ")"
    print(str_get_nl() + tostring(title_str))
    var_name = charactertostring(cdf_file->sfcvar_name(VAR_PS,:))
    var_unit  = charactertostring(cdf_file->sfcvar_unit(VAR_PS,:))
    transform_unit(var_name, var_unit, cdf_file->sfcvalues(itime, VAR_PS, istation), str_u, str_v)
    title_str = "      PS (" + str_u + "):     " + str_v \
    +           str_get_nl()
    print(tostring(title_str))

    ; ------------------------------------------------
    ; Atmospheric (3D) variables, main levels:

    ; build a matrix containing values of sampled variables
    print(tostring(column_str))
    print(tostring(unit_str))

    values = new( (/ max_main_levs, (dimsizes(print_list)+3) /), double)
    do ilev=0,(max_main_levs-1)
      values(ilev, 0) = ilev+1
    end do
    do ivar=0,(dimsizes(print_list)-1)
      do ilev=0,(cdf_file->var_nlevs(print_list(ivar))-1)
        var_name = charactertostring(cdf_file->var_name(print_list(ivar),:))
        var_unit = charactertostring(cdf_file->var_unit(print_list(ivar),:))
        transform_unit(var_name, var_unit, cdf_file->values(itime, ilev, print_list(ivar), istation), str_u, str_v)
        values(ilev, ivar+1) = str_v
      end do ; ilev
    end do ; ivar
    ; add SPEED (from U,V)
    do ilev=0,(cdf_file->var_nlevs(VAR_U)-1)
      u = cdf_file->values(itime, ilev, VAR_U, istation)
      v = cdf_file->values(itime, ilev, VAR_V, istation)
      speed = sqrt(u*u + v*v)
      var_name = charactertostring(cdf_file->var_name(VAR_U,:))
      var_unit = charactertostring(cdf_file->var_unit(VAR_U,:))
      transform_unit(var_name, var_unit, speed, str_u, str_v)
      values(ilev, dimsizes(print_list)+1) = str_v
    end do ; ilev
    ; add heights (HML)
    do ilev=0,(cdf_file->var_nlevs(VAR_PMAIN)-1)
      values(ilev, dimsizes(print_list)+2) = cdf_file->heights(ilev, VAR_PMAIN, istation)
    end do ; ilev
    ; write matrix
    opt        = True
    opt@row    = False
    write_matrix (values, tostring(dimsizes(print_list)+3)+"f"+tostring(colwidth)+".3 ", opt)

    ; ------------------------------------------------
    ; Atmospheric (3D) variables, half levels:

    ; build a matrix containing values of sampled variables
    print(tostring(column_str_half))
    print(tostring(unit_str_half))

    values_half = new( (/ max_half_levs, (dimsizes(print_list_half)+2) /), double)
    do ilev=0,(max_half_levs-1)
      values_half(ilev, 0) = ilev+1
    end do
    do ivar=0,(dimsizes(print_list_half)-1)
      do ilev=0,(cdf_file->var_nlevs(print_list_half(ivar))-1)
        var_name = charactertostring(cdf_file->var_name(print_list_half(ivar),:))
        var_unit = charactertostring(cdf_file->var_unit(print_list_half(ivar),:))
        transform_unit(var_name, var_unit, cdf_file->values(itime, ilev, print_list_half(ivar), istation), str_u, str_v)
        values_half(ilev, ivar+1) = str_v
      end do ; ilev
    end do ; ivar
    do ilev=0,(cdf_file->var_nlevs(VAR_PHALF)-1)
      values_half(ilev, dimsizes(print_list_half)+1) = cdf_file->heights(ilev, VAR_PHALF, istation)
    end do ; ilev
    ; write matrix
    opt        = True
    opt@row    = False
    write_matrix (values_half, tostring(dimsizes(print_list_half)+2)+"f"+tostring(colwidth)+".3 ", opt)

    ; ------------------------------------------------
    ; Surface (2D) variables:

    ; Surface variables
    sfc_get_unit_val(cdf_file, VAR_TCM, itime, istation, str_u, str_v)
    print("   Surface variables:        " + lalign("TCM", colwidth) + ":"  + sprintf("%8.3f",str_v))
    sfc_get_unit_val(cdf_file, VAR_TCH, itime, istation, str_u, str_v)
    print(" " + ralign("("+str_u+")", 21) + "       " + lalign("TCH", colwidth) + ":"  + sprintf("%8.3f",str_v))
    sfc_get_unit_val(cdf_file, VAR_Z0, itime, istation, str_u, str_v)
    print(" " + ralign("("+str_u+")", 21) + "       " + lalign("Z0", colwidth) + ":"  + sprintf("%8.3f",str_v))
    sfc_get_unit_val(cdf_file, VAR_SHFL, itime, istation, str_u, str_v)
    print(" " + ralign("("+str_u+")", 21) + "       " + lalign("SHFL", colwidth) + ":"  + sprintf("%8.3f",str_v))
    sfc_get_unit_val(cdf_file, VAR_LHFL, itime, istation, str_u, str_v)
    print(" " + ralign("("+str_u+")", 21) + "       " + lalign("LHFL", colwidth) + ":"  + sprintf("%8.3f",str_v))
    sfc_get_unit_val(cdf_file, VAR_QV_S, itime, istation, str_u, str_v)
    print(" " + ralign("("+str_u+")", 21) + "       " + lalign("QV_S", colwidth) + ":"  + sprintf("%8.3f",str_v))
    if (VAR_RUNOFF_S .ge. 0) then
      sfc_get_unit_val(cdf_file, VAR_RUNOFF_S, itime, istation, str_u, str_v)
      print(" " + ralign("("+str_u+")", 21) + "       " + lalign("RUNOFF_S", colwidth) + ":"  + sprintf("%8.3f",str_v))
    end if
    if (VAR_RUNOFF_G .ge. 0) then
      sfc_get_unit_val(cdf_file, VAR_RUNOFF_G, itime, istation, str_u, str_v)
      print(" " + ralign("("+str_u+")", 21) + "       " + lalign("RUNOFF_G", colwidth) + ":"  + sprintf("%8.3f",str_v) + str_get_nl())
    end if

    ; Plants + Ozone
    if (VAR_LAI .ge. 0) then
      sfc_get_unit_val(cdf_file, VAR_LAI,  itime, istation, str_u,  str_v)
    end if
    sfc_get_unit_val(cdf_file, VAR_VIO3, itime, istation, str_u2, str_v2)
    print("   Plants:                   " + lalign("LAI", colwidth) + ":"  + sprintf("%8.3f",str_v) +  "          "   \
    +     "   Ozone:                   " + lalign("VIO3", colwidth) + ":"  + sprintf("%8.3f",str_v2))
    if (VAR_PLCOV .ge. 0) then
      sfc_get_unit_val(cdf_file, VAR_PLCOV, itime, istation, str_u,  str_v)
    end if
    sfc_get_unit_val(cdf_file, VAR_HMO3,  itime, istation, str_u2, str_v2)
    print(" " + ralign("("+str_u+")", 21) + "       " + lalign("PLCOV", colwidth) + ":"  + sprintf("%8.3f",str_v) +  "          "   \
    +     " " + ralign("("+str_u2+")", 19) + "        " + lalign("HMO3", colwidth) + ":"  + sprintf("%8.3f",str_v2))
    if (VAR_PLCOV .ge. 0) then
      sfc_get_unit_val(cdf_file, VAR_ROOTDP, itime, istation, str_u,  str_v)
    end if
    print(" " + ralign("("+str_u+")", 21) + "       " + lalign("ROOTDP", colwidth) + ":"  + sprintf("%8.3f",str_v) + str_get_nl())

    ; soil temperatures + moistures
    if ((VAR_TSNOW .ge. 0) .and. (VAR_WSNOW .ge. 0)) then
      sfc_get_unit_val(cdf_file, VAR_TSNOW, itime, istation, str_u,  str_v)
      sfc_get_unit_val(cdf_file, VAR_WSNOW, itime, istation, str_u2, str_v2)
      print("   Soil temperatures:        " + lalign("T_SNOW", colwidth) + ":"  + sprintf("%8.3f",str_v) +  "          "   \
      +     "   Soil moistures/snow      " + lalign("W_SNOW", colwidth) + ":"  + sprintf("%8.3f",str_v2))
    end if
    if ((VAR_TS .ge. 0) .and. (VAR_WI .ge. 0)) then
      sfc_get_unit_val(cdf_file, VAR_TS, itime, istation, str_u,  str_v)
      sfc_get_unit_val(cdf_file, VAR_WI, itime, istation, str_u2, str_v2)
      print(" " + ralign("("+str_u+")", 21) + "       " + lalign("T_S", colwidth) + ":"  + sprintf("%8.3f",str_v) +  "          "   \
      +     " " + ralign("("+str_u2+")", 19) + "        " + lalign("W_I", colwidth) + ":"  + sprintf("%8.3f",str_v2))
    end if
    if (VAR_FRESHSNW .ge. 0) then
      sfc_get_unit_val(cdf_file, VAR_TG, itime, istation, str_u,  str_v)
      sfc_get_unit_val(cdf_file, VAR_FRESHSNW, itime, istation, str_u2, str_v2)
      print(" " + ralign("("+str_u+")", 21) + "       " + lalign("T_G", colwidth) + ":"  + sprintf("%8.3f",str_v) +  "          "   \
      +     " " + ralign("("+str_u2+")", 19) + "        " + lalign("FRESHSNW", colwidth) + ":"  + sprintf("%8.3f",str_v2))
    end if
    if (VAR_RHO_SNOW .ge. 0) then
      sfc_get_unit_val(cdf_file, VAR_RHO_SNOW, itime, istation, str_u,  str_v)
      print(" " + ralign("("+str_u+")", 80) + "        " + lalign("RHO_SNOW", colwidth) + ":"  + sprintf("%8.3f",str_v))
    end if
    if (VAR_H_SNOW .ge. 0) then
      sfc_get_unit_val(cdf_file, VAR_H_SNOW, itime, istation, str_u,  str_v)
      print(" " + ralign("("+str_u+")", 80) + "        " + lalign("H_SNOW", colwidth) + ":"  + sprintf("%8.3f",str_v) + str_get_nl())
    end if

    do isoil=1,max_soil_levs
      sfc_get_unit_val(cdf_file, VAR_T_SO,     itime, istation, str_u,   str_v)
      sfc_get_unit_val(cdf_file, VAR_W_SO,     itime, istation, str_u2,  str_v2)
      sfc_get_unit_val(cdf_file, VAR_W_SO_ICE, itime, istation, str_u3,  str_v3)
      print(" " + ralign("("+str_u+")", 21) + "  " + lalign("T_SO("+tostring(isoil)+")", colwidth) + ":"  + sprintf("%8.3f",str_v) \
      + ralign("("+str_u2+")", 21) + "  " + lalign("W_SO("+tostring(isoil)+")", colwidth) + ":"  + sprintf("%8.3f",str_v2) \
      + ralign("("+str_u3+")", 21) + "  " + lalign("W_SO_ICE("+tostring(isoil)+")", colwidth) + ":"  + sprintf("%8.3f",str_v3))
    end do
    print("")

    ; temperatures + winds
    sfc_get_unit_val(cdf_file, VAR_T2M,     itime, istation, str_u,   str_v)
    sfc_get_unit_val(cdf_file, VAR_U10M,    itime, istation, str_u2,  str_v2)
    print("   Temperatures:             " + lalign("T2M", colwidth) + ":"  + sprintf("%8.3f",str_v) +  "          "   \
    +     "   Winds:                   " + lalign("U10M", colwidth) + ":"  + sprintf("%8.3f",str_v2))
    sfc_get_unit_val(cdf_file, VAR_TD2M,    itime, istation, str_u,   str_v)
    sfc_get_unit_val(cdf_file, VAR_V10M,    itime, istation, str_u2,  str_v2)
    print(" " + ralign("("+str_u+")", 21) + "       " + lalign("TD2M", colwidth) + ":"  + sprintf("%8.3f",str_v) +  "          "   \
    +     " " + ralign("("+str_u2+")", 19) + "        " + lalign("V10M", colwidth) + ":"  + sprintf("%8.3f",str_v2) + str_get_nl())

    ; solar + thermal radiation
    sfc_get_unit_val(cdf_file, VAR_SOBT,    itime, istation, str_u,   str_v)
    print("   Solar radiation:          " + lalign("SOBT", colwidth) + ":"  + sprintf("%8.3f",str_v) )
    sfc_get_unit_val(cdf_file, VAR_SOBS,    itime, istation, str_u,   str_v)
    sfc_get_unit_val(cdf_file, VAR_THBS,    itime, istation, str_u2,  str_v2)
    print(" " + ralign("("+str_u+")", 21) + "       " + lalign("SOBS", colwidth) + ":"  + sprintf("%8.3f",str_v) +  "          "   \
    +     " " + ralign("("+str_u2+")", 19) + "        " + lalign("THBS", colwidth) + ":"  + sprintf("%8.3f",str_v2))
    sfc_get_unit_val(cdf_file, VAR_ALB,    itime, istation, str_u2,   str_v2)
    print("   Surface albedo:           " + lalign("ALB", colwidth) + ":"  + sprintf("%8.3f",str_v2) + str_get_nl())

    ; precipitation + cloud cover
    print("   Precipitation:            " + center("amount", colwidth))
    sfc_get_unit_val(cdf_file, VAR_RAIN_GSP,    itime, istation, str_u,   str_v)
    sfc_get_unit_val(cdf_file, VAR_SNOW_GSP,    itime, istation, str_u2,  str_v2)
    print(ralign("",29) + center("("+str_u+")", colwidth))
    print("       " + lalign("RAIN_GSP:", 22) + sprintf("%8.3f",str_v))
    print("       " + lalign("SNOW_GSP:", 22) + sprintf("%8.3f",str_v2))
    sfc_get_unit_val(cdf_file, VAR_RAIN_CON,    itime, istation, str_u,   str_v)
    sfc_get_unit_val(cdf_file, VAR_SNOW_CON,    itime, istation, str_u2,  str_v2)
    print("       " + lalign("RAIN_CON:", 22) + sprintf("%8.3f",str_v))
    print("       " + lalign("SNOW_CON:", 22) + sprintf("%8.3f",str_v2))

;  end do ; itime

end do ; istation

