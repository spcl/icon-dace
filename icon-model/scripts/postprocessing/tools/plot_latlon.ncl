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
; plot_latlon.ncl
;----------------------------------------------------------------------
; Contains
;   - Reading a NetCDF file
;   - Drawing interpolated variable on lat/lon grid
;     defined in NetCDF file
;----------------------------------------------------------------------
; 08/2011 : F. Prill, DWD
;----------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  

begin

;---------------------------------------------------------------
; user defined parameters
;---------------------------------------------------------------

  DataFileName = "plot_latlon01.nc"  ; input
  outfilename  = "plot_latlon01"
  Varname      = "T_LONLAT"

;---------------------------------------------------------------
; define map
;---------------------------------------------------------------

  pi_180 = atan(1.)/45.

  wks_type = "ps"        ; or "pdf"
  
; wks_type@wkOrientation = "landscape"

  wks  = gsn_open_wks (wks_type, outfilename)       ; open workstation
 
  gsn_define_colormap(wks,"BlWhRe")
  res  = True
  res@gsnDraw            = True
  res@gsnFrame           = True

; res@mpDataBaseVersion      = "MediumRes"          ; Default is LowRes
  res@mpProjection           = "Orthographic"
  res@mpProjection           = "CylindricalEquidistant"

; res@mpCenterLonF           = -50
; res@mpCenterLatF           =  20
; res@mpGridLatSpacingF      =  20
; res@mpGridLonSpacingF      =  20
  res@mpPerimOn              = False                ; don't draw box around map    

  res@mpGridLineDashPattern  = 2                    ; make lat/lon lines dash
  res@mpGridAndLimbOn        = True                 ; turn on lat/lon grid lines

  ; switch avoiding wrap-around globe
  res@mpGreatCircleLinesOn   = True

  res@mpFillOn               = False                ; turn on/off gray continents
  res@mpOutlineOn            = True                 ; turn on/off continental outlines

  res@tiMainFontHeightF      = 0.02


;---------------------------------------------------------------
; read in the meteorological field and grid information
;---------------------------------------------------------------

  ; open netCDF file
  cdf_file = addfile(DataFileName,"r") 

  startlev = "grid_nml:start_lev"
  grid_category = "R2B" + (cdf_file@$startlev$)
  res@tiMainString    = ":F22:ICON " + grid_category + \
  " RBF interpolation on lon-lat grid" ; add title, Font 22

  ; read (unrotated) lon-lat grid from NetCDF file:
  dim_lon = dimsizes(cdf_file->lon)
  dim_lat = dimsizes(cdf_file->lat)

  lon = cdf_file->lon(:)*pi_180
  lat = cdf_file->lat(:)*pi_180

  res@tiMainString    = res@tiMainString + \
  ", res. " + (lon(2)-lon(1))/pi_180 + "/" + (lat(2)-lat(1))/pi_180 + \
  " deg." ; add title, Font 22

  ; compute coordinates of lon-lat grid points:
  lonlat_coords = new( (/ dim_lat, dim_lon, 2 /), double)
  do i = 0,(dim_lat-1)
    do j = 0,(dim_lon-1)
      lonlat_coords(i,j,0) = (/ lon(j) /)
      lonlat_coords(i,j,1) = (/ lat(i) /)
    end do
  end do

;---------------------------------------------------------------
; draw tracer
;---------------------------------------------------------------

  ; interpolated data on lon-lat grid:
  var_lonlat = cdf_file->$Varname$(1,0,:,:)      ; dims: (time,lev,lon,lat)

  res@cnFillOn              = True             ; turn on color
  res@cnLinesOn             = False            ; turn off contour lines
  res@cnFillDrawOrder       = "PreDraw"        ; fill first (show continents)

; res@cnLevelSelectionMode  = "ManualLevels"   ; set manual contour levels
; res@cnMinLevelValF        = -100             ; set min contour level
; res@cnMaxLevelValF        = 100              ; set max contour level
; res@cnLevelSpacingF       =  20.             ; set contour spacing

; res@cnFillMode            = "CellFill"

  res@txFontHeightF = 0.015
  res@txFont = 22

; res@cnLevelSpacingF       = 0.1              ; contour spacing
  res@gsnSpreadColors       = True             ; use full colormap
; res@gsnSpreadColorStart   = 10               ; start at color 10
; res@gsnSpreadColorEnd     = 185              ; end at color 185
; res@cnLevelSpacingF       = 0.2              ; set contour spacing
; res@lbLabelStride         = 2

  res@sfXCStartV      = min(lon)/pi_180        ; where to start x
  res@sfXCEndV        = max(lon)/pi_180        ; where to end x etc
  res@sfYCStartV      = min(lat)/pi_180
  res@sfYCEndV        = max(lat)/pi_180

  res@lbBoxLinesOn    = False                  ; turn off box between colors
  res@lbOrientation   = "Vertical"

  plot = gsn_csm_contour_map(wks,var_lonlat,res)

end



