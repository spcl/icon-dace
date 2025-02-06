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
; plot_torus_sfc.ncl 
;----------------------------------------------------------------------
; Plot torus time series.
;
; Basic Usage:
; oFile=${plotdir}/${var}_${case}"_"${turb}
; ncl -n plot_torus_sfc.ncl iFile=\"${iFile}\" oFile=\"${oFile}\" oType=\"${oType}\" \
;    varName=\"${var}\" expnum=\"${title}\"
;
; Required Parameters:
; iFile         : input file 
; varName       : name of the variable to plot
;
; Optional Parameter:
; oFile         : plot file without extension (set by the output type: oType)
; oType         : output graphics format (ps, eps, png, default: eps)
;
;----------------------------------------------------------------------
; 10/2013 : M. Koehler, DWD (martin.koehler@dwd.de)
;----------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;----------------------------------------------------------------------
; argument list and defaults
;----------------------------------------------------------------------

  if(.not. isvar("oType")) oType = "eps" end if
  if(.not. isvar("oFile")) then
    ext   = get_file_suffix(iFile,0)
    oFile = ext@fBase
  end if

;----------------------------------------------------------------------
; read data
;----------------------------------------------------------------------
  
  File  = addfile( iFile, "r" )
  times = File->time                 ; time      [time steps]
  var   = File->$varName$
  varlongname = var@long_name
  varunit     = var@units

;----------------------------------------------------------------------
; plot setup
;----------------------------------------------------------------------

  print("")
  print("Outputfile: "+str_concat((/oFile,".",oType/)))

  res = True
  res@gsnCenterString        = expnum
  res@gsnCenterStringOrthogonalPosF = 0.08
  res@gsnLeftString          = varlongname 
  res@tiXAxisString          = "Time [h]" 
  res@tiYAxisString          = "[" +  varunit + "]"
  res@tiXAxisFontHeightF     = 0.018
  res@tiYAxisFontHeightF     = 0.018
  res@gsnStringFontHeightF   = 0.020
  res@gsnStringFontHeightF   = 0.018
  res@tmXBLabelFontHeightF   = 0.016
  res@tmYLLabelFontHeightF   = 0.016
  res@tmYROn                 = False
  res@gsnMaximize            = True
  res@gsnPaperOrientation    = "portrait"
  res@xyLineThicknessF       = 3.0
  res@xyLineColor            = "blue"       ; "red","green","blue","yellow"

  wks = gsn_open_wks(oType,oFile) 

;---------------------------------------------------------------------
; Make the plot
;---------------------------------------------------------------------

  if (varName.eq."t_2m") .or. (varName.eq."u_10m") .or. (varName.eq."v_10m") then
    vals = var(:, 0, 0, 0)
  else
    vals = var(:, 0, 0)
  end if
  if ( vals(0) .eq. 0 ) then
    vals1  = vals(1:)
    times1 = (times(1:)-times(0)) * 24.0   ; [h]
  else
    vals1  = vals(0:)
    times1 = (times(0:)-times(0)) * 24.0   ; [h]
  end if

;----------------------------------------------------------------------
; line plot
;----------------------------------------------------------------------

  plot = gsn_csm_xy(wks,times1,vals1,res)

end
