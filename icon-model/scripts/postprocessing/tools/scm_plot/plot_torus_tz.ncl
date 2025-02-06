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
; plot_torus_tz.ncl 
;----------------------------------------------------------------------
; Plot torus profile over time.
;
; Basic Usage:
; oFile=${plotdir}/${var}_${case}"_"${turb}
; ncl -n plot_torus_tz.ncl iFile=\"${iFile}\" oFile=\"${oFile}\" oType=\"${oType}\" \
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
; setup
;----------------------------------------------------------------------

  ntop   = 75
  nbot   = 90

;----------------------------------------------------------------------
; argument list and defaults
;----------------------------------------------------------------------

  if(.not. isvar("oType")) oType = "eps" end if
  if(.not. isvar("oFile")) then
    ext   = get_file_suffix(iFile,0)
    oFile = ext@fBase
  end if
  if (.not. isvar("colormap"))    colormap    = "BlueDarkRed18" end if

;----------------------------------------------------------------------
; read data
;----------------------------------------------------------------------
  
  File = addfile( iFile, "r" )
  if (varName.eq."w") then
    lev = File->height_2             ;lev_2 or ilev
    zz  = File->z_ifc
  else   
    lev = File->height               ;lev or mlev]
    zz  = File->z_mc
  end if
  nlevs = doubletoint(max(lev))
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
  res@trYReverse             = True		; reverse the Y-axis
  res@trYMinF                = ntop
  res@trYMaxF                = nbot
  res@gsnCenterString        = expnum
  res@gsnCenterStringOrthogonalPosF = 0.08
  res@gsnLeftString          = varlongname 
  res@tiXAxisString          = "Time [h]" 
  res@tiYAxisString          = "Model Level"
  res@tiXAxisFontHeightF     = 0.018
  res@tiYAxisFontHeightF     = 0.018
  res@gsnStringFontHeightF   = 0.020
  res@gsnStringFontHeightF   = 0.018
  res@tmXBLabelFontHeightF   = 0.016
  res@tmYLLabelFontHeightF   = 0.016
  res@tmYROn                 = False
  res@gsnMaximize            = False
  res@gsnPaperOrientation    = "portrait"
  res@cnLinesOn              = False             ; optional lines
  res@cnFillMode             = "RasterFill"      ; optional raster
  res@cnFillOn               = True              ; turn on color fill
  res@gsnSpreadColors        = True              ; use full colormap
  ; label bar resources
  res@lbOrientation          = "vertical"        ; vertical label bar
  res@pmLabelBarWidthF       = 0.24              ; horizontal shift
  res@lbBoxMinorExtentF      = 0.15              ; horizontal width
  res@vpXF                   = 0.13              ; viewport left margin (def: 0.2)
  res@vpWidthF               = 0.62              ; viewport width (def: 0.6)
  
  wks = gsn_open_wks(oType,oFile) 
  gsn_define_colormap(wks,colormap)  

;----------------------------------------------------------------------
; create an independent axis with level heights
;----------------------------------------------------------------------

  ystride      = 1 ; subdivisions of right Y-Axis labels
  expl_labels2 = new(nlevs, string)
  do i=0,(nlevs-1)
    expl_labels2 = sprintf("%3.0f m", zz)
  end do
  xy3 = create "xy3" xyPlotClass wks
      "trYReverse"              : "True"
      "trYMaxF"                 :  nbot
      "trYMinF"                 :  ntop
      "tmYROn"                  : "True"
      "tmYLOn"                  : "False"
      "tmYUseLeft"              : "False"
      "tmYRLabelsOn"            : "True"
      "tmYRMajorLengthF"        : .01
      "tmYRMajorOutwardLengthF" : .01
      "tmXBOn"                  : "False"
      "tmXTOn"                  : "False"
      "tmYLOn"                  : "False"
      "tmYROn"                  : "True"
      "tmYRMode"                : "Explicit"
      "tmYRValues"              : lev(0:(nlevs-1):ystride)
      "tmYRLabels"              : expl_labels2(0:(nlevs-1):ystride)
      "tmYRLabelJust"           : "CenterRight"
      "tiYAxisSide"             : "Right"
      "tmYRLabelFontHeightF"    : 0.01
      "tiXAxisFontHeightF"      : 0.01
      "tiYAxisFontHeightF"      : 0.03
      "tmYRMinorOn"             : "False"
      "tmYLMinorOn"             : "False"
      "tmXBMinorOn"             : "False"
      "tmXTBorderOn"            : "False"
      "tmXBBorderOn"            : "False"
      "tmYLBorderOn"            : "False"
      "vpXF"                    : 0.13    
      "vpWidthF"                : 0.62    
  end create
  draw(xy3)

;---------------------------------------------------------------------
; Make the plot
;---------------------------------------------------------------------
  
  var0 = var(:, res@trYMinF-1:res@trYMaxF-1, 0, 0)
  z0   = lev(   res@trYMinF-1:res@trYMaxF-1)
  var1 = transpose(var0)

  res@sfXArray  = (times-times(0)) * 24.0      ; [h]
  res@sfYArray  = z0

;----------------------------------------------------------------------
; contour plot
;----------------------------------------------------------------------

  plot = gsn_csm_contour(wks,var1,res)

end
