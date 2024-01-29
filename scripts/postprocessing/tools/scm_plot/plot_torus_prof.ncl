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
; plot_torus_prof.ncl 
;----------------------------------------------------------------------
; Plot torus profile.
;
; Basic Usage:
; oFile=${plotdir}/${var}_${case}"_"${turb}
; ncl -n plot_torus_prof.ncl iFile=\"${iFile}\" oFile=\"${oFile}\" oType=\"${oType}\" \
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
  isteps = (/0,3,6,8/)
  nsteps = 4

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
  res@tiXAxisString          = "[" +  varunit + "]"
  res@tiYAxisString          = "Model Levels"
  res@tiXAxisFontHeightF     = 0.018
  res@tiYAxisFontHeightF     = 0.018
  res@gsnStringFontHeightF   = 0.020
  res@gsnStringFontHeightF   = 0.018
  res@tmXBLabelFontHeightF   = 0.016
  res@tmYLLabelFontHeightF   = 0.016
  res@tmYROn                 = False
  res@gsnMaximize            = False
  res@gsnPaperOrientation    = "portrait"
  res@xyLineThicknessF       = 3.0
  res@xyComputeXMin          = True
  res@gsnDraw                = False
  res@gsnFrame               = False
  
  wks = gsn_open_wks(oType,oFile) 

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
  end create
  draw(xy3)

;---------------------------------------------------------------------
; Make the plot
;---------------------------------------------------------------------
  
  var0 = var(isteps(0), res@trYMinF-1:res@trYMaxF-1, 0, 0)
  var1 = var(isteps(1), res@trYMinF-1:res@trYMaxF-1, 0, 0)
  var2 = var(isteps(2), res@trYMinF-1:res@trYMaxF-1, 0, 0)
  var3 = var(isteps(3), res@trYMinF-1:res@trYMaxF-1, 0, 0)
  z0 = lev(res@trYMinF-1:res@trYMaxF-1)

  colors = (/"red","green","blue","yellow"/)

  min1=(/min(var0),min(var1),min(var2),min(var3)/)
  max2=(/max(var0),max(var1),max(var2),max(var3)/)
  res@trXMinF = min(min1)
  res@trXMaxF = max(max2) 

  res@xyLineColor      = colors(0)
  plot0 = gsn_csm_xy(wks,var0,z0,res)
  res@xyLineColor      = colors(1)
  plot1 = gsn_csm_xy(wks,var1,z0,res)
  res@xyLineColor      = colors(2)
  plot2 = gsn_csm_xy(wks,var2,z0,res)
  res@xyLineColor      = colors(3)
  plot3 = gsn_csm_xy(wks,var3,z0,res)

  overlay(plot0,plot1)
  overlay(plot0,plot2)
  overlay(plot0,plot3)

;----------------------------------------------------------------------
; Attach a legend
;----------------------------------------------------------------------

  lgres                    = True
  lgres@lgLineColors       = colors
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = 0.08           ; legend label font thickness
  lgres@vpWidthF           = 0.15            ; width of legend (NDC)
  lgres@vpHeightF          = 0.10           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgLineThicknessF   = 3.0 
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0

  labels = " t = "+ ( ( times(isteps) - times(0) )* 24.0 ) + "h "
  legend = gsn_create_legend (wks, nsteps, labels,lgres)

  amres                  = True
  amres@amJust           = "TopRight"       ; box location
  amres@amParallelPosF   =  0.5             ; Move legend to right
  amres@amOrthogonalPosF = -0.5             ; Move legend down.
  annoid = gsn_add_annotation(plot0,legend,amres)  ; add legend to plot

;----------------------------------------------------------------------
; line plot
;----------------------------------------------------------------------

  draw(plot0)

end
