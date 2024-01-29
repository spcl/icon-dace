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

#  Calculate meridional overturning stream function from icon ocean model output
#
#  Author: Stephan Lorenz, MPIfMet, 07/2012 - 07/2013
#
#  Method:
#   - calculation of global/Atlantic/Pacific+Indian MOC is done in the model
#   - written in extra format
#   - simple splitting of MOC result file and plotting using cdo and ncl
#   - options: give masked file, which can be computed from the input like this:
#     cdo -setmisstoc,0.0 -div timmean_moc_0185-0190.ext timmean_moc_0185-0190.ext mocMask.nc
#   - now timmean is done in this script - no attention yet spent for incomplete years,
#     this might result in wrong MOC due to overestimation of a season

ifile=moc.bliz.r9368.2420-2470ym
ifile=moc.loc.r9558.2046.2y
ifile=moc.loc.r9513.2042.2y
ifile=moc.loc.r9558.2061.2y
ifile=$1
maskfile=$2

# check intput:
echo "Input file is '$ifile'"

basename=$(basename $ifile .ext) # ext files expected

if [ -z "$maskfile" ]; then
  cdo -f nc timmean -selvar,var777 $ifile scr_moc_glb.nc 
  cdo -f nc timmean -selvar,var778 $ifile scr_moc_atl.nc
  cdo -f nc timmean -selvar,var779 $ifile scr_moc_pac.nc
else
  cdo -f nc timmean -selvar,var777 $ifile scr_moc_glb.nc 
  cdo -div scr_moc_glb.nc -selname,var777 $maskfile tmp.nc
  mv tmp.nc scr_moc_glb.nc

  cdo -f nc timmean -selvar,var778 $ifile scr_moc_atl.nc
  cdo -div scr_moc_atl.nc -selname,var778 $maskfile tmp.nc
  mv tmp.nc scr_moc_atl.nc

  cdo -f nc timmean -selvar,var779 $ifile scr_moc_pac.nc
  cdo -div scr_moc_pac.nc -selname,var779 $maskfile tmp.nc
  mv tmp.nc scr_moc_pac.nc
fi

# run ncl-script

cat >scr_plot_moc_my.ncl << EOF
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;************************************************

begin

; if (isvar("infilename")) then
;   prefix=infilename + "_"
; else
; prefix=""
; end if
         
 f = addfile("scr_moc_atl.nc","r")
 g = addfile("scr_moc_glb.nc","r")
 h = addfile("scr_moc_pac.nc","r")

; ncl counts from SH: global 0:179; atl 30S to 90N
 moc_atl     = f->var778(0,:,60:179)  ; 120 values - 30S to 90N
 moc_glo     = g->var777(0,:, 0:179)  ; 180 values - 90S to 90N
 moc_pac_ind = h->var779(0,:,60:179)  ; 120 values - 30S to 90N

 moc_atl     = moc_atl/1000000000
 moc_glo     = moc_glo/1000000000
 moc_pac_ind = moc_pac_ind/1000000000

;************************************************************************
;  more TBD:
;   - 60S instead of -60 -> needs correct grid description
;************************************************************************

  wks = gsn_open_wks("ps","$basename")
; wks = gsn_open_wks("pdf","$basename")

  aplot = new(3,graphic);  resources for 3 plots on one page

; gsn_define_colormap(wks,"BlWhRe")            ; colormap
; gsn_define_colormap(wks,"BlGrYeOrReVi200")   ; colormap
; gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
; gsn_define_colormap(wks,"temp_diff_18lev")
; gsn_define_colormap(wks,"temp_19lev")
; gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
; gsn_define_colormap(wks,"BlueWhiteOrangeRed")
; gsn_define_colormap(wks,"testcmap")
  gsn_define_colormap(wks,"ViBlGrWhYeOrRe")

  res                      = True
  res@gsnDraw              = False             ; False: don't draw on own page
  res@gsnFrame             = False             ; False: don't advance frame

; res@cnFillMode           = "RasterFill"
  res@gsnPaperOrientation  = "Portrait"
; res@gsnPaperOrientation  = "Landscape"
; res@gsnMaximize          = True              ; no effect for part of plots

  res@gsnYAxisIrregular2Linear = True          ; True: axis-scaling is linear in Y, False: following resolution

  res@vpWidthF             = 0.85              ; change aspect ratio of plot
  res@vpHeightF            = 0.32              ; change aspect ratio of plot

  res@cnFillOn             = True              ; color for contours
; res@cnLinesOn            = False             ; contour lines
  res@lbLabelBarOn         = False             ; color label bar
  res@cnInfoLabelOn        = False             ; contour info label
  res@gsnSpreadColors      = True              ; subset of the whole colormap
; res@gsnSpreadColorStart  =  25               ; start color
; res@gsnSpreadColorEnd    = -25               ; end color
; res@cnFillMode           = "CellFill"        ; filled cells


; choice: range of levels
; res@cnLevelSelectionMode  = "ManualLevels"    ; set manual contour levels
; res@cnMinLevelValF        = -30               ; set min contour level
; res@cnMaxLevelValF        = 30                ; set max contour level
; res@cnLevelSpacingF       = 2                 ; set contour spacing

; choice: explicit levels
  res@cnLevelSelectionMode  = "ExplicitLevels"  ; set explicit contour levels
; res@cnLevels              = (/ -100, -20, -10, -5,0,5,10,20,100 /)
; res@cnLevels              = (/ -200, -150, -100, -50, -30, -20, -15, -10, -7, -5, -2, 0, \
;                                   2, 5, 7, 10, 15, 20, 50, 100, 150, 200 /)
; res@cnLevels              = (/ -100, -50, -30, -25, -20, -16, -13, -10, -7, -5, -3, -1, 0, \
;                                   1, 3, 5, 7, 10, 13, 16, 20, 25, 30, 50, 100 /)
; res@cnLevels              = (/ -100, -80, -60, -50, -40, -35, -30, -25, -22.5, -20, -18, -16, \
;                                 -14, -12, -10, -8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8, 10, \
;                                12, 14, 16, 18, 20, 22.5, 25, 30, 35, 40, 50, 60, 80, 100 /)
  res@cnLevels              = (/ -50, -40, -35, -30, -25, -22.5, -20, -18, -16, -14, -12, -10, \
     -8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22.5, 25, 30, 35, 40, 50 /)
; res@cnLevels              = (/ -27.5, -25, -22.5, -20, -18, -16, -14, -12, -10, \
;    -8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22.5, 25, 27.5 /)

  res@gsnContourZeroLineThicknessF   = 2.2      ; factor for thickness of zero contour
  res@gsnContourLineThicknessesScale = 1.5      ; factor for other lines

  res@cnLineLabelsOn        = True              ; contour line labels
  res@cnLineLabelFontHeightF= 0.012             ; contour line label height
  res@cnLineLabelInterval   = 2                 ; line label interval (default=2)
  res@cnLabelMasking        = True		; masking contour lines below line labels
; res@cnLineLabelBackgroundColor = "transparent"
  res@cnLineLabelBackgroundColor = "white"

; res@lbLabelFontHeightF    = 0.014             ; label bar font
  res@lbLabelStride         = 2                 ; labelling interval at color bar and contour lines
; res@lbLabelPosition       = "Left"  
; res@pmLabelBarOrthogonalPosF = 0.1            ; move label bar closer

  res@cnMissingValFillColor        = "gray30"    ; not defined?
  res@cnMissingValPerimOn          = True
  res@cnMissingValFillPattern      = -1          ; set the missing value fill pattern
  res@cnMissingValFillScaleF       = 0.9         ; increase the density of the fill pattern (default   = 1.0)
  res@cnMissingValPerimColor       = "black"     ; change the missing value perimeter to black
  res@cnMissingValPerimDashPattern = 1           ; set the dash pattern of the missing value perimeter to 1
  res@cnMissingValPerimThicknessF  = 2.0         ; factor for thickness of the missing value perimeter 3X

; res@gsnRightString	        = "[Sv]"         ; set below differently
  res@gsnLeftString	        = "$basename"    ; filename string
                            
  res@trYReverse                = True           ; reverse the Y-axis

  res@tiXAxisOn                 =  True        ;  false: x-axis title removed - cannot be redrawn?
; res@tiXAxisString             = "latitude"
  res@tiXAxisString             = " "          ;  text set to empty
  res@tiXAxisFontHeightF        = 0.025

  res@tiYAxisString             = "depth [m]"
  res@tiYAxisFontHeightF        = 0.025
; resources used by KMF - to be checked
; res@tmXTOn                    =  True        ;  false: turns off top ticks
; res@tmXTLabelsOn              =  True        ;  false: turns off top labels
; res@tmXBOn                    =  False       ;  false: x-axis numbers and ticks removed
  res@tmXBMode 		        = "Explicit"   ;  switch from dimension to latitude
; res@tmXBLabelsOn              =  False       ;  false: x-axis labels removed
  res@tmXBLabelFontHeightF      = 0.020        ;  x-axis labels font height
  res@tmYLLabelFontHeightF      = 0.018        ;  y-axis labels font height

; res@tiMainString	       = "Meridional Overturning Stream Function"
; res@tiMainOn                 = False       ;  default: Main title on is true

; res@vpXF = 0.1  ; not when using gsn_panel
; res@vpYF = 0.7

; *****************************************
; Plotting the 3 panels
; *****************************************

  res@tmXBValues 	       = ispan(0,120,15)   ; 0-120 values (1 deg)
  res@tmXBLabels 	       = ispan(-30,90,15)  ; spans 30S to 90N
  res@tiMainOn                 = False       ;  False: Main title removed
; res@tiMainString	       = "MOC Atlantic"
  res@gsnRightString	       = "Atlantic [Sv]"

  plot     = gsn_csm_contour(wks,moc_atl,res)
  aplot(0) = plot                            ;  creates the uppermost plot

; ==========  ========== ==========

; wks = gsn_open_wks("ps","$basename.atl")   ;  new resources are not complete for new filename

;  tmxb values and labels must correspond to dimension declared to the variables above
  delete(res@tmXBValues)
  delete(res@tmXBLabels)
  res@tmXBValues 	       = ispan(0,120,15)    ;   30S to 90N = 120x1 deg
  res@tmXBLabels 	       = ispan(-30,90,15)   ;   latitude labels
  res@gsnRightString	       = "Pacific+Ind [Sv]"

  plot     = gsn_csm_contour(wks,moc_pac_ind,res)
  aplot(1) = plot

; ==========  ========== ==========

  delete(res@tmXBValues)
  delete(res@tmXBLabels)
  res@tmXBValues 	       = ispan(0,180,30)
  res@tmXBLabels 	       = ispan(-90,90,30)
  res@gsnRightString	       = "Global [Sv]"

; delete(res@tiXAxisOn)  
; delete(res@tmXBLabelsOn)
; delete(res@tmXBOn)
; res@tiXAxisOn                = True         ;  false: x-axis title removed
  res@tiXAxisOffsetYF          =  0.010 
  res@tiXAxisString            = "latitude"
  plot     = gsn_csm_contour(wks,moc_glo,res)
  aplot(2) = plot

; *****************************************
;  Resources for 3 panels
; *****************************************

; Attention
;  - the maximize function (and other carackteristics of the whole page using gsn_panel) need own resources (resP)
  resP                           = True
  resP@gsnPaperOrientation       = "Portrait"
; resP@gsnPaperOrientation       = "Landscape"
  resP@gsnMaximize               =  True
; resP@txString                  = "MOC"
  resP@txString                  = "ICON: Meridional Overturning Stream Function"
  resP@txFontHeightF             = 0.020
; resP@tiMainString	         = "ICON: Meridional Overturning Stream Function"
; resP@tiMainOffsetYF            =  -0.3       ;  wo ist der title
  resP@tiMainOn                  = True        ;  False: Main title removed

; resP@gsnPanelYWhiteSpacePercent =  0.1       ;  <-- KMF
; resP@gsnFrame                  =  False
; resP@gsnPaperHeight            =  11.69      ;  <-- KMF
; resP@gsnPaperWidth             =   8.27      ;  <-- KMF

  resP@gsnPanelBottom            =  0.02       ;  <-- KMF
  resP@gsnPanelTop               =  0.95       ;  <-- KMF
; resP@tiXAxisOn                 =  False      ;  <-- KMF
; resP@vpXF                      =  0.05
; resP@vpYF                      =  0.95

  resP@gsnPanelLabelBar          = True        ; add common colorbar
  resP@lbLabelFontHeightF        = 0.014       ; size of labelbar font
  resP@lbAutoManage              = False
  resP@lbLabelStride             = 2                 ; skip every other label
  resP@pmLabelBarWidthF          = 0.1  
  resP@pmLabelBarHeightF         = 0.94
  resP@lbOrientation             = "vertical"
; resP@lbTopMarginF              = 0.0         ; no effect
; resP@lbBottomMarginF           = 0.1         ; no effect

; drawNDCGrid(wks)                       ;  - verkleinert Plot?
  gsn_panel(wks,(/aplot/),(/3,1/),resP)  ; now draw as one plot


end
EOF


ncl scr_plot_moc_my.ncl
rm scr_plot_moc_my.ncl
rm scr_moc_???.nc

exit

