; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

; Hui Wan (MPI-M, 2009-03)
;------------------------------------------------------------------------------ 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;------------------------------------------------------------------------------ 
begin

  PlotFormat = getenv("plot_file_format")
  DataPath    = getenv("DataPath")
  PlotPath    = getenv("PlotPath")
  ExpName     = getenv("ExpName")
  RightString = getenv("ConfigStr")

  day = (/10,30/)
  ndays = dimsizes(day)
  plot  = new( ndays,graphic )

  FontHeight = 0.02
  dashps = (/2,1,0/)
  thknes = (/2,2,2/)
  colors = (/2,17,30/)

;------------------------------
; the k^{-3} reference line
;------------------------------

  duml  = plot   ; ref. line
  dumt  = plot   ; label of the ref. line

  xl = (/15.,300./)
  yl = 5*(xl/15)^(-3) ; make it 6 @ k=15

  ResPL = True
  ResPL@gsLineColor           = 40  ; gray
  ResPL@gsLineThicknessF      = 2
  ResPL@gsLineDashPattern     = 7 

  ResT               = True
  ResT@txFontHeightF = FontHeight 
  ResT@txFontColor   = 40

  Text = "n~S~-3~N~"
  xt = 100   ; location of the text
  yt = 0.1   ; location of the text

;------------------------------
; inquire array shape
;------------------------------

  SpecDataFile = DataPath+ExpName+"_KE-spectrum_day"+day(0)+".nc" 
  File = addfile( SpecDataFile,"r" )

  nspc = dimsizes(File->ke(0,:)) -1   ; exclude the global mean
  nlev = dimsizes(File->ke(:,0))
   lev = File->lev

  ke   = new( (/3,nspc/),double )
  wn   = ispan(1,nspc,1)

;------------------------------
; xy-plot: resources
;------------------------------

  ResL = True

  ResL@gsnDraw            = False
  ResL@gsnFrame           = False

  ResL@tmYLLabelFontHeightF    = FontHeight
  ResL@tmXBLabelFontHeightF    = FontHeight
  ResL@tmXBMinorOn             = True 
  ResL@tmXBMinorPerMajor       = 8
  ResL@tmXTMinorPerMajor       = 8
  ResL@tmYRMinorPerMajor       = 8
  ResL@tmYLMinorPerMajor       = 8
  ResL@tmXTOn                  = True 

  ResL@gsnLeftStringFontHeightF  = FontHeight
  ResL@gsnRightStringFontHeightF = FontHeight
  ResL@tiXAxisFontHeightF        = FontHeight
  ResL@tiYAxisFontHeightF        = FontHeight

  ResL@tiXAxisString = "Spectral wavenumber n"
  ResL@tiYAxisString = "Kinetic energy (m~S~2~N~ s~S~-2~N~)"

  ResL@pmLegendDisplayMode    = "Always"
  ResL@pmLegendSide           = "Top"
  ResL@pmLegendParallelPosF   = 0.2
  ResL@pmLegendOrthogonalPosF = -1.
  ResL@pmLegendWidthF         = 0.06
  ResL@pmLegendHeightF        = 0.1
  ResL@lgPerimOn              = False 
  ResL@lgLabelFontHeightF     = FontHeight 
  ResL@xyExplicitLegendLabels = "  "+(/"rotational","divergent","total"/)

  ResL@xyLineColors         = colors
  ResL@xyMonoDashPattern    = False
  ResL@xyDashPatterns       = dashps
  ResL@xyMonoLineThickness  = False
  ResL@xyLineThicknesses    = thknes
 
  ResL@trYLog        = True 
  ResL@trYMinF       = 1e-7 
  ResL@trYMaxF       = 1e3 
  ResL@tmYROn        = True
  ResL@tmYRLabelsOn  = False 
  ResL@trXLog        = True 
  ResL@trXMinF       = 1 
  ResL@trXMaxF       = 500 

;----------------------------------------------------------------------------
; Draw the kinetic energy spectra: each vertical level is in a separate file
;----------------------------------------------------------------------------

do ilev=0,nlev-1

  wks = gsn_open_wks(PlotFormat,PlotPath+ExpName+"_KE-spectrum_"+lev(ilev)/100.+"hPa")
  gsn_merge_colormaps(wks,"StepSeq25","gsdtol")

  do jday = 0,ndays-1 

     SpecDataFile = DataPath+ExpName+"_KE-spectrum_day"+day(jday)+".nc" 
     File = addfile( SpecDataFile,"r" )

   ; read data

     ke(2,:) = File->ke (ilev,1:)
     ke(1,:) = File->div(ilev,1:)
     ke(0,:) = File->rot(ilev,1:)

   ; draw curves 

     ResL@gsnLeftString = lev(ilev)/100. +" hPa, " +"day "+ day(jday) 
     ResL@gsnRightString = RightString
     ResL@gsnLeftStringOrthogonalPosF  = 0.01
     ResL@gsnRightStringOrthogonalPosF = 0.01

     plot(jday) = gsn_csm_xy(wks,wn,ke,ResL)
     duml(jday) = gsn_add_polyline(wks,plot(jday),xl,yl,ResPL)
     dumt(jday) = gsn_add_text(wks,plot(jday),Text, xt,yt,ResT)
  
  end do

  ResP = True
  ResP@gsnPanelXWhiteSpacePercent = 5
  ResP@gsnPanelYWhiteSpacePercent = 5

  gsn_panel(wks,plot,(/ndays/2,2/),ResP)

end do

end
