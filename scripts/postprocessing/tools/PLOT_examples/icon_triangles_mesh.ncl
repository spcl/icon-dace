; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;*************************************************
; icon_xy.ncl
;
; Concepts illustrated:
;   - Plotting ICON model data
;   - Contouring one-dimensional X, Y, Z data
;   - Using triangular meshes to create contours
;   - overlaying the refined grid onto the coarser one
;
;*************************************************
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;---------------------------------------------------------------
 begin
  wcStrt     = systemfunc("date")
  Model      = "ICOHAM" 
  Resolution1 = "R2B04L31" 
  Resolution2 = "R2B06L31" 
  ConfigStr  = "D1 spr0.90"

  Diri      =  "~dreinert/NCL_scripts/icon_workshop/ICON_data/"
  FileName1 = "MRWB4N6_DOM01_"+Resolution1+"_0002.nc"  ; input
  FileName2 = "MRWB4N6_DOM03_"+Resolution2+"_0002.nc"  ; input

  GraphicFile  = "icon_triang_mesh"
  
  colormap  = "BlueDarkRed18"

  rad2deg = 45./atan(1.)   ; radians to degrees

;---------------------------------------------------------------
; read in the meteorological field and grid information
;---------------------------------------------------------------

  File  = addfile(Diri+FileName1, "r" )
  File2 = addfile(Diri+FileName2, "r" )

  print("File1="+FileName1)
  print("File2="+FileName2)
  
  VarName1 = "U"
  VarName2 = "OMEGA"
  VarName3 = "OMEGA"

  level    = 5
  time     = 19
  scale    = 1e3
   
  var = File->$VarName1$(time,level,:)   ; dims: (time,lev,cell)
  var = var

  var2 = File->$VarName2$(time,level,:)   ; dims: (time,lev,cell)
  var2 = var2*scale
  var2@units="mPa"
  printVarSummary(var2)
  
  var3 = File2->$VarName3$(time,level,:)   ; dims: (time,lev,cell)
  var3 = var3*scale
  var3@units="mPa"
  printVarSummary(var3)
 
  print("max="+max(var2))
  print("min="+min(var2))
  
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

; Labeling
  
  LeftString   = Resolution1
  RightString  = ""
  CenterString = Model
  TitleString  = var2@long_name+" at Hybrid Level "+level+" at day "+time
  FontHeight   = 0.02
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
; coordinates of global coarse grid
  x = File->clon *rad2deg   ; cell center, lon
  y = File->clat *rad2deg   ; cell center, lat

; coordinates of refined grid  
  x_rf = File2->clon *rad2deg   ; cell center, lon
  y_rf = File2->clat *rad2deg   ; cell center, lat

;---------------------------------------------------------------
; make plot
;---------------------------------------------------------------

  wks = gsn_open_wks("ps",GraphicFile)
  gsn_define_colormap(wks,colormap)

; Set up resources for contour/map plot.
  ResC                = True

  ResC@mpProjection         = "Orthographic"
  ResC@mpPerimOn            = False             ; turn off box around plot
  ResC@mpFillOn             = False
  ResC@mpCenterLonF         = 90.
  ResC@mpMinLatF            = -90.
  ResC@mpMaxLatF            =  90.
  ResC@mpGreatCircleLinesOn = True
   
  ResC@gsnFrame             = False
  ResC@gsnMaximize          = True
  ResC@gsnSpreadColors      = True
  ResC@gsnStringFontHeightF = FontHeight 
  ResC@gsnLeftString        = LeftString 
  ResC@gsnCenterString      = CenterString
  ResC@gsnRightString       = RightString 
  ResC@cnFillOn             = True
  ResC@cnLinesOn            = False
  ResC@cnInfoLabelOn        = False
  ResC@cnFillMode           = "Areafill"
  ResC@cnRasterSmoothingOn  = True

  ResC@sfXArray       = x  ; These are 1D arrays, so a triangular
  ResC@sfYArray       = y  ; mesh will be created internally.

  ResC@tiMainFont           = 4
  ResC@tiMainFontHeightF    =  FontHeight   
  ResC@tiMainFontThicknessF = 4
  ResC@tiMainOffsetYF       = 0.02
  ResC@tiMainString         = TitleString
   
  ResC@pmLabelBarHeightF        = 0.07
  ResC@pmLabelBarWidthF         = 0.7
  ResC@pmLabelBarOrthogonalPosF = 0.05
  ResC@lbLabelBarOn             = True 
  ResC@lbLabelAutoStride        = True
  ResC@lbTitlePosition          = "Bottom"
  ResC@lbTitleString            = var2@long_name+" in "+var2@units
  ResC@lbTitleFontHeightF       = 0.015
  
; automatically create nice min/max/ci values for blue/red colortable
  symMinMaxPlt (var2,20,False,ResC)
  
; Create and draw the plot, but don't advance the frame.
; This is necessary in order to get the tickmarks and labelbar.

; plot omega
  plot = gsn_csm_contour_map(wks,var2,ResC)

  frame(wks)
  
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ; plot high resolution field 
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
; ATTENTION  
; if you want to overlay plots, the overlaying plot must be ONLY a contour plot

; deleting ressources which do not belong to the "gsn_csm_contour"
  delete( ResC@sfXArray )
  delete( ResC@sfYArray )
  delete(ResC@mpMinLatF)
  delete(ResC@mpMaxLatF)
  delete(ResC@mpCenterLonF)
  delete(ResC@mpFillOn)
  delete(ResC@mpPerimOn)
  delete( ResC@mpProjection)
  delete(ResC@mpGreatCircleLinesOn) 
  
  ResC2                             = ResC
  ResC2@sfXArray                    = x_rf  ; These are 1D arrays, so a triangular
  ResC2@sfYArray                    = y_rf  ; mesh will be created internally.
  ResC2@cnFillMode                  = "Areafill" 
  ResC2@gsnDraw                     = False 
  ResC2@cnMissingValFillColor       = -1
  ResC2@cnMissingValPerimOn         = True
  ResC2@cnMissingValPerimThicknessF = 2

  ResC2@gsnLeftString   =""
  ResC2@gsnRightString  = Resolution2
  ResC2@gsnCenterString = ""
  
  plotOMEGA_rf = gsn_csm_contour(wks,var3,ResC2)
  overlay(plot,plotOMEGA_rf)
  draw(plot)
  frame(wks)   ; Advance the frame
  wallClockElapseTime(wcStrt, "triangle plot", 0)
 
  system("gv "+GraphicFile +" &")
end
