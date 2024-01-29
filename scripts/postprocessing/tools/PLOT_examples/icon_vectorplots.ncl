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
; icon_2.ncl
;
; Concepts illustrated:
;   - Plotting ICON model data
;   - Contouring one-dimensional X, Y, Z data
;   - Using triangular meshes to create contours
;   - regridding ICON-data to lat/lon grid
;   - drawing vector plots over map projection
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

  colormap1  = "testcmap"
  colormap2  = "StepSeq25"
  colormap3  = "rainbow+gray"

  scale      = 1
  scaleOMEGA = 1E3

  varMin = -30.  ; minimum contour level
  varMax =  30  ; maximum contour level
  varInt =   3  ; interval between contours

  rad2deg = 45./atan(1.)   ; radians to degrees

;---------------------------------------------------------------
; read in the meteorological field and grid information
;---------------------------------------------------------------

  Diri      =  "~dreinert/NCL_scripts/icon_workshop/ICON_data/"
  Filename1 = "MRWB4N6_DOM01_"+Resolution1+"_0002.nc"  ; input
  Filename2 = "MRWB4N6_DOM03_"+Resolution2+"_0002.nc"  ; input

  File  = addfile(Diri+Filename1, "r" )
  File2 = addfile(Diri+Filename2, "r" )

  GraphicFile  = "vector_plots"
  
  VarName1    = "U"
  VarName2    = "V"
  VarName3    = "OMEGA"
  VarName3_rf = "OMEGA"  ; rf means 'refinement'
  VarName4    = "T"

  level    = 5
  time     = 19

  var1 = File->$VarName1$(time,level,:)   ; dims: (time,lev,cell)
  var1 = var1*scale

  var2 = File->$VarName2$(time,level,:)   ; dims: (time,lev,cell)
  var2 = var2*scale

  var3 = File->$VarName3$(time,level,:)   ; dims: (time,lev,cell)
  var3 = var3*scaleOMEGA

  var3_rf = File2->$VarName3_rf$(time,level,:)   ; dims: (time,lev,cell)
  var3_rf = var3_rf*scaleOMEGA
  
  var4 = File->$VarName4$(time,level,:)   ; dims: (time,lev,cell)
  var4 = var4*scale

  ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

; Labeling  
  LeftString   = "U"
  RightString  = Resolution1 
  CenterString = ""
  TitleString  = Model+ " at level "+level+" and day "+time

  print("max="+max(var4))
  print("min="+min(var4))

; old coordinates  
  x = File->clon *rad2deg   ; cell center, lon
  y = File->clat *rad2deg   ; cell center, lat

; old coordinates rf  
  x_rf = File2->clon *rad2deg   ; cell center, lon
  y_rf = File2->clat *rad2deg   ; cell center, lat

  print("max X="+max(x))
  print("min X="+min(x))

  print("max X="+max(x_rf))
  print("min X="+min(x_rf))


; new coordinates
  
  lat= fspan(-90,90,90)
  lon= fspan(-180,180,180)

  lat@units="degrees_north"
  lon@units="degrees_east"

; new coordinates rf
  
  lat_rf= fspan(-89,89,181)
  lon_rf= fspan(-180,180,720)

  lat_rf@units="degrees_north"
  lon_rf@units="degrees_east"

  lat@units="degrees_north"
  lon@units="degrees_east"

; regridding mesh refinement
  newOMEGA_rf =  triple2grid(x_rf,y_rf,  var3_rf ,lon_rf, lat_rf, False)

; copy the attributes of the original variable  
  copy_VarAtts(var3_rf,newOMEGA_rf)

  newOMEGA_rf!0="lat"
  newOMEGA_rf!1="lon"
  newOMEGA_rf&lat=lat_rf
  newOMEGA_rf&lon=lon_rf


; regridding the data onto latlon
  newU =  triple2grid(x,y, var1 ,lon, lat, False)
  print("max newu="+max(newU))

  newU!0="lat"
  newU!1="lon"
  newU&lat=lat
  newU&lon=lon
; copy the attributes of the original variable  
  copy_VarAtts(var1,newU)
  
  printVarSummary(newU)
;  printVarSummary(lon)
  
  newV =  triple2grid(x,y,  var2 ,lon, lat, False)
  copy_VarMeta( newU, newV)
  copy_VarAtts(var2,newV)

  newOMEGA =  triple2grid(x,y,  var3 ,lon, lat, False)
  print("max newomega="+max(newOMEGA))
  copy_VarMeta( newU, newOMEGA)
  
; copy the attributes of the original variable  
  copy_VarAtts(var3,newOMEGA)
  newOMEGA@units="mPa/s"
  
  newT =  triple2grid(x,y,  var4 ,lon, lat, False)
  print("max newT="+max(newT))
  copy_VarMeta( newU, newT)
; copy the attributes of the original variable  
  copy_VarAtts(var4,newT)

; something magic to get rid of the missing values  
  guess =1
  is_cyclic = True             ; cyclic [global]
  nscan     = 1500             ; usually much less than this
  eps       = 1.e-2            ; variable dependent
  relc      = 0.6              ; relaxation coefficient
  opt       = 0 
  poisson_grid_fill( newU, is_cyclic, guess, nscan, eps, relc, opt)
  poisson_grid_fill( newV, is_cyclic, guess, nscan, eps, relc, opt)
  poisson_grid_fill( newT, is_cyclic, guess, nscan, eps, relc, opt)
;  poisson_grid_fill( newOMEGA, is_cyclic, guess, nscan, eps, relc, opt)

  print("Missingvalues="+num(ismissing(newU)))
;---------------------------------------------------------------
; make plot
;---------------------------------------------------------------

  wks = gsn_open_wks("ps",GraphicFile)
  gsn_define_colormap(wks,colormap1)
;
; Merge two pre-defined colormaps;
;
  
; Set up resources for contour/map plot.
  ResC                      = True
  FontHeight = 0.018  
  ResC@gsnFrame             = False
  ResC@gsnMaximize          = True
  ResC@gsnSpreadColors      = True
  ResC@gsnLeftString        = LeftString 
  ResC@gsnCenterString      = CenterString
  ResC@gsnRightString       = RightString
  ResC@tiMainFont           = 4
  ResC@tiMainFontHeightF    =  FontHeight   
  ResC@tiMainFontThicknessF = 4
  ResC@tiMainOffsetYF       = 0.02
  ResC@tiMainString         = TitleString
  
  ResC@cnFillOn             = True
  ResC@cnLinesOn            = False
  ResC@cnInfoLabelOn        = False
  ResC@cnLevelSelectionMode = "ManualLevels"
  ResC@cnMinLevelValF       = varMin
  ResC@cnMaxLevelValF       = varMax
  ResC@cnLevelSpacingF      = varInt
  ResC@cnFillMode           = "rasterfill"

  ResC@tiXAxisFontHeightF   = FontHeight
  ResC@tiYAxisFontHeightF   = FontHeight
  ResC@tmXBLabelFontHeightF = FontHeight
  ResC@tmYLLabelFontHeightF = FontHeight
  
  ResC@gsnStringFontHeightF = FontHeight - 0.002
  ResC@gsnAddCyclic         = False        ; add cyclic point
  ResC@gsnPolar             = "NH"         ; specify the hemisphere
  
  ResC@tmXBLabelJust  = "CenterCenter"

  ResC@mpPerimOn         = False             ; turn off box around plot
  ResC@mpFillOn          = False
  ResC@mpLimitMode       = "LatLon"
  ResC@mpCenterLonF      =   0.
  ResC@mpMinLatF         =  0.
  ResC@mpMaxLatF         =  90.
  
  ResC@lbLabelBarOn       = True 
  ResC@lbLabelAutoStride  = True
  ResC@lbLabelFontHeightF = FontHeight
  ResC@lbTitleFontHeightF = 0.015
  ResC@lbTitlePosition    = "Bottom"
    
  ResC@pmLabelBarHeightF        = 0.07
  ResC@pmLabelBarWidthF         = 0.7
  ResC@pmLabelBarOrthogonalPosF = 0.06

; Create and draw the plot, but don't advance the frame.
; This is necessary in order to get the tickmarks and labelbar.

; Plot the wind

; zonal component

  ResC@lbTitleString       = newU@long_name+" in "+newU@units
  plotU = gsn_csm_contour_map_polar(wks,newU,ResC)

  print("newU plotted")
  frame(wks)

; meridional component
  ResC@gsnLeftString        = "V"
  ResC@lbTitleString       = newV@long_name+" in "+newV@units
  plotV = gsn_csm_contour_map_polar(wks,newV,ResC)

  print("newV plotted")
  frame(wks)

; vertical component
  delete(ResC@cnMinLevelValF  )
  delete(ResC@cnMaxLevelValF )     
  delete(ResC@cnLevelSpacingF )    

 
; automatically create nice min/max/ci values for blue/red colortable
;  res= True
  symMinMaxPlt (newOMEGA,20,False,ResC)

; plot without Poisson-Fill
  ResC@gsnLeftString = "without Poisson Fill"
  ResC@lbTitleString = newOMEGA@long_name+" in "+newOMEGA@units

  plotOMEGA = gsn_csm_contour_map_polar(wks,newOMEGA,ResC)
  frame(wks)
  
  ResC@gsnLeftString = ""
  ResC@gsnRightString = ""
  ResC@lbTitleString       = newOMEGA@long_name+" in "+newOMEGA@units
  
  ; plot high resolution field without Poisson-Fill

  delete(ResC@gsnPolar) 
  delete(ResC@mpLimitMode)
  delete(ResC@mpMinLatF)
  delete(ResC@mpMaxLatF)
  delete(ResC@mpCenterLonF)
  delete(ResC@mpFillOn)
  delete(ResC@mpPerimOn)
  
  ResC2               = ResC
  ResC2@cnFillMode    = "cellfill" 
  ResC2@gsnDraw       = False 
;  print(newOMEGA_rf@_FillValue)
;  print(dimsizes(ind(ndtooned(ismissing(newOMEGA_rf)))))
  ResC2@cnMissingValFillColor = -1
  ResC2@cnMissingValPerimOn = True
  ResC2@cnMissingValPerimThicknessF = 3
  plotOMEGA_rf = gsn_csm_contour(wks,newOMEGA_rf,ResC2)
; draw(plotOMEGA_rf)
;  print("newOMEGA plotted")
;  frame(wks)

  ; plot with Poisson-Fill
  poisson_grid_fill( newOMEGA, is_cyclic, guess, nscan, eps*(scale/1000.), relc, opt)

  ResC@gsnLeftString  = "Poisson OMEGA"
  ResC@gsnRightString = Resolution1+" & "+Resolution2
  ResC@lbTitleString  = newOMEGA@long_name+" in "+newOMEGA@units
  plotOMEGAP = gsn_csm_contour_map_polar(wks,newOMEGA,ResC)

  overlay(plotOMEGAP,plotOMEGA_rf)
  draw(plotOMEGAP)
  frame(wks)

  print("now plot different hemispheres")
  
  
; Vector Plots of SH and NH
  
  vcres                     = True      ; plot mods desired
  vcres@gsnFrame            = False
  vcres@gsnDraw             = False
  vcres@gsnAddCyclic        = False        ; add cyclic point
  vcres@gsnSpreadColors     = True      ; use full colormap
  vcres@gsnSpreadColorEnd   = 193       ; last color to use
  vcres@gsnSpreadColorStart = 6         ; first color to use
  vcres@lbLabelStride       = 2         ; plot every other colar bar label
  i = NhlNewColor(wks,0.7,0.7,0.7)      ; add gray to colormap for continents
                                        ; w/o this, they are white, and you
                                        ; can not see them.

  vcres@vcRefMagnitudeF         = 5.0              ; make vectors larger
  vcres@vcRefLengthF            = 0.040            ; ref vec length
  vcres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  vcres@vcMinDistanceF          = 0.030            ; thin out vectors

  vcres@tiMainString            = ""
  vcres@gsnLeftString           = "wind vectors"
  vcres@gsnRightString          = " "

 ; show vertical velocity field and overlay velocity vectors 
 ; vcres@gsnLeftString        = "vector plot"
  
  vplot=gsn_csm_vector(wks,newU,newV,vcres)  ; create plot

  overlay(plotOMEGAP,vplot)
  draw(plotOMEGAP)           ; draw plot
  frame(wks)            ; advance frame

  vcres@gsnPolar       = "NH"                          ; specify the hemisphere
  vcres@tiMainString   = "NH"
  plot=gsn_csm_vector_map_polar(wks,newU,newV,vcres)  ; create plot
 
  draw(plot)
  frame(wks)
 
  t= newT

  ; merge colormap 1 and 2
  gsn_define_colormap(wks,colormap3)
;  gsn_merge_colormaps(wks,colormap3,colormap2)
  vcres@gsnSpreadColors     = True
  vcres@gsnSpreadColorStart = 30      ; last color to use
  vcres@gsnSpreadColorEnd   = 235     ; first color to use

  vcres@lbLabelBarOn        = True 
  vcres@lbLabelAutoStride   = True
  vcres@lbLabelFontHeightF  = FontHeight
  vcres@lbTitleFontHeightF  = FontHeight  ;0.016
  vcres@lbTitlePosition     = "Bottom"
  vcres@lbTitleString       = newT@long_name+" in "+newT@units

  vcres@pmLabelBarHeightF        = 0.08 
  vcres@pmLabelBarWidthF         = 0.7
  vcres@pmLabelBarOrthogonalPosF = 0.06
  plot=gsn_csm_vector_scalar_map_polar(wks,newU,newV,t,vcres)  ; create plot
 
  wallClockElapseTime(wcStrt, "triangle plot", 0)

  draw(plot)
  frame(wks)   ; Advance the frame
  system("gv "+GraphicFile +" &")
end
