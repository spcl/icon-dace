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
; icon_2triangle_polygons.ncl
;
; Concepts illustrated:
;   - Plotting ICON model data
;   - Contouring 2-dimensional X, Y, Z data
;   - Using triangular meshes to create contours
;   - Drawing filled polygons on a map
;   - plot each contour level separately
;   - Turning on edges for polygons
;   - Using "getvalues" to retrieve resource values
;
;*************************************************
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;---------------------------------------------------------------
 begin
  wcStrt     = systemfunc("date")
  Model      = "ICOHAM" 
  Resolution = "R2B04L31"

  Diri      =  "~dreinert/NCL_scripts/icon_workshop/ICON_data/"
  Filename  = "MRWB4N6_DOM01_"+Resolution+"_0002.nc"

  File = addfile(Diri+Filename, "r" )
  print("File="+Filename)
  
  VarName = "DIV"                                               
  scale   = 1.e6
  level   = 20
  time    = 9 ;! attention on time stepping, here 24h
  
  var = File->$VarName$(time,level,:)   ; dims: (time,lev,cell)
  var = var*scale
  printVarSummary(var)

  print("Varscale="+scale)  
  print("max="+max(var))
  print("min="+min(var))
  print("mean="+avg(var))

;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

; Labeling  
  LeftString   = var@long_name+"  at Hybrid Level "+level+" at day "+time
  RightString  = Model+" "+Resolution
  CenterString = ""

  GraphicFile  = "icon_triangle"
  colormap  = "testcmap"

;  For use with "Manual Levels"
;  varMax = floattoint(max(var)-(max(var)-avg(var))*0.1)  ; maximum contour level
;  varMin = floattoint(min(var)+(avg(var)-min(var))*0.1 ) ; maximum contour level
;  varInt = (varMax-varMin)/10.   ; interval between contours
;  print("varmax="+varMax)
;  print("varmin="+varMin)
;  print("varint="+varInt)
  
  rad2deg = 45./atan(1.)   ; radians to degrees

;---------------------------------------------------------------
; read in the meteorological field and grid information
;---------------------------------------------------------------

  x = File->clon *rad2deg   ; cell center, lon
  y = File->clat *rad2deg   ; cell center, lat
  vlon = File->clon_vertices * rad2deg
  print("longitude min/max: " + min(vlon) + " " + max(vlon))  
  vlon = where(vlon.lt.0, vlon + 360, vlon)
  print("longitude min/max: " + min(vlon) + " " + max(vlon))  
  vlat = File->clat_vertices * rad2deg

; Note: clon and clat are longitude and latitude of triangle centers.
;       Locations of the cell corners are given by 
;       clon_vertices and clat_vertices in the nc file.

;---------------------------------------------------------------
; make plot
;---------------------------------------------------------------

  wks = gsn_open_wks("ps",GraphicFile)
  gsn_define_colormap(wks,colormap)

; Set up resources for contour/map plot.
  ResC                = True

  ResC@gsnFrame       = False
  ResC@gsnMaximize    = True

  ResC@cnFillOn       = True
  ResC@cnLinesOn      = False
  ResC@cnInfoLabelOn  = False

  FontHeight = 0.01

  ResC@tiXAxisFontHeightF        = FontHeight
  ResC@tiYAxisFontHeightF        = FontHeight 
  ResC@tmXBLabelFontHeightF      = FontHeight
  ResC@tmYLLabelFontHeightF      = FontHeight
  ResC@tmXBMajorLengthF          = FontHeight*.5
  ResC@gsnStringFontHeightF      = FontHeight - 0.002

  ResC@tmXBLabelJust  = "CenterCenter"

;  ResC@mpProjection   = "Satellite"
  ResC@mpProjection   = "CylindricalEquidistant"
  ResC@mpFillOn          = False

  ResC@mpLimitMode    = "LatLon"
  ResC@mpCenterLonF   = 110.
  ResC@mpCenterLatF   = 30.               ; choose center lat
  ResC@mpMinLonF      = 70.
  ResC@mpMaxLonF      = 150.
  ResC@mpMinLatF      =  0.
  ResC@mpMaxLatF      =  60.

  ResC@gsnMajorLonSpacing =  30. 
  ResC@gsnMinorLonSpacing =  10.
  ResC@gsnMajorLatSpacing =  15. 
  ResC@gsnMinorLatSpacing =   5. 

  ResC@mpFillOn               = False

  ResC@sfXArray                 = x  ; These are 1D arrays, so a triangular
  ResC@sfYArray                 = y  ; mesh will be created internally.

  ResC@lbLabelBarOn             = True 
  ResC@pmLabelBarHeightF        = 0.04
  ResC@pmLabelBarWidthF         = 0.5
  ResC@lbLabelAutoStride        = True
  ResC@lbTitleOn                = True         ; turn on title
  ResC@lbLabelFontHeightF       = FontHeight
  ResC@lbTitleString            = VarName+" in "+scale+" "+var@units
  ResC@lbTitlePosition          = "Bottom"
  ResC@lbTitleFontHeightF       = 0.01

;  ResC@cnLevelSelectionMode = "ManualLevels"
;  ResC@cnMinLevelValF       = varMin
;  ResC@cnMaxLevelValF       = varMax
;  ResC@cnLevelSpacingF      = varInt

  ResC@gsnSpreadColors      = True

  ResC@gsnLeftString        = LeftString 
  ResC@gsnCenterString      = CenterString
  ResC@gsnRightString       = RightString 

  ResC@cnFillMode           = "rasterfill"
  ResC@cnRasterSmoothingOn  = True
  ResC@mpGreatCircleLinesOn = True

; Create and draw the plot, but don't advance the frame.
; This is necessary in order to get the tickmarks and labelbar.
  plot = gsn_csm_contour_map(wks,var,ResC)

; Retrieve the contour levels and colors used. This information
; will be used to draw the filled triangles.
  getvalues plot@contour
   "cnLevels"     : levels
   "cnFillColors" : colors
  end getvalues
;
; Go through the vertices and create a logical array that
; indicates if the vertices are w/in the lat/lon area we're
; interested in.
;
 
  flags = new(dimsizes(var),logical);,"No_FillValue")
  flags = True

; more complicated version
;  print("var="+dimsizes(var))
;  do i = 0,dimsizes(var) - 1
;    flags(i) = where(all(vlon(i,:) .gt. ResC@mpMaxLonF) .or. \
;                     all(vlon(i,:) .lt. ResC@mpMinLonF) .or. \
;                     all(vlat(i,:) .gt. ResC@mpMaxLatF) .or. \
;                     all(vlat(i,:) .lt. ResC@mpMinLatF), \
;                     False, True)
;  end do
;  print ("out-of-bounds triangles: " + ind(flags .eq. False))

; Set up a resource list for the polygons.
  pres             = True
  pres@gsEdgesOn   = True   ; Turn on edges
  pres@gsFillIndex = 0      ; Solid fill, the default

; First draw the triangles associated with the lowest level.
  vlow = ind(var .lt. levels(0))
  do i = 0, dimsizes(vlow)-1
    pres@gsFillColor = colors(0)
    if (.not. flags(vlow(i))) then
      continue
    end if
    gsn_polygon(wks,plot,vlon(vlow(i),:),vlat(vlow(i),:),pres)
  end do
  print ("finished level 0 -- " + dimsizes(vlow) + " triangles considered")


; Now draw the triangles associated with the rest of the levels.
  print("dim of levels"+dimsizes(levels))
  do i = 0, dimsizes(levels) -2
    vind = ind(var .ge. levels(i) .and. var .lt. levels(i+1))
    do j = 0, dimsizes(vind)-1
      if (.not. flags(vind(j))) then
        continue
      end if
      pres@gsFillColor = colors(i+1)
      gsn_polygon(wks,plot,vlon(vind(j),:),vlat(vind(j),:),pres)
   end do
   print ("finished level " + i + "  -- " + dimsizes(vind) + \
          " triangles considered")
     delete(vind)
  end do

  ; last - maximum color

  i=dimsizes(levels) -2
    vind = ind( var .ge. levels(i+1))
    do j = 0, dimsizes(vind)-1
      if (.not. flags(vind(j))) then
        continue
      end if
      pres@gsFillColor = colors(i+2)
      gsn_polygon(wks,plot,vlon(vind(j),:),vlat(vind(j),:),pres)
   end do
   print ("finished level " + i + "  -- " + dimsizes(vind) + \
          " triangles considered")
     delete(vind)

  wallClockElapseTime(wcStrt, "triangle plot", 0)

  frame(wks)   ; Advance the frame

   system("gv "+GraphicFile+ " &")
end
