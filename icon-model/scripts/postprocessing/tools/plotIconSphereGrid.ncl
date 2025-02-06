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
;
; plots one grid, coloring the edges depending on given variable
;
; usage:
; ncl plotIconSphereGrid.ncl 'grid="<grid file name>"' 'var="<edge var name>"' 'minValue=<min value>' 'maxValue=<max value>' 
;
; Leonidas Linardakis, MPI-M 2010-12-23
;
; TODO: 1. plot nested grids, 2. plot colored cells,  3. integrate the plotOceanGrid.ncl here
;
;*************************************************
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;---------------------------------------------------------------
 
;---------------------------------------
; some parameters
 rad2deg = 57.295779513082323110286093015487255497
;---------------------------------------

;************************************************
; set min-max cenetered at Lat=0, res@mpCenterLonF
;************************************************
procedure getCenteredMap(res:logical, latRadius:float, lonRadius:float)
begin 
  res@mpMinLatF      = res@mpCenterLatF - latRadius
  res@mpMaxLatF      = res@mpCenterLatF + latRadius
  res@mpMinLonF      = res@mpCenterLonF - lonRadius
  res@mpMaxLonF      = res@mpCenterLonF + lonRadius
end
;************************************************



;************************************************
; checks if point is inside min-max drawing area
; Note: This assumes that the min-max drawing area is correct
;************************************************
function pointIn(inlon:double,inlat:double,res:logical)
begin
  lon=inlon
  lat=inlat
;  if (lon .gt. 180.) then
;     lon = lon - 360.0
;  end if
;  if (lon .lt. 0.) then
;     lon = lon + 360.0
;  end if
;  inPointShift = 1 * abs(res@mpMaxLatF) + 1.0
  inPointShift = 2.0
 
  if (lon .lt. (res@mpMinLonF - inPointShift)  ) then
     return(False)
  end if
  if (lon .gt. (res@mpMaxLonF + inPointShift)) then
     return(False)
  end if
  
  if (lat .lt. (res@mpMinLatF - inPointShift) ) then
     return(False)
  end if
  if (lat .gt. (res@mpMaxLatF + inPointShift) ) then
     return(False)
  end if
  
  return(True)

end
;************************************************



 begin

 setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize": 132000000
 end setvalues

  wcStrt     = systemfunc("date")
  colormap  = "BlueDarkRed18"
  colormap  = "rainbow"
  colormap  = "ncview_default"
  edgeThickness = 2.5 
  gridFolder=""
  
  ; graphicsFormat="png"
  ; graphicsFormat@wkWidth = 1200
  ; graphicsFormat@wkHeight =1200

  ; "pdf"
  graphicsFormat="ps"
  if (isvar("format")) then
    graphicsFormat = format
  end if

  ; "Satellite" "Robinson" "Orthographic" "WinkelTripel"
  if (.not. isvar("projection")) then
    projection = "Satellite"
  end if
     
  if( isvar("folder")) then
     gridFileName = folder
  end if

  if( isvar("grid")) then
     gridFileName = grid
  else
     print("grid is not given")
     exit     
  end if
  
  gridFile = addfile(gridFolder+gridFileName+".nc", "r" )
  
  if( isvar("stats")) then
     statsFileName = stats
     statsFile = addfile(gridFolder+statsFileName+".nc", "r" )
  else
     statsFileName=gridFileName
     statsFile=gridFile
  end if
  
  if( isvar("var")) then
     varName = var
  else
     varName="edge_length"
  end if
  
  if( isvar("thickness")) then
     edgeThickness = thickness
  end if

  graphicsFile  = statsFileName+"."+varName 
  
  dims  = getfilevardims(statsFile,varName)
  sizes = filevardimsizes(statsFile,varName)
 ; plotVar = statsFile->$varName$(:)   ; dims: (cell)
  
  ; Print variable information.
  print( "Plotting " + varName + ": " + dimsizes(sizes) + " Dimensions:")
  if(.not.any(ismissing(dims))) then
    do j = 0, dimsizes(dims) -1
    print( j + ") " + dims(j) + ": " + sizes(j))
    end do
  end if
   
;---------------------------------------------------------------
; make plot
;---------------------------------------------------------------
; Labeling

  LeftString   = ""
  RightString  = ""
  TitleString  = gridFileName
  CenterString = varName
  FontHeight   = 0.01
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  wks = gsn_open_wks(graphicsFormat,graphicsFile)
  gsn_define_colormap(wks,colormap)
;  gsn_reverse_colormap(wks)
  
  setvalues wks            
;    "wkColorMap"        : "gsltod"  
;    "wkForegroundColor" : (/0.,0.,0./)  
    "wkBackgroundColor" : (/1.,1.,1./)
;    "wkBackgroundColor" : (/-1.,-1.,-1./)
  end setvalues

; Set up resources for contour/map plot.
  res                = True

  res@mpProjection         = projection
  res@mpPerimOn            = False             ; turn off box around plot
  res@mpFillOn             = False
  centerLon = -25.0  ;  view for new Stommel basin
  centerLat =  40.0   ;  view for new Stommel basin
  centerLon =  -40.0
  centerLat =  0.0
  latRadius = 90.0
  lonRadius = 90.0    
  res@mpCenterLonF   = centerLon         ; choose center lon
  res@mpCenterLatF   = centerLat         ; choose center lon
  getCenteredMap(res,latRadius,lonRadius)
  res@mpGreatCircleLinesOn = True
   
  res@gsnFrame             = False
  res@gsnMaximize          = True
  ; res@gsnSpreadColors      = True
  res@gsnStringFontHeightF = FontHeight
  res@gsnLeftString        = LeftString
  res@gsnCenterString      = CenterString
  res@gsnRightString       = RightString
  
;  res@cnFillOn             = True
;  res@cnLinesOn            = False
;  res@cnInfoLabelOn        = False
;  res@cnFillMode           = "Areafill"
;  res@cnRasterSmoothingOn  = True
;  res@cnRasterSmoothingOn  = False
;  res@cnFillBackgroundColor= "Transparent"
;  res@cnFillColor          = "Transparent"
;  res@cnFillBackgroundColor= "White"
;  res@cnFillColor          = "White"
;  res@cnMissingValFillColor = -1

;  res@sfXArray       = x  ; These are 1D arrays, so a triangular
;  res@sfYArray       = y  ; mesh will be created internally.

  res@tiMainFont           = 4
  res@tiMainFontHeightF    =  FontHeight
  res@tiMainFontThicknessF = 4
  res@tiMainOffsetYF       = 0.02
  res@tiMainString         = TitleString
   
  res@pmLabelBarHeightF        = 0.07
  res@pmLabelBarWidthF         = 0.7
  res@pmLabelBarOrthogonalPosF = 0.05
  
  res@gsnDraw        = False             ; don't draw the plots
  res@gsnFrame       = False             ; don't advance the frame
    
  plot = gsn_csm_map(wks,res)
      
;************************************************
; prepare plot
  edgeVert1 = gridFile->edge_vertices(0,:)-1
  edgeVert2 = gridFile->edge_vertices(1,:)-1
  vertLon = gridFile->vlon * rad2deg
  vertLat = gridFile->vlat * rad2deg
  plotVar = gridFile->$varName$(:)
  res@cnLevelSelectionMode = "ManualLevels"

  
  min_value =  new((/1/),"double")
  max_value =  new((/1/),"double")
  var_value =  new((/1/),"double")
  if( .not. isvar("minValue")) then
    min_value = min(plotVar)
  else
    min_value = minValue
  end if
  if( .not. isvar("maxValue")) then
    max_value = max(plotVar)
  else
    max_value = maxValue
  end if

  edgeDim = dimsizes(edgeVert1)
  noOfEdges = edgeDim(0)
  vertDim = dimsizes(vertLon)
  noOfVertices = vertDim(0)

  colorList=gsn_retrieve_colormap(wks)
  colorDims= dimsizes(colorList)
  numOfColors = colorDims(0)
 ; print("numOfColors:"+numOfColors(0))
  numOfPlotColors = numOfColors-2
  if (max_value .eq. min_value) then
     print("Warning: max_value .eq. min_value")
     max_value = max_value + 0.5
     min_value = min_value - 0.5
  end if
  colorsDistr =  numOfPlotColors/ (max_value - min_value)
;************************************************
; create plot

  resp                  = True                      ; polyline mods desired
  resp@gsLineThicknessF = edgeThickness             ; thickness of lines
;  resp@gsnMaximize      = True

  drawVertLon = new((/2/),"double")
  drawVertLat = new((/2/),"double")
  edgecolor=1
  print("Writing "+noOfEdges(0)+" edges...") 
  do i = 0 , noOfEdges-1
        
    drawVertLon(0) = vertLon(edgeVert1(i))
    drawVertLon(1) = vertLon(edgeVert2(i))
    drawVertLat(0) = vertLat(edgeVert1(i))
    drawVertLat(1) = vertLat(edgeVert2(i))      
    
    vert1in = pointIn(drawVertLon(0),drawVertLat(0),res)
    vert2in = pointIn(drawVertLon(1),drawVertLat(1),res)
    if (vert1in .or. vert2in) then
      pdum = unique_string("pdum")
      var_value=plotVar(i)
      if (var_value .gt. max_value) then
        var_value=max_value
      else if (var_value .lt. min_value) then
        var_value=min_value
      end if 
      end if 
      edgecolor = doubletointeger((var_value - min_value) * colorsDistr) + 2
;      if (edgecolor .lt. 2) then
;        edgecolor=2
;      end if
      resp@gsLineColor      = edgecolor                 ; color of lines
     ; print (plotVar(i))
     ; print (edgecolor)
      plot@$pdum$ = gsn_add_polyline(wks,plot,drawVertLon,drawVertLat,resp)
    end if
  end do
 
;************************************************
;  clean-up  
   delete(edgeVert1)
   delete(edgeVert2)
   delete(vertLon)
   delete(vertLat)
   delete(drawVertLon)
   delete(drawVertLat)
   delete(resp)
;************************************************


;************************************************
; create label
  x = 0.1 ; X position of labelbar
  y = 0.4 ; Y position of labelbar

  lbres = True ; Set up labelbar resources.

  lbres@vpWidthF = 0.7 ; Change width and height
  lbres@vpHeightF = 0.1 ; of labelbar.

  lbres@lbOrientation = "Horizontal"
  lbres@lbPerimOn = False

  lbres@lbFillColors = colorList(2:numOfColors-1,:)
  lbres@lbMonoFillPattern = True ; Use solid fill, not patterns

  lbres@lbLabelAlignment = "InteriorEdges" ; default is "BoxCenters"
  lbres@lbLabelFontHeightF = 0.008

  ; lbres@lbTitleString = plotVar@long_name
  lbres@lbTitleString = ""
  lbres@lbTitleFontHeightF = 0.0
  lbres@lbLabelStride      = numOfColors/5 + 1
  ; lbres@lbLabelAngleF      = 45                 ; angle labels

  levels = new(numOfPlotColors+1,double)
  ; maxMinDiff = doubletofloat(max_value-min_value)
  maxMinDiff = (max_value-min_value)
  do i=0,numOfPlotColors
;   levels(i) = i+1
    levels(i) = (min_value) + (maxMinDiff/int2flt(numOfPlotColors)) * int2flt(i)
  end do
  
  labels = "" + levels
  lbid = gsn_create_labelbar(wks,dimsizes(levels),labels,lbres)
  
;
; Default was that labelbar appeared in center of plot.
;
; Change the zone to 2 so the labelbar is outside the
; plot then center it under plot and move it down a smidge.
;
  amres = True
  amres@amZone           = 1
  amres@amParallelPosF   = 0.5            ; Center labelbar.
  amres@amOrthogonalPosF = 0.0075            ; Move down, away from plot
  
;  annoid = gsn_add_annotation(plot,lbid,False)
  annoid = gsn_add_annotation(plot,lbid,amres)
  
;************************************************
  
  draw(plot)
  frame(wks)
    
  wallClockElapseTime(wcStrt, "create ocean plot", 0)

   print (graphicsFile+"."+graphicsFormat +" is created.")
  if (graphicsFormat .eq. "png" ) then
    system("eog "+graphicsFile+"."+graphicsFormat +" &")
  end if
  if (graphicsFormat .eq. "pdf" ) then
    system("acroread "+graphicsFile+"."+graphicsFormat +" &")
  end if
  if (graphicsFormat .eq. "ps" ) then
;   system("kghostview "+graphicsFile+"."+graphicsFormat +" &")
    system("gv "        +graphicsFile+"."+graphicsFormat +" &")
  end if

end
