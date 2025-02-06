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
; plots the grid and the sea-land-mask (optionally the bathymetry)
;
; usage:
; ncl plotOceanGird.ncl 'grid=<grid file name>' ['folder=<grid folder>' 'var=<variable>]
;
;*************************************************
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;---------------------------------------------------------------
 
;---------------------------------------
; some parameters
 rad2Dgr = 57.295779513082323110286093015487255497
;---------------------------------------

;--------------------------------------- 
; parameters for interpolating along a geodesic 
 maxInterpPoints = 150
 noInterpPoints=2
 minLineDist = 5.0d0 ; line length ~5km near the poles
 interSize = 2.0d0
;--------------------------------------- 
;************************************************
; set min-max cenetered at Lat=0, res@mpCenterLonF
;************************************************
procedure getCenteredMap(res:logical, latRadius:double, lonRadius:double)
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


;************************************************
; the main drawing routine 
;************************************************
function iconDrawGrid(gridfile:file, wks:graphic, plot:graphic, edgecolor:string,edgethick:float,res:logical)
;************************************************

begin 


  edgeVert1 = gridfile->edge_vertices(0,:)-1
  edgeVert2 = gridfile->edge_vertices(1,:)-1
  vertLon = gridfile->vlon * rad2Dgr
  vertLat = gridfile->vlat * rad2Dgr
    
  edgeDim = dimsizes(edgeVert1)
  noOfEdges = edgeDim(0)
  vertDim = dimsizes(vertLon)
  noOfVertices = vertDim(0)
 
  ; get edge size in km
  edgeDist = gc_latlon(vertLat(edgeVert1(0)),vertLon(edgeVert1(0)), vertLat(edgeVert2(0)), vertLon(edgeVert2(0)), 0, 4)
  
  ; calculate geodesic interpolation size 
  interSize = edgeDist / minLineDist
   if (interSize .gt. maxInterpPoints-2) then
     interSize = maxInterpPoints-2
  end if

;************************************************
; create plot

  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = edgecolor                 ; color of lines
  resp@gsLineThicknessF = edgethick                 ; thickness of lines
;  resp@gsnMaximize      = True

  drawVertLon = new((/2/),"double")
  drawVertLat = new((/2/),"double")
  
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
  
  
    return(noOfEdges)
end
; iconDrawGrid ends here
;************************************************


 begin
  wcStrt     = systemfunc("date")
  
  gridFolder="../../grids/"
  
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
     statsFile=gridFile
  end if
  
  if( isvar("var")) then
     varName = var
  else
     varName="cell_sea_land_mask"
  end if

  
  GraphicFile  = "test"  
  graphicsFormat="pdf"
  colormap  = "BlueDarkRed18"
  rad2deg = 45./atan(1.)   ; radians to degrees
  
  plotVar = statsFile->$varName$(:)   ; dims: (cell)
   
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

; Labeling
  
  LeftString   = ""
  RightString  = ""
  CenterString = gridFileName
  TitleString  = plotVar@long_name
  FontHeight   = 0.02
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
; coordinates of global coarse grid
  x = gridFile->clon *rad2deg   ; cell center, lon
  y = gridFile->clat *rad2deg   ; cell center, lat

;---------------------------------------------------------------
; make plot
;---------------------------------------------------------------

  wks = gsn_open_wks(graphicsFormat,GraphicFile)
  gsn_define_colormap(wks,colormap)

; Set up resources for contour/map plot.
  res                = True

  res@mpProjection         = "Orthographic"
  res@mpPerimOn            = False             ; turn off box around plot
  res@mpFillOn             = False
  centerLon = 90
  centerLat = 0
  latRadius = 90
  lonRadius = 90    
  res@mpCenterLonF   = centerLon         ; choose center lon
  res@mpCenterLatF   = centerLat         ; choose center lon
  getCenteredMap(res,latRadius,lonRadius)
  res@mpGreatCircleLinesOn = True
   
  res@gsnFrame             = False
  res@gsnMaximize          = True
  res@gsnSpreadColors      = True
  res@gsnStringFontHeightF = FontHeight
  res@gsnLeftString        = LeftString
  res@gsnCenterString      = CenterString
  res@gsnRightString       = RightString
  res@cnFillOn             = True
  res@cnLinesOn            = False
  res@cnInfoLabelOn        = False
  res@cnFillMode           = "Areafill"
  res@cnRasterSmoothingOn  = True
  res@cnRasterSmoothingOn  = False

  res@sfXArray       = x  ; These are 1D arrays, so a triangular
  res@sfYArray       = y  ; mesh will be created internally.

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
      
  plot = gsn_csm_contour_map(wks,plotVar,res)
     
  edges=iconDrawGrid(gridFile,wks,plot,"green",1.0,res)

  draw(plot)
  frame(wks)
    
  wallClockElapseTime(wcStrt, "create ocean plot", 0)
 
  system("acroread "+GraphicFile+"."+graphicsFormat +" &")
end
