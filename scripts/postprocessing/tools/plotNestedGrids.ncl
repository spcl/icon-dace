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
; Plots samples of sphere grids used by ICON
; Leonidas Linardakis (MPI-M 2009-11-05)
;
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
  
;--------------------------------------- 
;  change to the graphicsFormat you prefer
 graphicsFormat="ps"
 myColorMap=(/"white","black","tan","LightGoldenrod1", "LightBlue","Blue","Red", "coral3", "cyan", "darkolivegreen1",\
    "darkgreen"/)
;--------------------------------------- 

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
; checks if point is inside min-max drawing area
; Note: This assumes that the min-max drawing area is correct
;************************************************
function pointIn(inlon:double,inlat:double,res:logical)
begin
  
  tlon=inlon
  tlat=inlat
  if (tlon .gt. 180.) then
     tlon = tlon - 360.0
  end if
;  inPointShift = 2 * abs(res@mpCenterLatF) + 1.0
  inPointShift = 1.0
 
  if (tlon .lt. (res@mpMinLonF - inPointShift)  ) then
     return(False)
  end if
  if (tlon .gt. (res@mpMaxLonF + inPointShift)) then
     return(False)
  end if
  
  if (tlat .lt. (res@mpMinLatF - inPointShift) ) then
     return(False)
  end if
  if (tlat .gt. (res@mpMaxLatF + inPointShift) ) then
     return(False)
  end if
  
  return(True)

end
;************************************************

;************************************************
; set min-max for the whole globe
;************************************************
procedure getWholeGlobe(res:logical)
begin 
  res@mpMinLatF      = -90.0
  res@mpMaxLatF      = 90.0
  res@mpMinLonF      = -180
  res@mpMaxLonF      = 180
end
;************************************************

;************************************************
; set min-max cenetered at Lat=0, res@mpCenterLonF
;************************************************
procedure getCenteredGlobe(res:logical)
begin 
  res@mpCenterLatF   = 0.0              ; choose center lat
  res@mpMinLatF      = -90.0
  res@mpMaxLatF      = 90.0
  res@mpMinLonF      = res@mpCenterLonF - 90.0 
  res@mpMaxLonF      = res@mpCenterLonF + 90.0 
end
;************************************************

;************************************************
; set min-max cenetered at Lat=0, res@mpCenterLonF
;************************************************
procedure getNorthPole(res:logical)
begin 
  res@mpCenterLonF   = 0.0              ; choose center lat
  res@mpCenterLatF   = 90.0              ; choose center lat
  res@mpMinLatF      = 0.0
  res@mpMaxLatF      = 90.0
  res@mpMinLonF      = -180.0
  res@mpMaxLonF      =  180.0
end
;************************************************

;************************************************
; the main drawing routine 
;************************************************
function iconDrawGrid(fileName:string, wks:graphic, plot:graphic, edgecolor:string,edgethick:float,res:logical)
;************************************************

begin 

; read in grid netCDF file
  a = addfile(fileName,"r")
  print("Processing "+fileName(0)+"...")

; fill in edges, verts
  edgeVert1 = a->edge_vertices(0,:)-1
  edgeVert2 = a->edge_vertices(1,:)-1
  vertLon = a->vlon * rad2Dgr
  vertLat = a->vlat * rad2Dgr
    
  edgeDim = dimsizes(edgeVert1)
  noOfEdges = edgeDim(0)
  vertDim = dimsizes(vertLon)
  noOfVertices = vertDim(0)
 
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
    
;    vert1in = pointIn(drawVertLon(0),drawVertLat(0),res)
;    vert2in = pointIn(drawVertLon(1),drawVertLat(1),res)
;    if (vert1in .or. vert2in) then
      pdum = unique_string("pdum")
      plot@$pdum$ = gsn_add_polyline(wks,plot,drawVertLon,drawVertLat,resp)
;    end if
  end do
        
 
;************************************************
;  clean-up  
   delete(drawVertLon)
   delete(drawVertLat)
   delete(resp)
;************************************************
  
  
    return(noOfEdges)
end
; iconDrawGrid ends here
;************************************************
 

 
;************************************************
; main program
; sets the projection and plot parameters,
; and just calls one or more of the above procedures
;************************************************
; main program
begin
  wcStrt     = systemfunc("date")

  res                = True              ; plot mods desired
  res@mpProjection = "Orthographic"
;  res@mpProjection = "Satellite"
;  res@mpProjection = "LambertEqualArea"
;  res@mpProjection = "Stereographic"
;  res@gsnCenterString = res@mpProjection+" Projection"
 
  res@mpLimitMode    = "LatLon"          ; use lat/lon coordinates to limit area
  res@mpGreatCircleLinesOn = True

  res@gsnMaximize = True
  res@gsnPaperOrientation  = "portrait"  ; force portrait
  res@gsnDraw        = False             ; don't draw the plots
  res@gsnFrame       = False             ; don't advance the frame
  res@mpPerimOn      = False             ; no box around map
   
   res@mpDataBaseVersion    = "MediumRes"   ; Medium resolution database
 
   res@mpFillOn              =  True           ; turn on map fill
;  res@mpOutlineDrawOrder    = "PostDraw"       ; draw continental outline last
;  res@mpAreaMaskingOn        = True                 ; draw only parts of map 
;   res@mpOutlineDrawOrder          = "PostDraw"
;   res@mpFillDrawOrder             = "PreDraw"
;   res@mpFillAreaSpecifiers   = (/"Water","Land"/)   ; what parts to fill
;  res@mpOutlineOn            = True
;  res@mpGridAndLimbOn        = True
;  res@mpGridLineDashPattern  = 2
;  res@mpGridMaskMode         = "MaskLand"

  res@mpLandFillColor        = "LightGoldenrod1"
  res@mpOceanFillColor       = "LightBlue"
  res@mpInlandWaterFillColor = "LightBlue"
      
;************************************************
  cntLon = 0.0
  cntLat = 0.0
  if( isvar("lon")) then
    cntLon=lon
  end if
  if( isvar("lat")) then
    cntLat=lat
    print(cntLat)
  end if
;************************************************
  if( .not. isvar("plotName")) then
     plotName="plot"
  end if
  wks = gsn_open_wks(graphicsFormat,plotName)             ; open a graphicsFormat file
  gsn_define_colormap(wks,myColorMap)

  if (cntLat .eq. 90) then
    getNorthPole(res)
  else if (cntLat .eq. 0) then
    res@mpCenterLonF   = cntLon
    getCenteredGlobe(res)
  else
    res@mpCenterLonF   = cntLon         ; choose center lon
    res@mpCenterLatF   = cntLat         ; choose center lon
;    res@mpMinLatF      = minLat
;    res@mpMaxLatF      = maxLat
;    res@mpMinLonF      = minLon
;    res@mpMaxLonF      = maxLon
  end if
  end if

; plot earth  
  plot = gsn_csm_map(wks,res) 
 
; plot grid
  if( .not. isvar("gridFolder")) then
     gridFolder=""
  end if
  
  if( isvar("grid1")) then
    edges=iconDrawGrid(gridFolder+grid1,wks,plot,"blue",1.0,res)
  end if
  if( isvar("grid2")) then
    edges=iconDrawGrid(gridFolder+grid2,wks,plot,"green",1.0,res)
  end if
  if( isvar("grid3")) then
    edges=iconDrawGrid(gridFolder+grid3,wks,plot,"cyan",0.8,res)
  end if
  if( isvar("grid4")) then
    edges=iconDrawGrid(gridFolder+grid4,wks,plot,"coral3",0.8,res)
  end if
  
  draw(plot)
  frame(wks)
  delete(plot)
  delete(wks)
 
  wallClockElapseTime(wcStrt, "create ocean plot", 0)
  print(plotName(0)+"."+graphicsFormat(0)+" is created.")
  
  if (graphicsFormat .eq. "png" ) then
    system("eog "+plotName+"."+graphicsFormat +" &")
  end if
  if (graphicsFormat .eq. "pdf" ) then
    system("acroread "+plotName+"."+graphicsFormat +" &")
  end if
  if (graphicsFormat .eq. "ps" ) then
    system("kghostview "+plotName+"."+graphicsFormat +" &")
  end if
        
 
end
;************************************************
