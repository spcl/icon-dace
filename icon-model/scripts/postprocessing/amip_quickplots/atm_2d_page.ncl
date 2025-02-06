; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


undef("define_draw_strings")
 procedure define_draw_strings (Code,cArea,comment,wks,plot)
 
;************************************************
;  Define min max mean
;  Draw text on plot using NDC coordinates.
;************************************************
local txres

begin

;***compute mean
  meanV=dim_avg(Code*cArea)/dim_avg(cArea)

  minV=min(Code)
  maxV=max(Code)

  statV="min="+minV+"   mean="+meanV +"  max="+maxV

  pltdate = systemfunc("date +%d-%m-%Y")
;;;  print("actual_date: "+pltdate)

;--KMF create text object and add it to the plot
  txres                  =  True                 ; text mods desired
  txres@txFontHeightF    =  0.014                 ;-- smaller font
  txid = gsn_create_text(wks, statV, txres)
  
  amres                  =  True
  amres@amParallelPosF   =  0.0                  ;-- center position
  amres@amOrthogonalPosF =  0.7                  ;-- move downward
  amres@amSide           = "Bottom"
  annoid = gsn_add_annotation(plot, txid, amres)

end
;***end procedure 

;---------------------------------------
;-- function attach_labelbar
;---------------------------------------
undef("attach_labelbar")
function attach_labelbar(wks,map,labels,colors)
local lbres, nlevels, amres
begin
  nlevels = dimsizes(labels)                     ; number of labels

  lbres                      =  True
  lbres@lbPerimOn            =  False            ; no label bar box
  lbres@lbOrientation        = "Horizontal"      ; orientation
  lbres@vpWidthF             =  0.7              ; width of labelbar
  lbres@vpHeightF            =  0.10             ; height of labelbar
  lbres@lbLabelFontHeightF   =  0.015            ; label font height
  lbres@lbLabelAlignment     = "InteriorEdges"   ; where to label
  lbres@lbMonoFillPattern    =  True             ; fill sold
  lbres@lbFillColors         =  colors           ; use colors
  lbid = gsn_create_labelbar (wks,nlevels+1,labels,lbres)

  amres                      =  True
  amres@amJust               = "TopCenter"       ; annotation alignment
  amres@amOrthogonalPosF     =  0.7              ; move annotation downward
;  amres@gsnMaximize          =  True             ; be sure to remaxmize output
  map@annoid                 =  gsn_add_annotation(map,lbid,amres)

  return(map)
end
;***end procedure 
undef("fast_cell")
 function  fast_cell  (wks,var,x,y,cArea,vlon,vlat,cnres_in,missColor,\
                        mainTitle,subTitle,comment )
;************************************************
;  drawing triangle for ICON FASTER!!
;************************************************
local plot, cnres, mpres, pres, res

begin

  nv   = dimsizes(vlon(0,:))                      ; no of points in polygon

  lonMin       = -180                             ; longitude minimum
  lonMax       =  180                             ; longitude maximum
  latMin       =  -90                             ; latitude minimum
  latMax       =   90                             ; latitude maximum
  mapCenter    =    0                             ; center of map

;  print("")
;  print("Data longitude min/max: " + min(vlon) + "   " + max(vlon)\
;      + "   cell points = " + nv)
;  print("")
;  print("Plot area: "+lonMin+","+lonMax+" "+latMin+","+latMax)
;  print("")

;***set contour resources
  res=cnres_in
  res@gsnDraw              =  False        ;-- don't draw the plot
  res@gsnFrame             =  False        ;-- don't advance the frame
  res@cnLinesOn            =  False        ;-- don't draw contour lines
  res@cnInfoLabelOn        =  False        ;-- switch off contour info label
  res@cnFillOn             =  True         ;-- contour fill on

  res@cnFillMode           = "CellFill"    ;-- set fill mode
  res@sfXArray             =  x            ;-- transform x to mesh scalar field
  res@sfYArray             =  y            ;-- transform y to mesh scalar field
  res@sfXCellBounds        =  vlon         ;-- needed if cnFillMode= "CellFill"
  res@sfYCellBounds        =  vlat         ;-- needed if cnFillMod = "CellFill"
  
  res@lbLabelBarOn         =  False        ;-- don't draw a labelbar yet

  res@mpFillOn             =  True         ;?-- fill map grey
  res@mpOutlineOn          =  True 
  res@mpLandFillColor      =  18

  res@mpMinLonF            =  lonMin         ;-- sub-region minimum longitude
  res@mpMaxLonF            =  lonMax         ;-- sub-region maximum longitude
  res@mpMinLatF            =  latMin         ;-- sub-region minimum latitude
  res@mpMaxLatF            =  latMax         ;-- sub-region maximum latitude
  res@mpGreatCircleLinesOn =  False          ;-- important: v6.2.0 False !!
  res@gsnSpreadColors      =  False  
    
  res@tiMainString         =  mainTitle
  res@gsnCenterString      =  subTitle
  res@gsnLeftString        = " " 
  res@gsnRightString       = " " 


;***create the contour plot, but don't draw it. 
;***We need it to get the colors of the data values
  contour_object = gsn_csm_contour_map(wks,var,res)


;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (var,cArea,comment,wks,contour_object)

;***to plot the filled polygons we need the data levels and their colors
  getvalues contour_object@contour
      "cnLevels"     :  levels      ;-- retrieve the levels from plot@contour
      "cnFillColors" :  colors      ;-- retrieve the colors from plot@contour
  end getvalues
  nlevels = dimsizes(levels)        ;-- number of levels

;***clear plot, but keep all the information
  contour_object = setColorContourClear(contour_object,min(var),max(var)) 
print("min: "+min(var)+" max: "+max(var))

;***attach a labelbar to the plot
    contour_object = attach_labelbar(wks,contour_object,levels+"",colors)

;***create color array for triangles
  ntri     = dimsizes(y)          ;-- number of triangles
  gscolors = new(ntri,string)
  gscolors = "gray"               ;-- gray for missing colors

;***set resources for the triangles (polygons)
  pres                         =  True
  pres@gsEdgesOn               =  False            ;-- turn on edges
  pres@gsFillIndex             =  0                ;-- solid fill

  ntri_calc = 0      ;-- number of triangles

;*** set color for data less than given minimum value var_min
  vlow = ind(var .lt. levels(0))   ;-- get the indices of values less levels(0)
  if (.not. all(ismissing(vlow))) then
     gscolors(vlow) = colors(0)       ;-- choose color
     ntri_calc = dimsizes(vlow)       ;-- number of triangles
  end if
;***set colors for all cells in between var_min and var_max
  do i = 1, dimsizes(levels) - 1
;    get the indices of 'middle' values
     vind := ind(var .ge. levels(i-1) .and. var .lt. levels(i))  
     if (.not. all(ismissing(vind))) then
        gscolors(vind) = colors(i)                  ;-- choose the colors
        ntri_calc = ntri_calc + dimsizes(vind)      ;-- number of triangles
     end if
  end do

;***set color for data greater than given maximum var_max
  nc=dimsizes(colors)-1              ; get the number of colors minus one
  nl=dimsizes(levels)-1              ; get the number of levels minues one
  vhgh := ind(var .gt. levels(nl))   ; get indices of values greater levels(nl)
  if (.not. all(ismissing(vhgh))) then
     gscolors(vhgh) = colors(nc)                  ; choose color
     ntri_calc = ntri_calc + dimsizes(vhgh)       ; number of triangles
  end if

;--NEW set color for missing values
  ismis = ind(ismissing(var))
  if(.not. all(ismissing(ismis))) then
     gscolors(ismis) = missColor                     ; set missing value color
  end if
;--NEW

;  print("--> triangles calculated:     "+ ntri_calc)

;***attach all the triangles using the list of colors
  pres@gsColors        = gscolors       ; set colors for polygons
  pres@gsSegments      = ispan(0,dimsizes(var) * 3,3)  ; assign segments array
;***draw all triangles
  polygon = gsn_add_polygon(wks,contour_object,ndtooned(vlon),ndtooned(vlat),pres)  

;-- define attribute polygon to hold the polygon ids
  contour_object@polygon = polygon

  return(contour_object)

end
;***end function fast_cell
;--------------

undef("plot_cont")
 function plot_cont  ( wks,\
                        var,x,y,cArea,res_in,\
                        mainTitle,subTitle,comment )
;************************************************
;  contour plots for ICON
;************************************************
local plot, res

begin

;***resource list for color contour levels
    res=res_in
;    res                      = True
    res@gsnFrame             = False             ; don't advance frame
    res@gsnDraw              = False             ; don't draw plot yet
    res@cnFillOn       = True
    res@cnLinesOn      = False
    res@cnInfoLabelOn  = False

    res@tiMainString         = mainTitle
    res@gsnCenterString      = subTitle
    res@gsnLeftString        = " " 
    res@gsnRightString       = " " 

    res@mpCenterLonF         = 0         ; center of map
    res@gsnSpreadColors      = False     ; use full range of colors

    res@pmLabelBarOrthogonalPosF = 0.2           ; move farther to plot

    res@sfXArray = x      ; These are 1D arrays, so a triangular mesh
    res@sfYArray = y      ; will be used to create the contours.

    contour_object = gsn_csm_contour_map_ce(wks,var,res)

;***define statistical Values and draw Values, comment, pltdate
    define_draw_strings (var,cArea,comment,wks,contour_object)

  return(contour_object)

end
;***end function 

;************************************************
; main
;************************************************
begin

  values    = asciiread("var.txt",7,"string")
  name = values(0)
  typ = values(1)
  run = values(2)
  meantime = values(3)
  comment = values(4)
  workdir = values(5)        ; plot & input directory
  eratime = values(6)

;  print(eratime)
;  print(workdir)
;  print(pltdate)



;***define plot and open
   pltType = "pdf"                        ; x11, eps, ps, pdf, png
   pltPath = workdir
   pltName = pltPath+"/"+run+"_atm_2d_"+typ+"_map"        ; plot name
   wks     = gsn_open_wks(pltType, pltName) 
;***plot array for panel
   listplots = new(12,graphic)
   noPlot=-1  ;number of listplots
;***define colors
    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    i = NhlNewColor(wks,0.8,0.8,0.8)              ; add gray to color map 
    i = NhlNewColor(wks,0.6,0.6,0.6)               
    i = NhlNewColor(wks,0.4,0.4,0.4)  
    missColor = NhlNewColor(wks,0.5,0.5,0.5)     ; add gray to color map

;************************************************
;  Read grid information
;  define the x-, y-values and the polygon points
;************************************************
  GridInfoFileName = systemfunc("echo $GrdInfoFile")

  GridInfoFile= addfile( GridInfoFileName, "r" )

  rad2deg = 45./atan(1.)            ; radians to degrees

  x = GridInfoFile->clon *rad2deg   ; cell center, lon
  y = GridInfoFile->clat *rad2deg   ; cell center, lat
  cArea = GridInfoFile->cell_area   ; area of grid cell

; note: clon and clat are longitude and latitude of triangle centers.
;       Locations of the cell corners are given by 
;       clon_vertices and clat_vertices in the nc file.


  x!0     = "lon"
  x@units = "degrees_east"
  y!0     = "lat"
  y@units = "degrees_north"

  vlon = GridInfoFile->clon_vertices * rad2deg    ; cell longitude vertices
  vlon = where(vlon.lt.0, vlon + 360, vlon)       ; longitude: 0-360
  vlat = GridInfoFile->clat_vertices * rad2deg    ; cell lattitude vertices



;************************************************
;   Zonal wind stress CMOR: tauu only sea
;************************************************

;***setings
    Cvar="tauu"
    fili    = workdir+"/Ubusy_"+Cvar+".nc"

    mainTitle = "Zonal wind stress [mN/m~S~2~N~] "+typ
    subTitle  = run+meantime

    if (default.eq.1.or.cell.eq.1) then   ; set key and appendix for 
       fillKey="cell"                     ; contour or trigangels
    else
       fillKey="cont"                      
    end if

if (isfilepresent(fili)) then

;***read code tauu
    File       = addfile (fili , "r") 
    var     = File->$Cvar$(time|0,ncells|:)  ; dims: (time,ncells)
    print("ploted: "+mainTitle)

;***number of listplots
    noPlot=noPlot+1


;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/-200,-100,-50,-20,0,20,50,100,200/)
    res@cnFillColors         = (/2,4,6,8,9,10,11,12,14,15,16/)
    res@cnMissingValFillColor = missColor

 if (fillKey.eq."cont") then
     plot0 = plot_cont ( wks,var,x,y,cArea,res,\
                 mainTitle,subTitle,comment )
 else 
    plot0 = fast_cell (  wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                mainTitle,subTitle,comment)
 end if
 listplots(noPlot)=plot0

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if

;************************************************
;   Zonal wind stress CMOR: tauu  model - ERAinterim
;************************************************

;***setings
    Cvar="tauu"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"

    mainTitle = "Zonal wind stress [mN/m~S~2~N~] "+typ
    subTitle  =  run+meantime+" - ERAinterim"+eratime

    if (default.eq.1.or.cell.eq.0) then   ; set key and appendix for 
       fillKey="cont"                     ; contour or trigangels
    else
       fillKey="cell"                      
    end if

if (isfilepresent(fili)) then

;***read code tauu
    File       = addfile (fili , "r") 
    var     = File->$Cvar$(time|0,ncells|:)  ; dims: (time,ncells)
    print("ploted: "+mainTitle)

;***number of listplots
    noPlot=noPlot+1

    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/-75,-60,-45,-30,-15,15,30,45,60,75/)
    res@cnFillColors         = (/2,4,6,8,9,0,10,11,12,14,15,16/)
    res@cnMissingValFillColor = missColor


;***create plot
if (fillKey.eq."cont") then
    plot1 = plot_cont ( wks,var,x,y,cArea,res,\
                mainTitle,subTitle,comment )
else 
    plot1 = fast_cell ( wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                            mainTitle,subTitle,comment)
end if
    listplots(noPlot)=plot1

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if


;************************************************
;   Total precipitation code4  CMOR: pr
;************************************************
    Cvar="pr"
    fili    = workdir+"/Ubusy_"+Cvar+".nc"

    mainTitle = "Total precipitation [mm/d] "+typ
    subTitle  = run+meantime

    if (default.eq.1.or.cell.eq.1) then   ; set key and appendix for 
       fillKey="cell"                     ; contour or trigangels
    else
       fillKey="cont"                      
    end if          

if (isfilepresent(fili)) then

;***read code pr
    File      = addfile (fili , "r")  
    var     = File->$Cvar$(time|0,ncells|:)        
;   printVarSummary(var)
    print("ploted: "+mainTitle)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/1,2,4,8,12/)
    res@cnFillColors         = (/0,18,9,6,3,17/)
    res@cnMissingValFillColor = missColor

if (fillKey.eq."cont") then
    plot2 = plot_cont ( wks,var,x,y,cArea,res,\
                mainTitle,subTitle,comment )
else 
    plot2 = fast_cell ( wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                mainTitle,subTitle,comment )
end if
    listplots(noPlot)=plot2

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if
;************************************************
;   Total cloud cover CMOR: clt
;************************************************

;***setings
    Cvar="clt"
    fili     =  workdir+"/Ubusy_"+Cvar+".nc"
  
    mainTitle = "Total cloud cover [%] "+typ
    subTitle  = run+meantime

    if (default.eq.1.or.cell.eq.1) then   ; set key and appendix for 
       fillKey="cell"                     ; contour or trigangels
    else
       fillKey="cont"                      
    end if
                
if (isfilepresent(fili)) then
;***read code clt
    File    = addfile (fili , "r")  
    var     = File->$Cvar$(time|0,ncells|:)   ; dims: (time,ncells)    
;   printVarSummary(var)
    print("ploted: "+mainTitle)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/10,30,50,70,90/)
    res@cnFillColors         = (/3,6,9,18,19,20/)
    res@cnMissingValFillColor = missColor


if (fillKey.eq."cont") then
    plot3 = plot_cont ( wks,var,x,y,cArea,res,\
                mainTitle,subTitle,comment )
else 
    plot3 = fast_cell ( wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                mainTitle,subTitle,comment )
end if
    listplots(noPlot)=plot3

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if

;************************************************
;   Column water vapour CMOR: prw 
;************************************************

;***setings
    Cvar="prw"
    fili    = workdir+"/Ubusy_"+Cvar+".nc"

    mainTitle = "Column water vapour [kg/m~S~2~N~] "+typ
    subTitle  = run+meantime

    if (default.eq.1.or.cell.eq.1) then   ; set key and appendix for 
       fillKey="cell"                     ; contour or trigangels
    else
       fillKey="cont"                      
    end if

if (isfilepresent(fili)) then

    File       = addfile (fili , "r")   
    var     = File->$Cvar$(time|0,ncells|:)         
;   printVarSummary(var)
    print("ploted: "+mainTitle)

;***number of listplots
   noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/10,20,30,40,50/)
    res@cnFillColors         = (/0,18,9,6,3,17/)
    res@cnMissingValFillColor = missColor

if (fillKey.eq."cont") then
    plot4 = plot_cont ( wks,var,x,y,cArea,res,\
                mainTitle,subTitle,comment )
else 
    plot4 = fast_cell ( wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                mainTitle,subTitle,comment )
end if
    listplots(noPlot)=plot4

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if

;************************************************
;   Column water vapour CMOR: prw (model-ERA)
;************************************************

;***setings
    Cvar="prw"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"

    mainTitle = "Column water vapour [kg/m~S~2~N~] "+typ
    subTitle  = run+meantime+" - ERAinterim"+eratime

    if (default.eq.1.or.cell.eq.0) then   ; set key and appendix for 
       fillKey="cont"                     ; contour or trigangels
    else
       fillKey="cell"                      
    end if

if (isfilepresent(fili)) then

    File       = addfile (fili , "r") 
    var     = File->$Cvar$(time|0,ncells|:)         
    print("ploted: "+mainTitle)

;***number of listplots
   noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-10,-5,-2,-1,1,2,5,10/)
    res@cnFillColors         = (/4,6,7,9,0,11,12,14,15,16/)
    res@cnMissingValFillColor = missColor

if (fillKey.eq."cont") then
    plot5 = plot_cont ( wks,var,x,y,cArea,res,\
                mainTitle,subTitle,comment )
else 
    plot5 = fast_cell ( wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                mainTitle,subTitle,comment )
end if
    listplots(noPlot)=plot5

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if

;************************************************
;   Sea Level Pressure code151  CMOR: psl
;************************************************

    Cvar="psl"
    fili     =  workdir+"/Ubusy_"+Cvar+".nc"

    mainTitle = "Sea Level Pressure (hPa-1000) "+typ
    subTitle  = run+meantime

    if (default.eq.1.or.cell.eq.1) then   ; set key and appendix for 
       fillKey="cell"                     ; contour or trigangels
    else
       fillKey="cont"                      
    end if

if (isfilepresent(fili)) then

;***read code 151
    File       = addfile (fili , "r") 
    var     = File->$Cvar$(time|0,ncells|:)   ; dims: (time,ncells)    
    var     = var - 1000.
    print("ploted: "+mainTitle)

;***number of listplots
   noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/-25,-20,-15,-10,-5,0,5,10,15,20,25/)
    res@cnFillColors         = (/2,4,6,7,8,9,10,11,12,14,15,16/)
    res@cnMissingValFillColor = missColor

if (fillKey.eq."cont") then
    plot6 = plot_cont ( wks,var,x,y,cArea,res,\
                mainTitle,subTitle,comment )
else 
    plot6 = fast_cell ( wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                mainTitle,subTitle,comment )
end if
    listplots(noPlot)=plot6


    delete (res@cnLevels)
    delete (res@cnFillColors)
end if
;************************************************
;   Sea Level Pressure code151 (model-ERAin) CMOR: psl
;************************************************

    Cvar="psl"
    fili     =  workdir+"/Ubusy_"+Cvar+"_erain.nc"

    mainTitle = "Sea Level Pressure (hPa-1000) "+typ
    subTitle  = run+meantime+" - ERAinterim "+eratime

    if (default.eq.1.or.cell.eq.0) then   ; set key and appendix for 
       fillKey="cont"                     ; contour or trigangels
    else
       fillKey="cell"                      
    end if

if (isfilepresent(fili)) then

;***read code 151
    File       = addfile (fili , "r") 
    var     = File->$Cvar$(time|0,ncells|:)   ; dims: (time,ncells)    
    var     = var
    print("ploted: "+mainTitle)

;***number of listplots
   noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/-10,-5,-2,-1,1,2,5,10/)
    res@cnFillColors         = (/2,4,6,9,0,11,12,14,15,16/)

    res@cnMissingValFillColor = missColor

if (fillKey.eq."cont") then
    plot7 = plot_cont ( wks,var,x,y,cArea,res,\
                mainTitle,subTitle,comment )
else 
    plot7 = fast_cell ( wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                mainTitle,subTitle,comment )
end if
    listplots(noPlot)=plot7
print ("plot7: "+ noPlot)


    delete (res@cnLevels)
    delete (res@cnFillColors)
end if


;************************************************
;   2 m Temperature code167  Cmor: tas
;************************************************

    Cvar="tas"
    fili     =  workdir+"/Ubusy_"+Cvar+".nc"
    mainTitle = "2 m Temperature [C] "+typ
    subTitle  = run+meantime

    if (default.eq.1.or.cell.eq.1) then   ; set key and appendix for 
       fillKey="cell"                     ; contour or trigangels
    else
       fillKey="cont"                      
    end if


if (isfilepresent(fili)) then

;***read code 167
    File       = addfile (fili , "r") 
    var     = File->$Cvar$(time|0,ncells|:)   ; dims: (time,ncells)    
;   printVarSummary(var)
    print("ploted: "+mainTitle)

;***number of listplots
   noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30/)
    res@cnFillColors         = (/2,4,5,6,7,8,9,10,11,12,14,15,16,17/)
    res@cnMissingValFillColor = missColor

if (fillKey.eq."cont") then
    plot8 = plot_cont ( wks,var,x,y,cArea,res,\
                mainTitle,subTitle,comment )
else 
    plot8 = fast_cell ( wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                mainTitle,subTitle,comment )
end if
    listplots(noPlot)=plot8
print ("plot8: "+ noPlot)

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if
;************************************************
;   2 m Temperature code167 Cmor: tas model-era
;************************************************

    Cvar="tas"
    fili     =  workdir+"/Ubusy_"+Cvar+"_erain.nc"
    mainTitle = "2 m Temperature [C] "+typ
    subTitle  = run+meantime+" - ERAinterim "+eratime


    if (default.eq.1.or.cell.eq.0) then   ; set key and appendix for 
       fillKey="cont"                     ; contour or trigangels
    else
       fillKey="cell"                      
    end if

if (isfilepresent(fili)) then

;***read code 167
    File       = addfile (fili , "r") 
    var     = File->$Cvar$(time|0,ncells|:)   ; dims: (time,ncells)    
;   printVarSummary(var)
    print("ploted: "+mainTitle)

;***number of listplots
   noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-5,-3,-2,-1,-.5,.5,1,2,3,5/)
    res@cnFillColors         = (/2,4,6,8,9,0,10,11,12,14,15,16/)
    res@cnMissingValFillColor = missColor

if (fillKey.eq."cont") then
    plot9 = plot_cont ( wks,var,x,y,cArea,res,\
                mainTitle,subTitle,comment )
else 
    plot9 = fast_cell ( wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                mainTitle,subTitle,comment )
end if
    listplots(noPlot)=plot9
print ("plot9: "+ noPlot)

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if
;************************************************
;   Surface temperature code169 CMOR: ts
;************************************************

;***setings
    Cvar="ts"
    fili     =  workdir+"/Ubusy_"+Cvar+".nc"
    
    mainTitle = "Surface temperature [C] "+typ
    subTitle  = run+meantime

    if (default.eq.1.or.cell.eq.1) then   ; set key and appendix for 
       fillKey="cell"                     ; contour or trigangels
    else
       fillKey="cont"                      
    end if
                
if (isfilepresent(fili)) then
;***read 
    File    = addfile (fili , "r")  
    var     = File->$Cvar$(time|0,ncells|:)   ; dims: (time,ncells)    
   
print(fili+" fili  "+min(var)+"   "+max(var))
    print("ploted: "+mainTitle)

;***open plot
   noPlot=noPlot+1

             
;***number of listplots
    res                       = True             ; plot mods desired
    res@cnLevelSelectionMode  = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels              = (/-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30/)
    res@cnFillColors          = (/2,4,5,6,7,8,9,10,11,12,14,15,16,17/)
    res@cnMissingValFillColor = missColor

if (fillKey.eq."cont") then
    plot10 = plot_cont ( wks,var,x,y,cArea,res,\
                mainTitle,subTitle,comment )
else 
    plot10 = fast_cell ( wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                mainTitle,subTitle,comment )
end if
    listplots(noPlot)=plot10
print ("plot10: "+ noPlot)

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if

;************************************************
;   Liquid water + ice content code231_150 
;************************************************
;***setings
    Cvar="clwvi"
    fili     =  workdir+"/Ubusy_"+Cvar+".nc"
    
    mainTitle = "Liquid water + ice content [g/m~S~2~N~] "+typ
    subTitle  = run+meantime

    if (default.eq.1.or.cell.eq.1) then   ; set key and appendix for 
       fillKey="cell"                     ; contour or trigangels
    else
       fillKey="cont"                      
    end if
                
if (isfilepresent(fili)) then
;***read 
    File    = addfile (fili , "r")  
    var     = File->$Cvar$(time|0,ncells|:)   ; dims: (time,ncells)    
;   printVarSummary(var)
    print("ploted: "+mainTitle)

;***open plot
    noPlot=noPlot+1
             
;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/10,50,100,200,300/)
    res@cnFillColors         = (/3,6,8,18,19,20/)
    res@cnMissingValFillColor = missColor

if (fillKey.eq."cont") then
    plot11 = plot_cont ( wks,var,x,y,cArea,res,\
                mainTitle,subTitle,comment )
else 
    plot11 = fast_cell ( wks,var,x,y,cArea,vlon,vlat,res,missColor,\
                mainTitle,subTitle,comment )
end if
    listplots(noPlot)=plot11
print ("plot11: "+ noPlot)

    delete(res@cnLevels)
    delete (res@cnFillColors)
end if



;*** PANNEL
print("pannel:  "+noPlot)
;--KMF prepare plots for panel   
   pltres          =  True
   


   gsn_panel(wks, (/listplots(8),listplots(9),\
                    listplots(6),listplots(7),\
                    listplots(0),listplots(1)/), (/3,2/), pltres)   
   gsn_panel(wks, (/listplots(4),listplots(5),\
                    listplots(11),listplots(2),\
                    listplots(3),listplots(10)/), (/3,2/), pltres)   

exit
end
