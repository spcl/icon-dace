; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

; ICON vertical cross section
;
; Pilar Ripodas, DWD, Nov 2010
;
;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*************************************************


;*************************************************
; Main code
;************************************************* 

begin


  
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

  rad2deg = 45./atan(1.)   ; radians to degrees
  extra_level = 0      ;it is 1 only for var="W" and Model="ICONAM"
  vlev=31 ;"mlev"          ; it is ilev in case var="W"
  iu=1    ; if +1, then lon between 0,360  ; if  -1, then  lon between -180 and 180
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; get the variables from the driver
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  Model      = getenv("Model")
  File_Int   = getenv("file_int")
  var        = getenv("var")
  Rightlat   = getenv("rightlat")
  Leftlat    = getenv("leftlat")
  Rightlon   = getenv("rightlon")
  Leftlon    = getenv("leftlon")
  FileDir    = getenv("model_data_path")
  PlotDir    = getenv("plot_file_path")
  PlotFileN  = getenv("plot_file")
  PlotFormat = getenv("plot_file_format")
  IMan       = getenv("i_manual_levels")

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; read NAMELIST settings from directory
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Turbulence = systemfunc("grep inwp_turb "+FileDir+"NAMELIST_ICON")
;  Radiation = systemfunc("grep inwp_radiation "+FileDir+"NAMELIST_ICON")
;  print(Turbulence)



  
  if (IMan .eq. "T") then
   minlev    = stringtofloat(getenv("minlev"))
   maxlev    = stringtofloat(getenv("maxlev"))
   interval  = stringtofloat(getenv("interval"))
  end if

  
  leftlat    = stringtofloat(Leftlat)
  rightlat   = stringtofloat(Rightlat)
  leftlon    = stringtofloat(Leftlon)
  rightlon   = stringtofloat(Rightlon)


  if (Model .eq. "ICONAM" ) then
   if (var .eq. "TKE" .or. var .eq. "W" ) then
    height_var="ZH3"
    vlev="lev_2"
    extra_level=1
   else
    height_var="ZF3"
    vlev="lev"
    extra_level=0
   end if
  end if
 
  File = addfile(File_Int, "r")

  varnames = getfilevarnames(File)
  if(.not.any(ismissing(varnames))) then
    do i=0,dimsizes(varnames)-1
      printFileVarSummary (File,varnames(i))
    end do
  end if

  print(PlotDir)
  print(PlotFileN)
;  print(Fil)

;  print(var)
;  printVarSummary(ozone)

  var3dll = File->$var$(0,:,:,:)  ;time,lev,lat,lon
  ppmv2gg=1.e-6*47.9982/28.970
  var3dll = var3dll ;/ppmv2gg
  num_levels=dimsizes(var3dll&lev)-extra_level
  nlon=dimsizes(var3dll&lon)
  nlat=dimsizes(var3dll&lat)

  if (Model .eq. "ICONAM") then
   z3dll   = File->$height_var$(0,:,:,:) ; lev,lat,lon 
   topoll  = File->ZH3(0,num_levels,:,:)
  end if
  
  if (Model .eq. "ICOHDC") then
   z3dll   = File->ZF3(0,:,:,:)
   topoll  = File->PHIS(0,:,:)
   topoll  = topoll /9.80665
;   z3dll   = z3dll/9.80665
  end if 



; printVarSummary(var3dll)
; printVarSummary(z3dll)
; printVarSummary(topoll)

 
 ;File_Int has the variable and the height field in a gaussian grid, lon from 0 to 360


; Start cross section
;************************************
; calculate great circle along transect
;************************************

  npts     =   200                    ; number of points in resulting transect

  if (leftlat .ne. rightlat .and. leftlon .ne. rightlon) then

   if (leftlon .lt. 0.) then
    leftlon  = leftlon + 360.
    iu=1
   end if
   if (rightlon .lt. 0.) then
    rightlon = rightlon + 360.
   end if                             ; File_Int hat lon in the range 0 to 360, 
                                      ; because it comes from the interpolation to a 
                                      ; gaussian grid with cdo
   dist     = gc_latlon(leftlat,leftlon,rightlat,rightlon,npts,2*iu)
   points   = ispan(0,npts-1,1)*1.0 

   xlon = dist@gclon
   ylat = dist@gclat
   dx   = dist@spacing
   xunit_option=1     ;1 for km

  else

    if (leftlat .eq. rightlat .and.  leftlon .ne. rightlon) then

     ylat=new(npts,"float")
     ylat(:)=leftlat

     ; if (leftlon .gt. rightlon), better use the (180.,180.)range to 
     ;  calculate the longitudinal values of the points in the transect

     if(leftlon .gt. rightlon) then

      if (leftlon .gt. 180.) then
         leftlon = leftlon-360.
         iu=-1
         if (rightlon .gt. 180.)then
          rightlon = rightlon -360.
         end if
      end if
        
     end if

     minlon=min((/leftlon,rightlon/))
     maxlon=max((/leftlon,rightlon/))

      xlon=fspan(minlon,maxlon,npts)
      dx=(maxlon-minlon)/(npts-1.)
      xunit_option=0     ; 0  for degrees ; x axis will be in degrees
      xvarname="Lon"

    end if

    if (leftlon .eq. rightlon .and. leftlat .ne. rightlat) then
      xlon=new(npts,"float")
      xlon(:)=leftlon
      iu=-1
      minlat=min((/leftlat,rightlat/))
      maxlat=max((/leftlat,rightlat/))    

      ylat=fspan(minlat,maxlat,npts)
      dx=(maxlat-minlat)/(npts-1.)
      xunit_option=0     ; 0  for degrees
      xvarname="Lat"
    end if

    if (leftlon .eq. rightlon .and. leftlat .eq. rightlat) then
     print("it is not sense, give me two different (lat,lon) points")
     exit
    end if

  end if

;********************************
; interpolate data to great circle
;********************************
  
  xlon=where(xlon.lt.0.,xlon+360.,xlon)  ; transform again to (0.,360)

  var2dcs = linint2_points_Wrap(var3dll&lon,var3dll&lat(nlat-1:0),\
                                var3dll(:,nlat-1:0,:),True,\
                                xlon,ylat,2)

  zcoord  = linint2_points_Wrap(z3dll&lon,z3dll&lat(nlat-1:0),\
                                z3dll(:,nlat-1:0,:),True,\
                                xlon,ylat,2)

  topocs  = linint2_points_Wrap(topoll&lon,topoll&lat(nlat-1:0),\
                                topoll(nlat-1:0,:),False,\
                                xlon,ylat,2)

  zcoord  = zcoord/1000.
  topocs  = topocs/1000.   ; transform to km

  h = dimsizes(zcoord)
  zcoord_numz = h(0)
  zcoord_numx = h(1)

  ; x-Values
  xcoord = new( (/ zcoord_numz, zcoord_numx /),float)


  if (xunit_option .eq. 1) then
   x1=   fspan(0, 111.2*(npts-1)*dx, zcoord_numx)        ; 2*pi*Re/360 = 111.20 Km
  end if
  if (xunit_option .eq. 0) then
   if (xvarname .eq. "Lon")
    x1=xlon
    if (iu .eq. -1) then
     x1=where(x1 .gt. 180, x1-360.,x1)
    end if
   else
    x1=ylat
   end if
  end if

  do i=0, zcoord_numz-1
    xcoord(i,:) = x1
  end do

  ;plot
  plottingdir="/scratch/mpi/CC/mh0287/users/kristina/icon-dev/experiments/plots/"
  plotf="ozon_r2b4"

  wks = gsn_open_wks("ps",plottingdir+plotf)
  gsn_define_colormap(wks,"testcmap")

  ; Ressourcen-Variable :
  ; ---------------------

  res   = True   

  res@trGridType            = "TriangularMesh"
  res@sfXArray              = xcoord(:,:)
  res@sfYArray              = zcoord(:,:)

  ;res@trXMinF               =     ; in km  
  ;res@trXMaxF               =     ; in km  
  res@trYMinF               =   -0.3     ; in km  
  res@trYMaxF               = 50.     ; in km

  resT  = res

  res@gsnFrame = False ; do not advance frame
  res@gsnDraw = False ; do not draw until calling "draw"

  res@cnFillOn              = True     ; 'shaded plot' anschalten
  res@cnFillMode           = "CellFill"
  res@gsnSpreadColors = True
  res@cnLinesOn             = False     ; Turn on/off contour lines.
  res@cnLineLabelsOn        = False     ; Turn on/off line labels.
  res@cnFillOn        = True

; The levels are set manually if IMan is "T"
  if (IMan .eq. "T")
   res@cnLevelSelectionMode = "ManualLevels"
   res@cnMinLevelValF        = minlev
   res@cnMaxLevelValF        = maxlev
   res@cnLevelSpacingF       = interval 
  else
   res@cnLevelSelectionMode = "AutomaticLevels"
  end if

  res@lbLabelBarOn             = True
  res@lbLabelAutoStride = True
  res@lbOrientation     = "Vertical"

;  res@tiMainString  = "transect ("+leftlon+","+leftlat+") ("+rightlon+","+rightlat+")"
    res@tiMainString  = ""


  if (xunit_option .eq. 1) then
   res@tiXAxisString         = "x [km]"
  else
   res@tiXAxisString         = xvarname+" [degrees]"
  end if
  res@tiXAxisFontHeightF    = 0.02 ;
  res@tiXAxisFontThicknessF = 1.5

  res@tiYAxisString         = "z [km]"
  res@tiYAxisFontHeightF    = 0.02 ;
  res@tiYAxisFontThicknessF = 1.0
  res@tmYLLabelFontHeightF  = 0.02

  
;;;;;;;;;;;;;;;;;;;;;
; contour plot of var
;;;;;;;;;;;;;;;;;;;;;

  plot = gsn_csm_contour( wks, var2dcs, res )    ; 


  
; add the topography (taken from Guehnter)
; Stipple the terrain area and add a line at the surface
; ;;;;;;;;;;;;;;;;;;;;

  dimsr  = dimsizes(xcoord)
  npoints = dimsr(1)

  xpoly = new(npoints,typeof(xcoord))
  ypoly = new(npoints,typeof(topocs))

  xpoly(0:npoints-1) = (/xcoord(0,:)/)
  ypoly(0:npoints-1) = (/topocs/)

  resT@gsFillIndex      = 17           ; stipple
  resT@gsLineThicknessF = 2.0          ; twice as thick

  dum1 = gsn_add_polygon(wks,plot,xpoly,ypoly,resT)
  dum0 = gsn_add_polyline(wks,plot,xpoly,ypoly,resT)


;************************************************
; Draw text on plot using NDC coordinates
;************************************************
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.015                     ; font smaller. default big

;  gsn_text_ndc(wks,Turbulence,0.3,0.1,txres)
;  txres@txPosYF = 0.08                    ; font smaller. default big
;  gsn_text_ndc(wks,Radiation,0.3,0.1,txres)

  draw (plot)
  frame(wks) 

; system("gv "+PlotDir+PlotFileN+"."+PlotFormat+" &")

end
