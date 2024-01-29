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
  vlev="mlev"          ; it is ilev in case var="W"
  iu=1    ; if +1, then lon between 0,360  ; if  -1, then  lon between -180 and 180
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; get the variables from the driver
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  Model      = getenv("Model")
  File_Int   = getenv("file_int")
  var        = getenv("var")
  Rightlat    = getenv("rightlat")
  Rightlatp10 = getenv("rightlatp10")
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
   
  Turbulence = systemfunc("grep inwp_turb "+FileDir+"NAMELIST_ICON")
  Radiation = systemfunc("grep  inwp_radiation "+FileDir+"NAMELIST_ICON")
  
   if (var .eq. "TKE" ) then
    height_var="ZH3"
    vlev="ilev"
    extra_level=1
   else
    height_var="ZF3"
    vlev="mlev"
    extra_level=0
   end if
  end if
 
  File = addfile(File_Int, "r")

  var1dt = File->$var$(:,:,1)  ;tim,elev,lat,lon
  
  dNames=getvardims(File->$var$)
;  vlev=var1dt&'dNames(1)'
  
  print(dNames)
;  printVarSummary(var1dt)
;  print(var1dt&time)
;  print(var1dt&vlev)

  num_levels=dimsizes(var1dt&$vlev$)-extra_level
;  nlon=dimsizes(var1dt&lon)
; nlat=dimsizes(var1dt&lat)

  if (Model .eq. "ICONAM") then
   z3dll   = File->$height_var$(:,:,1) ; lev,lat
   topoll  = File->ZH3(0,num_levels,1)
  end if
  
  if (Model .eq. "ICOHDC") then

    Turbulence = systemfunc("grep lturb "+FileDir+"NAMELIST_ICON")
    Radiation = systemfunc("grep  lrad  "+FileDir+"NAMELIST_ICON")
  
   z3dll   = File->ZF3(:,:,1)
   topoll  = File->PHIS(:,1)
   topoll  = topoll /9.80665
  end if

  dim = dimsizes(var1dt)
  tdim = dim(0)
  tstart = 0
  tend   = tdim-1
  xtime = fspan(1,tend*600,tdim)
  print(xtime)
  var1dt&time = xtime

; printVarSummary(var1dt)
; printVarSummary(z3dll)
; printVarSummary(topoll)

 
 ;File_Int has the variable and the height field in a gaussian grid, lon from 0 to 360

  ;plot

  wks = gsn_open_wks("ps",PlotDir+PlotFileN)
  gsn_define_colormap(wks,"gui_default")

  ; Ressourcen-Variable :
  ; ---------------------

  res   = True

  res                          = True                ; plot mods desired
  res@cnInfoLabelOn            = False              ; Turn off informational
  res@cnFillOn                 = True                ; color on
  res@cnFillMode               = "CellFill" ; "AreaFill", "CellFill" or "RasterFill"
  res@gsnSpreadColors          = True 
  res@lbLabelStride            = 2                  ; every other label
  res@lbOrientation            = "Vertical"         ; vertical label bar
  res@pmLabelBarOrthogonalPosF = -0.01              ; move label bar closer to axis
  res@cnLinesOn                = False              ; turn off contour lines
  res@cnLineLabelsOn           = False              ; Turn off contour line labels
  res@lbAutoManage             = False
  res@lbLabelAutoStride        = True
  res@lbTitleFontHeightF       = 0.015
  res@lbTitleOn                = True               ; turn on title
  res@lbTitleString            = var1dt@units
  res@lbTitlePosition         = "Top"
  res@lbTitleOffsetF           = 0.01               ; more space between label bar and text
;  res@pmLabelBarHeightF        = 0.08
;  res@pmLabelBarWidthF         = 0.1
;  res@pmLabelBarOrthogonalPosF = 0.08
;  res@cnLevelSelectionMode     = "ExplicitLevels"

  resT  = res

  res@gsnFrame = False ; do not advance frame
  res@gsnDraw = False ; do not draw until calling "draw"


; The levels are set manually if IMan is "T"
  if (IMan .eq. "T")
   res@cnLevelSelectionMode = "ManualLevels"
   res@cnMinLevelValF        = minlev
   res@cnMaxLevelValF        = maxlev
   res@cnLevelSpacingF       = interval 
  else
   res@cnLevelSelectionMode = "AutomaticLevels"
  end if

    res@tiMainString  = ""

  res@tiXAxisFontHeightF    = 0.02 ;
  res@tiXAxisFontThicknessF = 1.5

  res@tiYAxisFontHeightF    = 0.02 ;
  res@tiYAxisFontThicknessF = 1.0
  res@tmYLLabelFontHeightF  = 0.02


;  --------------------------------------------------------------------------------
; plot-specific ressources
;---------------------------------------------------------------------------------

;    res@tiYAxisString           =  "Height (km)" ;plev@long_name
    res@trYReverse              = True                ; flip    the Y-axis
;    res@trYMinF               =   -0.3     ; in km  
;    res@trYMaxF               = 10.     ; in km
    res@tiXAxisString           = "time s" 
    res@sfXArray                = xtime
;;;;;;;;;;;;;;;;;;;;;
; contour plot of var
;;;;;;;;;;;;;;;;;;;;;

  plot = gsn_csm_contour( wks, var1dt(ilev|40:60,time|:), res )    ; 

;************************************************
; Draw text on plot using NDC coordinates
;************************************************
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.015                     ; font smaller. default big

  gsn_text_ndc(wks,Turbulence,0.3,0.1,txres)
  txres@txPosYF = 0.08                    ; font smaller. default big
  gsn_text_ndc(wks,Radiation,0.3,0.1,txres)

  draw (plot)
  frame(wks) 

 system("gv "+PlotDir+PlotFileN+"."+PlotFormat+" &")

end
