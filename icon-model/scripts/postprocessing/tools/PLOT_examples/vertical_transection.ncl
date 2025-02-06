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
; icon_4.ncl
;
; Concepts illustrated:
;   - Plotting ICON model data
;   - Contouring one-dimensional X, Y, Z data
;   - Using triangular meshes to create contours
;   - Attaching filled polygons to a map
;   - Attaching a custom labelbar to a map
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;*************************************************
; Main code
;************************************************* 
begin
;  scale   = 1e6
  rad2deg = 45./atan(1.)   ; radians to degrees
  p00 = 1.e3  ; reference pressure in hPa

  Resolution1 = "R2B04L47"
  Resolution2 = "R2B06L31"
  Model       = "ICOHAM"

  FileName1 = "atm-Cels.nc"
;  FileName1 = "atm.nc"

  File  = addfile(FileName1, "r" )

  ntimestep = 0 ; time step selected for plotting (starts with 0!)
  nlevel = 40    ; model level selected for plotting (starts with 0!)
  num_levels = 47 ; total number of model levels

  x = File->clon *rad2deg   ; cell center, lon
  y = File->clat *rad2deg   ; cell center, lat

  topo = File->PHIS(0,:)   ; dims: (time,cell)
  topo = topo/9.81 ; convert from geopotential to geometric height
  phis = File->PHIS(0,:)   ; dims: (time,cell)
  topo@units=""
  topo@long_name=""

  psfc = File->PS(ntimestep,:)  ; surface pressure, dims: (time,cell)
  psfc = psfc*1.e-2  ; convert to hPa
  hyam = dble2flt(File->hyam*1.e-2)  ; convert to hPa
  hybm = dble2flt(File->hybm)
  hyai = dble2flt(File->hyai*1.e-2)  ; convert to hPa
  hybi = dble2flt(File->hybi)
  klevi = dimsizes(hyai)                ;
  klevm = dimsizes(hyam)                ;

  temp3d = File->T(ntimestep,:,:)   ; dims: (time,lev,cell)
  u3d    = File->U(ntimestep,:,:)   ; dims: (time,lev,cell)
  v3d    = File->V(ntimestep,:,:)   ; dims: (time,lev,cell)

  pres3d  = new ((/num_levels,  dimsizes(psfc)/),float)
  pres3di = new ((/num_levels+1,dimsizes(psfc)/),float)
  theta3d = new ((/num_levels,  dimsizes(psfc)/),float)

  do k = 0, num_levels-1
    pres3d(k,:) = psfc(:)*hybm(k) + hyam(k)
  end do
  do k = 0, num_levels
    pres3di(k,:) = psfc(:)*hybi(k) + hyai(k)
  end do

  ; potential temp. at lowest level
  theta3d = temp3d;*(p00/pres3d)^(2./7.)
  theta3d@units="K"
  theta3d@long_name="Temperature"

;  print(max(theta))
  ;var = var*scale

  x    = File->clon *rad2deg   ; cell center, lon
  y    = File->clat *rad2deg   ; cell center, lat
  vlat = File->clat_vertices * rad2deg
  vlon = File->clon_vertices * rad2deg
  vlon = where(vlon.lt.0, vlon + 360, vlon)

  ylat = fspan(-90,90,91)
  ylat@units = "degrees_north"
  xlon = fspan(-180,180,181)
  xlon@units = "degrees_east"
  ylat!0 = "ylat"
  xlon!0 = "xlon"

  u3dll     = triple2grid(x,y,u3d,xlon,ylat,False)
  v3dll     = triple2grid(x,y,v3d,xlon,ylat,False)
  pres3dll  = triple2grid(x,y,pres3d,xlon,ylat,False)
  pres3dill = triple2grid(x,y,pres3di,xlon,ylat,False)
  theta3dll = triple2grid(x,y,theta3d,xlon,ylat,False)
  phisll    = triple2grid(x,y,phis,xlon,ylat,False)

  u3dll!1 = "ylat"
  u3dll!2 = "xlon"
  u3dll&ylat  = ylat
  u3dll&xlon  = xlon

  phisll!0 = "ylat"
  phisll!1 = "xlon"
  phisll&ylat = ylat
  phisll&xlon = xlon

  copy_VarMeta( u3dll,  pres3dill)
  copy_VarMeta( u3dll,  pres3dll)
  copy_VarMeta( u3dll,  theta3dll)

 ; now we try to interpolate from model levels to pressure levels

  plevs = (ispan(100,1000,25))*1.0

  temp_at_p = new ((/dimsizes(plevs),dimsizes(ylat),dimsizes(xlon)/),double)
  copy_VarMeta( u3dll,  temp_at_p)
  temp_at_p@units = theta3d@units
  temp_at_p@long_name =theta3d@long_name

; define options for the foloowinf function

  extrapolate = True ; switch for doing extrapolation below the ground:
  intmethod   = 1 ; 1: method for temperature, -1: method for geopotential, 0: other vars
  tlow =  theta3dll(num_levels-1,:,:)  ; temperature at lowest model level
                                       ; (use ground temperature if available)

  temp_at_p = vinth2p_ecmwf(theta3dll,hyam,hybm,plevs,pres3dill(num_levels,:,:)*100.,\
                            intmethod,1.,1,extrapolate,1,tlow,phisll)
  printVarSummary(temp_at_p)

; Start cross section
;************************************
; calculate great circle along transect
;************************************
  leftlat        = 80.
  leftlatorient  = "N"
  rightlat       = -80.0
  rightlatorient = "N"

  leftlon        = 0.0
  leftlonorient  = "E"
  rightlon       = 0.0
  rightlonorient = "E"

  npts     =   100                    ; number of points in resulting transect

  dist     = gc_latlon(leftlat,leftlon,rightlat,rightlon,npts,2)
  points   = ispan(0,npts-1,1)*1.0
  dist@gclon = where(dist@gclon.ge.180, dist@gclon-360, dist@gclon)
;********************************
; interpolate data to great circle
;********************************
  trans   = linint2_points_Wrap(temp_at_p&xlon,temp_at_p&ylat,temp_at_p,False,\
                                dist@gclon,dist@gclat,2)

;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

; Graphic Labeling
  TitleString  = "ICON coupled aqua planet, ATMO"
  LeftString   = Model+" "+Resolution1
  RightString  = theta3d@long_name
  CenterString = ""
  YAxisString  = "Pressure (hPa)"
  LabelString  = "["+theta3d@units+"]"
  LabelString  = "[C]"

  GraphicFile  = "transection"
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

; Start the graphics.
;
  wks = gsn_open_wks("png",GraphicFile)       ; open a ps file
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200") ; choose colormap

  res                     = True          ; plot mods desired
  res@tmXBMode            = "Explicit"    ; explicitly label x-axis
  res@tmXBValues          = (/points(0),points(npts-1)/) ; points to label
; label values
  res@tmXBLabels          = (/leftlat + leftlatorient+", "+leftlon+leftlonorient,rightlat+rightlatorient+", "+rightlon+rightlonorient/)

  res@cnFillOn            = True         ; turn on color
  res@lbLabelAutoStride   = True         ; nice label bar label stride
  res@gsnSpreadColors     = True         ; use full range of colormap
  res@cnLinesOn           = False        ; turn off countour lines
;  res@lbOrientation       = "horizontal" ; vertical label bar
  res@lbOrientation       = "vertical"   ; vertical label bar
; res@pmLabelBarHeightF        = 0.08
; res@pmLabelBarWidthF         = 0.5
 res@pmLabelBarOrthogonalPosF = -0.04    ; move label bar closer to plot

  res@tiMainString        = TitleString  ; add title
  res@trYReverse          = True         ; reverse y axis
  res@tiYAxisString       = YAxisString

  res@gsnLeftString     = LeftString
  res@gsnRightString    = RightString
  res@lbTitlePosition   = "Bottom"
  res@lbTitleString     = LabelString

  plot = gsn_csm_contour(wks,trans,res)  ; create plot

;********************************
; show transect on a map
;********************************
   i               = NhlNewColor(wks,0.8,0.8,0.8) ;add gray to continents
  mres             = True                ; plot mods desired
  mres@gsnFrame    = False               ; don't turn page yet
; mres@gsnDraw     = False               ; don't draw yet
  mres@mpFillOn    = True
  mres@mpFillColor = "Background"
  mres@sfXArray                 = x  ; These are 1D arrays, so a triangular
  mres@sfYArray                 = y  ; mesh will be created internally.

  mres@tiMainString = "Transect Location" ; title

  mapplot = gsn_csm_contour_map(wks,temp3d(30,:),mres)         ; create map

; add polyline
  pres                  = True           ; polyline mods desired
  pres@gsLineColor      = "red"          ; color of lines
  pres@gsLineThicknessF = 3.0            ; line thickness
  gsn_polyline(wks,mapplot,(/leftlon,rightlon/),(/leftlat,rightlat/),pres)

  frame(wks)
end

