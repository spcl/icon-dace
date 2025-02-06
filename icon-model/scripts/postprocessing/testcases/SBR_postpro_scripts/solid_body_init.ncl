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
; plot_solid_body.ncl
;
; Purpose:
; Visualizing data of ICON solid body rotation test case
;
; Concepts illustrated:
;   - Plotting ICON model data
;   - Contouring one-dimensional X, Y, Z data
;   - Using triangular meshes to create contours
;
;History:
;
;Initial version by Daniel Rainert
;Modified version by Constantin Junk (MPI-M) (2010-10-18)


;*************************************************
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;---------------------------------------------------------------
;

begin

;*****************************************************************
;******************** get and define variables *******************
;******************************************************************

;***** get environment variables *****

  fnam      =  getenv("FNAM") 
  gridfile  =  getenv("GRIDFILE")
  diri      =  getenv("DIRI")
  diro      =  getenv("DIRO")
  pltType   =  getenv("PFMT")
  Resolution=  getenv("RESOLUTION")
  outfile   =  diro+getenv("PNAME")
  ConfigStr =  getenv("CONFIG_STRING")
  Model     =  getenv("MODEL")
  nlev      = stringtointeger(getenv("VERTICAL_RES"))-1
  ivel      = stringtointeger(getenv("ivel"))




;***** get file suffix *****

  fsfx      = get_file_suffix(fnam,0)
  pltTitle  = fsfx@fBase


;***** define ncl script variables ****

  colormap  = "testcmap"
  colormap2 = "rainbow+white+gray"
  colormap3 = "WhBlGrYeRe"
  colormap4 = "nrl_sirkes"
  FontHeight= 0.015
  scale_u   = 1.
  scale_v   = 1.
  scale_q4  = 1.
  rad2deg   = 45./atan(1.)   ; radians to degrees

;*****************************************************************
;***** read meteorological field and grid information *****
;*****************************************************************

  GridFile  = addfile(gridfile, "r" )
  File      = addfile(diri+fnam, "r" )

  VarName1 = "U"
  VarName2 = "V"
  VarName4 = "Q4"  ; contains numeric solution

;***** read variable 'time' of data file *****
  if (isfilevar(File, "time")) then 
    time = File->time(:)
    dims_time = dimsizes(time)
  end if

;***** choose time step *****
  ntimes = dims_time(0)
  ntime  = ntimes-1

;***** read meteorological field *****

;***** a) velocity U
  if (isfilevar(File, VarName2)) then 
    u = File->$VarName1$(ntime,nlev,:)   ; dims: (time,lev,cell)
    u = u*scale_u
    max_u = max(u)
    min_u = min(u)
    print("reading " + VarName1 + " complete" )
    print("max " + VarName1 + " " + max_u )
    print("min " + VarName1 + " " + min_u )
  end if

;***** b) velocity V
  if (isfilevar(File, VarName2)) then 
    v = File->$VarName2$(ntime,nlev,:)   ; dims: (time,lev,cell)
    v = v*scale_v
    max_v = max(v)
    min_v = min(v)
    print("reading " + VarName2 + " complete" )
    print("max " + VarName2 + " " + max_v )
    print("min " + VarName2 + " " + min_v )
  end if

;***** c) read numerical solution Q4 
  if (isfilevar(File, VarName4)) then 
    q4 = File ->$VarName4$(ntime,nlev,:) ; dims: (time,lev,cell)
    q4 = q4*scale_q4
    max_q4 = max(q4)
    min_q4 = min(q4)
    q4@long_name = "Q4 tracer"
    print("reading " + VarName4 + " complete" )
    print("max " + VarName4 + " " + max_q4 )

    ; read initialization
    ; Q4 (t=0)
    q4_ini = File ->$VarName4$(0,nlev,:) ; dims: (time,lev,cell)
    q4_ini = q4_ini*scale_q4
    q4_ini@long_name = "Q4 tracer"
    max_q4_ini = max(q4_ini)
    min_q4_ini = min(q4_ini)
  end if

;***** read numerical solution


  
;***** global grid
  x = File->clon *rad2deg   ; cell center, lon
  y = File->clat *rad2deg   ; cell center, lat

;---------------------------------------------------------------

;***************************************************************
;**************************  make plots ************************
;***************************************************************
;
; Set up resources for contour/map plot.
  ResC                = True

  ResC@gsnFrame             = False
  ResC@gsnSpreadColors      = True
  ResC@gsnLeftString        = "Reconstructed zonal wind"
  ResC@gsnCenterString      = " "
  ResC@gsnRightString       = ConfigStr 
  ResC@cnFillOn             = True
  ResC@cnLinesOn            = True
  ResC@cnInfoLabelOn        = False
  ResC@cnFillMode           = "Rasterfill"
  ResC@cnRasterSmoothingOn  = True

  ResC@tiMainString         = " "
  ResC@tiXAxisFontHeightF   = FontHeight
  ResC@tiYAxisFontHeightF   = FontHeight
  ResC@tmXBLabelFontHeightF = FontHeight
  ResC@tmYLLabelFontHeightF = FontHeight
  ResC@gsnStringFontHeightF = FontHeight

  ResC@mpProjection         = "CylindricalEquidistant"
  ResC@mpPerimOn            = False
  ResC@mpFillOn             = False


  ResC@mpLimitMode          = "LatLon"
  ResC@mpCenterLonF         = -60.
  ResC@mpCenterLatF         = 0.
  ResC@mpMinLatF            = -90.
  ResC@mpMaxLatF            =  90.
  ResC@mpGreatCircleLinesOn = True
  ResC@gsnDraw              = False
  ResC@sfXArray             = x  ; These are 1D arrays, so a triangular
  ResC@sfYArray             = y  ; mesh will be created internally.

  ResC@lbLabelBarOn         = True 
  ResC@lbLabelFontHeightF   = 0.0145
  ResC@lbLabelAutoStride    = True
  ResC@lbTitlePosition      = "Bottom" 
  ResC@lbTitleFontHeightF   = 0.015

  ResC@pmLabelBarHeightF        = 0.08
  ResC@pmLabelBarWidthF         = 0.65
  ResC@pmLabelBarOrthogonalPosF = 0.13

  ResC@mpGridAndLimbOn          = True      ; turn on lat/lon lines
  ResC@mpGridLineDashPattern    = 2         ; make lat/lon lines dash
  ResC@mpLabelsOn               = False      ; labels identifying the equator ("EQ"),..
  ResC@mpLabelFontHeightF       =0.011      ; labels font height


; Annotation resources
  txres                = True       
  txres@txJust         = "CenterLeft"
  txres@txFontHeightF  = 0.014              ; font smaller. default big
;
; 1.) show velocity field in CylindricalEquidistant projection
;***************************************************************

  if (ivel .eq. 1) then
  wks = gsn_open_wks(pltType,outfile+"_velocity")
  gsn_define_colormap(wks,colormap)


  ResC@lbTitleString        = u@long_name+" in "+u@units






;******* plot zonal and meridional velocity component *********
  plot_v = new(2,graphic)

; automatically create nice min/max/ci values for blue/red colortable
  symMinMaxPlt (u,30,False,ResC)
  plot_v(0) = gsn_csm_contour_map_ce(wks,u,ResC)

  symMinMaxPlt (v,30,False,ResC)
  ResC@lbTitleString            = v@long_name+" in "+v@units
  ResC@gsnLeftString        = "Reconstructed meridional wind"
  plot_v(1) = gsn_csm_contour_map_ce(wks,v,ResC)

  resP                       = True
  resP@gsnPanelFigureStrings = (/"U","V"/)
  gsn_panel(wks,plot_v,(/1,2/),resP)

  delete(wks)
  end if

;***************************************************************

;****************************************************************
; 2.) show initial conditions, the numerically calculated tracer field 
;     in CylindricalEquidistant projection and the difference between 
;     numerical and analytical solution (difference between q4 and q4_ini)
;***************************************************************

  wks = gsn_open_wks(pltType,outfile+"_tracer")

  FontHeight= 0.017

;****** Zoom in ******
  ResC@mpCenterLonF         = -60.
  ResC@mpCenterLatF         = 0.
  ResC@mpMinLatF            = -30.
  ResC@mpMaxLatF            =  30.
  ResC@mpMinLonF            = -120.
  ResC@mpMaxLonF            = -60.

  gsn_merge_colormaps(wks,colormap3,colormap)
  ResC@gsnSpreadColorStart = 0       ; first color to use
  ResC@gsnSpreadColorEnd   = 101     ; last color to use
  ResC@cnMaxLevelValF      = 1.05
  ResC@cnMinLevelValF      = -0.05

  ResC@lbTitleString            = q4@long_name  +" in "+q4@units
  ResC@cnLevelSelectionMode     = "AutomaticLevels"
  ResC@cnLevelSpacingF          = 0.1

  plot_q = new(3,graphic)

  ResC@tiMainString         = " "
  ResC@gsnLeftString        = "Initial conditions"
;***** plot variable for initial conditions *****
  plot_q(0) = gsn_csm_contour_map_ce(wks,q4_ini,ResC)

  text1=gsn_add_text(wks,plot_q(0),"max: "+sprintf("%10.4e", max_q4_ini),-119.4,-28.0,txres)
  text2=gsn_add_text(wks,plot_q(0),"min: "+sprintf("%10.4e", min_q4_ini),-119.4,-25.0,txres)

  ResC@tiMainString         = " "
  ResC@gsnLeftString        = "Numerical solution"

;***** plot variable for numerically calculated tracer field *****
  plot_q(1) = gsn_csm_contour_map_ce(wks,q4,ResC)

  text3=gsn_add_text(wks,plot_q(1),"max: " + sprintf("%10.4e", max_q4),-119.4,-28.0,txres)
  text4=gsn_add_text(wks,plot_q(1),"min: " + sprintf("%10.4e", min_q4),-119.4,-25.0,txres)

;***** calculate and plot difference between numerical and analytical sol *****
  q4mq4_ini = q4 - q4_ini
  copy_VarMeta (q4, q4mq4_ini)
  q4mq4_ini@long_name = "Q4 - Q4_ini"
  q4mq4_ini@name     = "tracer concentration in kg/kg"

  max_err = max(q4mq4_ini)
  min_err = min(q4mq4_ini)

; go back to standard color map
; gsn_define_colormap(wks,colormap4)

  ResC@gsnLeftString        = "Error: q4(t=12d) - q4(t=0d)"
  ResC@tiMainString         = " "
  ResC@lbTitleString        = q4mq4_ini@name ;+" in "+q4m6@units
  ResC@cnLevelSpacingF      = 0.01
  ResC@gsnSpreadColorStart  = 146         ; first color to use
  ResC@gsnSpreadColorEnd    = 256       ; last color to use
  ResC@cnLevelSelectionMode = "ManualLevels"
  ResC@cnMaxLevelValF       = 0.065
  ResC@cnMinLevelValF       = -0.065

;***** plot variable for difference btw numerical and analytical sol *****
  plot_q(2) = gsn_csm_contour_map(wks,q4mq4_ini,ResC)

  text5=gsn_add_text(wks,plot_q(2),"max: " + sprintf("%10.4e", max_err),-119.4,-28.0,txres)
  text6=gsn_add_text(wks,plot_q(2),"min: " + sprintf("%10.4e", min_err),-119.4,-25.0,txres)

;**** plot inital conditions (plot_q(0)), result (plot_q(1)) 
;**** and error (plot_q(2)) on one panel

  panel_res                       = True
  panel_res@gsnPanelFigureStrings = (/"t=0d","t=12d","err"/)

  gsn_panel(wks,plot_q,(/1,3/),panel_res)

end

