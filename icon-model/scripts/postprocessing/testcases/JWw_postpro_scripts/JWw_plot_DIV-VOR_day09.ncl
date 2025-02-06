; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;------------------------------------
; Hui Wan (MPI-Met, 2008-07-05)
; Hui Wan (MPI-Met, 2008-07-13)
; Hui Wan (MPI-Met, 2008-07-15)
; Hui Wan (MPI-Met, 2009-04-15)
; Hui Wan (MPI-Met, 2009-04-28: Modification for the hex. model)
; Hui Wan (MPI-Met, 2009-05-12: Modification for ECHAM)
;------------------------------------
; Script type: visualization
;---------------------------------------------------------------
; This script makes contour plots of 850 hPa vorticity, divergence
; and omega for the Jablonowski-Williamson steady state test.
;---------------------------------------------------------------
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;---------------------------------------------------------------
 begin

  PlotFormat = getenv("plot_file_format")
  Model      = getenv("Model")
  DataPath   = getenv("DataPath")
  PlotPath   = getenv("PlotPath")
  Resolution = getenv("Resolution")
  ConfigStr  = getenv("ConfigStr")
  ExpName    = getenv("ExpName")
  plotBaseName = getenv("plotBaseName")
  VarName    = getenv("JWw_VarName")
  CellType   = getenv("cell_type")
  OutFrq     = stringtoint(getenv("output_frequency"))
  VorAtCorners = stringtoint(getenv("vorticity_at_corners"))

  rad2deg = 45./atan(1.)   ; radians to degrees

;---------------------------------------------------------------
; names of the data file and the plot
;---------------------------------------------------------------

  DataFileName = DataPath + ExpName +"_"+VarName+"850.nc"
  PlotFileName = PlotPath + plotBaseName +"_"+VarName+"850_day09"

;---------------------------------------------------------------
; time steps to be plotted
;---------------------------------------------------------------

  iday  = (/9/)               ; days to plot
  ndays = dimsizes(iday)
  plot  = new( ndays,graphic )
  dum   = plot

; time indices in the model output

  if (Model.eq."ECHAM5") then
     tid = OutFrq*iday -1  ; ECHAM doesn't write out initial data
  else
     tid = OutFrq*iday
  end if

;---------------------------------------------------------------
; lat./lon.
;---------------------------------------------------------------

  if (Model.eq."ICOHDC") then
     if (CellType.eq."3") then
        if ((VarName.eq."VOR").and.(VorAtCorners.eq.1)) then
           latName    = "vlat" 
           lonName    = "vlon"                                            
        else                                                              
           latName    = "clat"                                            
           lonName    = "clon"                                            
        end if                                                            
     else                                                                 
     if (CellType.eq."6") then                                            
        if ((VarName.eq."VOR").and.(VorAtCorners.eq.1)) then
           latName    = "clat"                                            
           lonName    = "clon"                                            
        else                                                              
           latName    = "vlat"
           lonName    = "vlon"                                            
        end if                                                            
     else                                                                 
        print("Wrong cell type. Should be 3 or 6")                        
     end if                                                               
     end if                                                               
  end if                                           

;---------------------------------------------------------------
; labels and some variable dependent settings
;---------------------------------------------------------------

  if (VarName.eq."VOR") then

     varString   = "850 hPa vorticity (10~S~-5~N~s~S~-1~N~) at day "+iday
     scale      = 1e5

     colormap   = "BlWhRe"
     colorStart = 38
     colorEnd   = 97

     varMin     = (/ -12./)
     varMax     = (/  30./)
     varInt     = (/  3./)
     clevels    = (/-3,0,3/)

;     if (Model.eq."ECHAM5") then
        vname = "svo"  ; var name given by afterburner 
;     else
;        vname = "ha_diag_rel_vort_c"
;     end if

  else 
  if (VarName.eq."DIV") then ;divergence

     varString   = "850 hPa divergence (10~S~-6~N~ s~S~-1~N~) at day "+iday 
     scale      = 1e6

     colormap   = "testcmap"
     colorStart = 2 
     colorEnd   = 200 

     varMin     = (/ -15./)
     varMax     = (/  15./)
     varInt     = (/  1.5/)
     clevels    = (/-1.5,1.5/)

;     if (Model.eq."ECHAM5") then
        vname = "sd"  ; var name given by afterburner
;     else
;        vname = "ha_diag_div"
;     end if

  else ; omega

     varString  = "850 hPa omega (Pa s~S~-1~N~) at day "+iday 
     scale      = 1.

     colormap   = "testcmap"
     colorStart = 2 
     colorEnd   = 200 

     varMin     = (/ -0.4/)
     varMax     = (/  0.2/)
     varInt     = (/  0.05/)
     clevels    = (/-0.05,0.05/)

     if (Model.eq."ECHAM5") then
        vname = "var135"  ; var name given by afterburner
     else
        vname = "wap"
     end if
  end if
  end if

  CenterString = " "
  expString  = ConfigStr ;+" "+Resolution 

;---------------------------------------------------------------
; open data file and get grid information
;---------------------------------------------------------------

  File = addfile( DataFileName, "r" )

  if (Model.eq."ICOHDC") then
     x = File->$lonName$ *rad2deg
     y = File->$latName$ *rad2deg
  end if

;---------------------------------------------------------------
; check whether the simulation is long enough to make the plots
;---------------------------------------------------------------

  nsteps_in_data_file = dimsizes(File->time)

  if (Model.eq."ECHAM5") then
     ndays_in_data_file =  nsteps_in_data_file/OutFrq
  else
     ndays_in_data_file = (nsteps_in_data_file -1.)/OutFrq
  end if

  nmax_needed         = max(iday)

  if ( ndays_in_data_file .lt. nmax_needed ) then
     print("Problem encountered: found only "+ ndays_in_data_file +\
           " day(s) of data in the model output")
     print("Can not plot the flow field at day "+ nmax_needed )
     exit
  end if

;------------------------------------------------------------------
; make plots
;------------------------------------------------------------------
; This time we overlay a color shading plot by a pure contour line plot

  wks = gsn_open_wks(PlotFormat,PlotFileName)
  gsn_define_colormap(wks,colormap)

;------------------------------------------------------------
; the commong setups

  ResC                = True

  ResC@gsnDraw        = False 
  ResC@gsnFrame       = False

  ResC@cnInfoLabelOn  = False
  ResC@cnLineLabelsOn = False

  FontHeight = 0.011

  ResC@tiXAxisFontHeightF    = FontHeight
  ResC@tiYAxisFontHeightF    = FontHeight
  ResC@tmXBLabelFontHeightF  = FontHeight
  ResC@tmYLLabelFontHeightF  = FontHeight
  ResC@gsnStringFontHeightF  = FontHeight +0.002

  ResC@tmXBLabelJust  = "CenterCenter"

  if (Model.eq."ICOHDC") then
     ResC@sfXArray = x
     ResC@sfYArray = y
  end if

  ResC2               = ResC

;----------------------------------
; shadings vs contours

  ResC@cnFillOn       = True
  ResC@cnLinesOn      = False 

  ResC2@cnFillOn       = False
  ResC2@cnLinesOn      = True 

;-----------------------------------
; label bar for the underlying plot 

  ResC@lbLabelBarOn             = True 
  ResC@pmLabelBarHeightF        = 0.04
  ResC@pmLabelBarOrthogonalPosF = 0.25 
  ResC@lbLabelAutoStride        = True
  ResC@lbLabelFontHeightF       = FontHeight +0.001
  ResC@pmLabelBarWidthF         = 0.7 

  ResC2@lbLabelBarOn            = False 

;-----------------------------------
; contour levels

  ResC@cnLevelSelectionMode = "ManualLevels"
  ResC@gsnSpreadColors      = True
  ResC@gsnSpreadColorStart  = colorStart
  ResC@gsnSpreadColorEnd    = colorEnd

  ResC2@cnLevelSelectionMode = "ExplicitLevels"
  ResC2@cnLevels             = clevels

;-----------------------------------
; strings only on the contour plot

  ResC@gsnCenterString = ""
  ResC@gsnLeftString  = ""
  ResC@gsnRightString = ""

  ResC2@gsnStringFontHeightF  = 7*FontHeight +0.002
  ResC2@gsnCenterString       = CenterString
  ResC2@gsnLeftString         = varString 
  ResC2@gsnRightString        = expString

  ResC2@gsnLeftStringOrthogonalPosF  = 0.06
  ResC2@gsnRightStringOrthogonalPosF = 0.075

;------------------------------------------------------------
; the color shading plot is a "map plot"

  ResC@mpProjection           = "CylindricalEquidistant"
  ResC@mpGeophysicalLineColor = "transparent"
  ResC@mpFillOn               = False
  ResC@mpLimitMode            = "LatLon"
  ResC@mpMinLatF              = 25. 
  ResC@mpMaxLatF              = 75. 
  ResC@mpMinLonF              = 0.
  ResC@mpMaxLonF              = 360.
  ResC@mpCenterLonF           = 180.

  ResC@gsnMajorLatSpacing = 15. 
  ResC@gsnMajorLonSpacing = 30. 
  ResC@gsnMinorLatSpacing =  5. 
  ResC@gsnMinorLonSpacing = 10.

;---------------------------------------------------------------
; read data at each time step and make contour plot
;---------------------------------------------------------------

  do it = 0,ndays-1

     if (Model.eq."ICOHDC") then
        var = File->$vname$(tid(it),0,:)
     else
        var = File->$vname$(tid(it),0,:,:)
     end if

     var = var*scale

     ResC@cnMinLevelValF  = varMin(it) 
     ResC@cnMaxLevelValF  = varMax(it) 
     ResC@cnLevelSpacingF = varInt(it) 

     if (VarName.eq."VOR") then
        delete(ResC@gsnSpreadColors)
        ResC@cnFillColors = (/28,32,36,44,50, 54,57,62,66,70,74,79,84,90,97,102/)
     end if

     if (VarName.eq."OMEGA") then
        delete(ResC@gsnSpreadColors)
        ResC@cnFillColors = (/ 2,14,26,38,50,62,74,86,95, 105,140,160,175,200/)
     end if

     plot(it) = gsn_csm_contour_map(wks,var,ResC)
     dum (it) = gsn_csm_contour    (wks,var,ResC2)
     overlay( plot(it),dum(it) )
     
  draw(plot(it))
  frame(wks)
  end do

end

