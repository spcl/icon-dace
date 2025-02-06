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
; Hui Wan (MPI-Met, 2009-04-28: Modification for the hex. model)
; Hui Wan (MPI-Met, 2009-05-12: Modification for ECHAM)
;------------------------------------
; Script type: visualization
;---------------------------------------------------------------
; This script makes contour plots of 850 hPa vorticity and 
; divergence for the Jablonowski-Williamson steady state test.
;---------------------------------------------------------------
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
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
  PlotFileName = PlotPath + plotBaseName +"_"+VarName+"850"

;---------------------------------------------------------------
; time steps to be plotted
;---------------------------------------------------------------

  iday  = (/7,9/)             ; days to plot
  ndays = dimsizes(iday)
  plot  = new( ndays,graphic )

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

     colormap   = "BlRe"
     scale      = 1e5

     varMin     = (/ -3., -10./)
     varMax     = (/  5.,  35./)
     varInt     = (/  1.,   5./)

     varString   = "850 hPa vor. (10~S~-5~N~s~S~-1~N~)"

;     if (Model.eq."ECHAM5") then
        vname = "svo"  ; var name given by afterburner 
;     else
;        vname = "ha_diag_rel_vort_c"
;     end if

  else 
  if (VarName.eq."DIV") then ;divergence

     VarString  = " D "
     colormap   = "ViBlGrWhYeOrRe"
     colormap   = "testcmap"
     scale      = 1e6

     varMin     = (/ -5., -15./)
     varMax     = (/  5.,  15./)
     varInt     = (/  1.,   3./)

     varString   = "850 hPa div. (10~S~-6~N~ s~S~-1~N~) at day "+iday 

;     if (Model.eq."ECHAM5") then
        vname = "sd"  ; var name given by afterburner 
;     else
;        vname = "ha_diag_div"
;     end if

  else ; omega

     colormap   = "ViBlGrWhYeOrRe"
     colormap   = "testcmap"
     scale      = 1.

     varMin     = (/ -0.2, -0.6/)
     varMax     = (/  0.1,  0.2/)
     varInt     = (/ 0.03,  0.1/)

     varString   = "850 hPa omega (Pa s~S~-1~N~) at day "+iday 

     if (Model.eq."ECHAM5") then
        vname = "var135"  ; var name given by afterburner 
     else
        vname = "wap"
     end if
  end if
  end if

  CenterString = " "
  expString  = ConfigStr     ;+" "+Resolution 

;---------------------------------------------------------------
; horizontal domain to plot
;---------------------------------------------------------------

  MinLon = (/ 90., 120./)
  MaxLon = (/210., 270./)
  CtrLon = (/150., 195./)

  vpX = (/0.06 ,0.52/)
  vpW = (/0.36,0.45/)
  vpY = (/0.6, 0.6 /)
  vpH = vpW/ (/4,5/) *5/3 

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

  nmax_needed = max(iday)

  if ( ndays_in_data_file .lt. nmax_needed ) then
     print("Problem encountered: found only "+ ndays_in_data_file +\
           " day(s) of data in the model output")
     print("Can not plot the flow field at day "+ nmax_needed )
     exit
  end if

;---------------------------------------------------------------
; make plots
;---------------------------------------------------------------

  wks = gsn_open_wks(PlotFormat,PlotFileName)
  gsn_define_colormap(wks,colormap)

  ResC                = True

  ResC@gsnDraw        = True 
  ResC@gsnFrame       = False

  ResC@cnFillOn       = True
  ResC@cnLinesOn      = False
  ResC@cnInfoLabelOn  = False

  FontHeight = 0.012

  ResC@tiXAxisFontHeightF        = FontHeight
  ResC@tiYAxisFontHeightF        = FontHeight
  ResC@tmXBLabelFontHeightF      = FontHeight
  ResC@tmYLLabelFontHeightF      = FontHeight
  ResC@gsnStringFontHeightF      = FontHeight +0.002

  ResC@tmXBLabelJust  = "CenterCenter"

  ResC@mpProjection   = "CylindricalEquidistant"
  ResC@mpLimitMode    = "LatLon"
  ResC@mpMinLatF      = 25. 
  ResC@mpMaxLatF      = 75. 

  ResC@mpGeophysicalLineColor = "transparent"
  ResC@mpFillOn               = False

  ResC@gsnMajorLatSpacing = 15. 
  ResC@gsnMajorLonSpacing = 30. 
  ResC@gsnMinorLatSpacing =  5. 
  ResC@gsnMinorLonSpacing = 10.

  if (Model.eq."ICOHDC") then
     ResC@sfXArray = x
     ResC@sfYArray = y
  end if

  ResC@lbLabelBarOn             = True 
  ResC@pmLabelBarHeightF        = 0.04
  ResC@pmLabelBarOrthogonalPosF = 0.2 
  ResC@lbLabelAutoStride        = True
  ResC@lbLabelFontHeightF       = FontHeight

  ResC@cnLevelSelectionMode = "ManualLevels"
  ResC@gsnSpreadColors      = True

  ResC@gsnCenterString = CenterString

;---------------------------------------------------------------
; read data of each time step and make contour plot
;---------------------------------------------------------------

  do it = 0,ndays-1

     if (Model.eq."ICOHDC") then
        var = File->$vname$(tid(it),0,:)
     else
        var = File->$vname$(tid(it),0,:,:)
     end if

     var = var*scale

     ResC@vpXF            = vpX(it)        
     ResC@vpWidthF        = vpW(it)           
     ResC@vpYF            = vpY(it) 
     ResC@vpHeightF       = vpH(it) 

     ResC@pmLabelBarWidthF  = ResC@vpWidthF *0.9

     ResC@gsnLeftString  = "day "+iday(it)
     ResC@gsnRightString = expString


     ResC@mpMinLonF      = MinLon(it) 
     ResC@mpMaxLonF      = MaxLon(it) 
     ResC@mpCenterLonF   = CtrLon(it) 

     ResC@cnMinLevelValF  = varMin(it) 
     ResC@cnMaxLevelValF  = varMax(it) 
     ResC@cnLevelSpacingF = varInt(it) 

     if (VarName.eq."VOR") then

        if (it.eq.0) then
           delete(ResC@gsnSpreadColors)
           ResC@cnFillColors = (/36,40,44,48, 50,54,59,66,78,90/)
        else
           delete(ResC@cnFillColors)
           ResC@cnFillColors = (/36,43,48, 50,54,59,65, 75,89,92, 97/)
        end if
     end if

     plot(it) = gsn_csm_contour_map(wks,var,ResC)

  end do

 ;drawNDCGrid(wks)
  frame(wks)

end

