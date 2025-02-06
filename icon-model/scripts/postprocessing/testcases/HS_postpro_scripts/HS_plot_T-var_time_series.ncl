; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;
; Hui Wan (MPI-M, 2005-07-27)
; Marco Giorgetta (MPI-M, 2017-04-26)
;-----------------------------------------------------------------------------
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;-----------------------------------------------------------------------------
 begin
  PlotFormat = getenv("plot_file_format")
  DataPath   = getenv("DataPath")
  DataID     = getenv("DataID")
  PlotPath   = getenv("PlotPath")
  ConfigStr  = getenv("ConfigStr")
  VarName    = getenv("VarName")

  title = "Held-Suarez Test"
  DataFile = DataPath+DataID+".nc" 
  PlotFile = PlotPath+DataID

;----------------------------------------------


  File = addfile( DataFile,"r" )
  lev  = doubletointeger(File->lev)
  lev  = lev/100
;--------------------------------------------------------------- 
; plot settings

     ResC = True

     ResC@gsnMaximize    = True
     ResC@vpHeightF      = 0.2
     ResC@vpWidthF       = 0.6
     ResC@vpXF           = 0.03
     ResC@vpYF           = 0.97

  ;  ResC@tiYAxisString          = "Latitude"
  ;  ResC@tiYAxisFontHeightF     = 0.016
     ResC@tmYLLabelFontHeightF   = 0.014
     ResC@tmYLMajorLengthF       = 0.01
     ResC@tmYLMinorLengthF       = 0.005
     ResC@tmYRMajorLengthF       = 0.01
     ResC@tmYRMinorLengthF       = 0.005
     ResC@tmYLMajorThicknessF    = 3. 
     ResC@tmYRMajorThicknessF    = 3.

     ResC@tiXAxisString          = "Time ( days )"
     ResC@tiXAxisFontHeightF     = 0.016

     ResC@tmXBLabelFontHeightF   = 0.014
     ResC@tmXBMajorLengthF       = 0.01
     ResC@tmXBMinorLengthF       = 0.005
     ResC@tmXBMajorThicknessF    = 3. 
     ResC@tmXTMajorThicknessF    = 3. 

     ResC@cnFillOn       = True
     ResC@cnLinesOn      = False
     ResC@cnInfoLabelOn  = False

     ResC@lbLabelBarOn             = True
     ResC@lbOrientation            = "vertical"
  ;  ResC@pmLabelBarOrthogonalPosF = -0.001
     ResC@pmLabelBarWidthF         = 0.07
     ResC@lbLabelFontHeightF       = .012
     ResC@lbLabelAutoStride        = True
     ResC@gsnSpreadColors          = True            ; use full range of colors
     ResC@gsnSpreadColorStart      = 2
     ResC@gsnSpreadColorEnd        = 176

     wks = gsn_open_wks( PlotFormat,PlotFile )
     gsn_define_colormap( wks,"BkBlAqGrYeOrReViWh200" )

     var = File->$VarName$(:,0,:,0)  ;(time,lev,lat,lon)
     ntime = dimsizes(var&time)
     delete(var&time)
     var&time = fspan(0,ntime-1,ntime)

     ResC@gsnLeftString   = ConfigStr
     ResC@gsnCenterString = ""
     ResC@gsnRightString  = "Zonal temperature variance (K~S~2~N~) at "+lev+" hPa"
     
     ResC@cnLevelSelectionMode  = "ManualLevels"
     ResC@cnMinLevelValF        =  0             
     ResC@cnMaxLevelValF        =  80 
     ResC@cnLevelSpacingF       =  5 

     contour = gsn_csm_lat_time( wks,var(lat|:,time|:),ResC )

end
