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
; Hui Wan (MPI, 2005-07-27)
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

  title    = "Held-Suarez Test"
  DataFile = DataPath+DataID+".nc" 
  PlotFile = PlotPath+DataID

;----------------------------------------------


  File = addfile( DataFile,"r" )

;--------------------------------------------------------------- 
; plot settings

     ResXY = True

     ResXY@gsnMaximize            = True
     ResXY@gsnPaperOrientation    ="portrait"

     ResXY@tiYAxisString          = "Surface Pressure Global Mean (hPa)"
     ResXY@tiYAxisFontHeightF     = 0.016
     ResXY@trYMinF                =  999.75
     ResXY@trYMaxF                = 1000.25

     ResXY@tiXAxisString          = "Time ( days )"
     ResXY@tiXAxisFontHeightF     = 0.016



     wks = gsn_open_wks( PlotFormat,PlotFile )

     var = File->$VarName$(:,0,0)  ;(time,lat,lon)
     var = var/100.

     ntime = dimsizes(var&time)
     
     time = fspan(0,ntime-1,ntime)

     ;ResXY@gsnLeftString   = ConfigStr
     ;ResXY@gsnCenterString = ""
     ;ResXY@gsnRightString  = "Surface Pressure global mean" 
     ResXY@tiMainString = ConfigStr
     ResXY@tiMainFontHeightF  = 0.016

     contour = gsn_csm_xy( wks,time,var,ResXY )

end
