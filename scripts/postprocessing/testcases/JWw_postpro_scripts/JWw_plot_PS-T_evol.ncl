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
; Hui Wan (MPI-Met, 2008-07-13)
; Hui Wan (MPI-Met, 2009-04-28: Modification for the hex. model)
;------------------------------------
; Script type: visualization
;---------------------------------------------------------------
; This script makes contour plots of 850 hPa temperature 
; for the Jablonowski-Williamson steady state test.
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

  rad2deg = 45./atan(1.)   ; radians to degrees

;---------------------------------------------------------------
; names of the data file and the plot
;---------------------------------------------------------------

  if (VarName.eq."T") then
     DataFileName = DataPath + ExpName +"_T850.nc"
     PlotFileName = PlotPath + plotBaseName +"_T850"
  else
     DataFileName = DataPath + ExpName +".nc"
     PlotFileName = PlotPath + plotBaseName +"_PS"
  end if

;---------------------------------------------------------------
; time steps to be plotted
;---------------------------------------------------------------

  iday  = (/4,6,7,8,10,  9/)  ; days to plot
  ndays = dimsizes(iday)
  plot1 = new( ndays-1,graphic )

; time indices in the model output

  if (Model.eq."ECHAM5") then
     tid = OutFrq*iday -1  ; ECHAM doesn't write out initial data
  else
     tid = OutFrq*iday 
  end if

;---------------------------------------------------------------
; labels
;---------------------------------------------------------------

  LeftString   = new( ndays,string )
  LeftString   = ""
  RightString  = LeftString
  CenterString = LeftString

  if (VarName.eq."T") then
     LeftString(0)       = "850 hPa T (K) at day "+ iday(0) 
     LeftString(ndays-1) = "850 hPa T (K) at day "+ iday(ndays-1) 
  else
     LeftString(0)       = "ps (hPa) at day "+ iday(0)
     LeftString(ndays-1) = "ps (hPa) at day "+ iday(ndays-1)
  end if

  CenterString          = "day "+iday
  CenterString(0)       = ""
  CenterString(ndays-1) = ""

  RightString(0)        = ConfigStr   ;+" "+Resolution 
  RightString(ndays-1)  = ConfigStr   ;+" "+Resolution 

;---------------------------------------------------------------
; lat./lon. of cell center (ICOHDC only)
;---------------------------------------------------------------

  if (Model.eq."ICOHDC") then
     if (CellType.eq."3") then
        latName    = "clat"
        lonName    = "clon"
     else
     if (CellType.eq."6") then
        latName    = "vlat"
        lonName    = "vlon"
     else
        print("Wrong cell type. Should be 3 or 6")
     end if
     end if
  end if

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
     print("Problem encountered: found only "+ ndays_in_data_file+ \
           " day(s) of data in the model output")
     print(nmax_needed +" days are needed for visualising "       \
           +"the evolution of temperature and surface pressure")  
     status_exit (1)
  end if

;---------------------------------------------------------------
; make plots
;---------------------------------------------------------------
; first the "evolution plot"

  wks = gsn_open_wks(PlotFormat,PlotFileName)
  gsn_define_colormap(wks,"JWw_"+VarName+"_colormap")

  ResC                = True
  ResC@gsnDraw        = False
  ResC@gsnFrame       = False

  ResC@tmXBLabelJust  = "CenterCenter"

  if (Model.eq."ICOHDC") then
     ResC@sfXArray = x
     ResC@sfYArray = y
  end if

  load "JWw_PS-T_evol_resources.ncl"

;---------------------------------------------------------------
; read data of each time step and make contour plot
;---------------------------------------------------------------

  do it = 0,ndays-2

     if (Model.eq."ICOHDC") then

        if (VarName.eq."T") then
            var = File->ta(tid(it),0,:)
        else
            var = File->ps(tid(it),:)
            var = var*0.01
        end if

     else ;ECHAM5
        if (VarName.eq."T") then
            var = File->st(tid(it),0,:,:)
        else
            var = File->ps(tid(it),:,:)
            var = var*0.01
        end if
     end if

     ResC@gsnLeftString   = LeftString(it)
     ResC@gsnCenterString = CenterString(it)
     ResC@gsnRightString  = RightString(it) 
     ResC@lbLabelBarOn    = False 

     if (VarName.eq."T") then
        ResC@cnMinLevelValF       = tempMin 
        ResC@cnMaxLevelValF       = tempMax 
        ResC@cnLevelSpacingF      = tempInt 
     else 

      if (iday(it).lt.7) then
        ResC@cnMinLevelValF       = psMin1 
        ResC@cnMaxLevelValF       = psMax1 
        ResC@cnLevelSpacingF      = psInt1 
        ResC@cnFillColors         = (/2,2,4,4,7,7,8,8,9,9,10,10,11,11/)
      else
        ResC@cnMinLevelValF       = psMin2 
        ResC@cnMaxLevelValF       = psMax2 
        ResC@cnLevelSpacingF      = psInt2 
       
        if (isatt(ResC,"cnFillColors")) then 
            delete(ResC@cnFillColors)
        end if
      end if

     end if


     if (iday(it).eq.6) then
        ResC@lbLabelBarOn  = True
     end if

     if (it.eq.ndays-2) then
        ResC@lbLabelBarOn  = True
     end if

     plot1(it) = gsn_csm_contour_map(wks,var,ResC)

  end do

;---------------------------------------------------------------
; put all panels on one page
;---------------------------------------------------------------

  ResP                            = True
  ResP@gsnPanelYWhiteSpacePercent = 9
  ResP@gsnFrame = False

  t1 = 1.
  b1 = 0.67

  ResP@gsnPanelTop    = t1
  ResP@gsnPanelBottom = b1

  gsn_panel (wks, plot1(0:1), (/2,1/), ResP)

  t2 = 0.62
  b2 = t2 - 1.5*(t1-b1)

  ResP@gsnPanelTop    = t2
  ResP@gsnPanelBottom = b2
  gsn_panel (wks, plot1(2: ), (/3,1/), ResP)

  frame(wks)
  delete(wks)

;---------------------------------------------------------------
; then the snapshot of day 9 in a sperate file
;---------------------------------------------------------------

  wks = gsn_open_wks(PlotFormat,PlotFileName+"_day09")
  gsn_define_colormap(wks,"JWw_"+VarName+"_colormap")

  it = ndays - 1

  if (Model.eq."ICOHDC") then

     if (VarName.eq."T") then
         var = File->ta(tid(it),0,:)
     else
         var = File->ps(tid(it),:)
         var = var*0.01
     end if

  else ;ECHAM5
     if (VarName.eq."T") then
         var = File->st(tid(it),0,:,:)
     else
         var = File->ps(tid(it),:,:)
         var = var*0.01
     end if
  end if

  ResC@gsnDraw        = True 
  ResC@gsnFrame       = True

  ResC@gsnLeftString   = LeftString(it)
  ResC@gsnCenterString = CenterString(it)
  ResC@gsnRightString  = RightString(it) 
  ResC@lbLabelBarOn    = True 

  plot2 = gsn_csm_contour_map(wks,var,ResC)
end

