; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;-----------------------------------------------------------------------------
; Hui Wan (MPI, 2005.07.27)
; Martin Koehler ( DWD, Nov 2010)
;-----------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;-----------------------------------------------------------------------------

begin

 PlotFormat = getenv("plot_file_format")
 Model      = getenv("Model")
 DataPath   = getenv("DataPath")
 DataID     = getenv("DataID")
 DataGrid   = getenv("trunc")
 PlotPath   = getenv("PlotPath")
 ConfigStr  = getenv("ConfigStr")
 Resolution = getenv("Resolution")
 evol_stat  = getenv("evol_stat")
 top_title  = getenv("top_title")
 expnum     = getenv("expnum")
 EXPNAME    = getenv("EXPNAME")
 dtout1     = getenv("dtout")
 dtout      = stringtofloat(dtout1)

 PlotFile = PlotPath+EXPNAME+"."+Resolution+"."+expnum+".timeseries"


;--------------------------------------------------------------- 
; setup (synchronize with batch file iconplot.zonal.s)
;--------------------------------------------------------------- 

 VarNames2 = (/"Cloud Cover [1]"            , \
               "Temperature [K]"            , \
               "QV grid-scale [g/kg]"       , \
               "QC grid-scale [mg/kg]"      , \
               "QI grid-scale [mg/kg]"      , \
               "QR grid-scale [mg/kg]"      , \
               "QS grid-scale [mg/kg]"      , \
               "U-wind (ms~S~-1~N~)"        , \
               "V-wind (ms~S~-1~N~)"        , \
               "OMEGA-wind (Pas~S~-1~N~)"    /)
 
 VarNames = (/"ACLC", "T", "Qv", "Qw", "Qi", "Q4", "Q5", "U", "V", "OMEGA"/)
 Scale =     (/1e2,    1.,  1e3,  1e6,  1e6,  1e6,  1e6,  1.,  1.,  1./)
;Levels   = (/80,69,52,38,11/)  ; 1000m, 5000m, 10000m, and 30000m on 80 levels
 Levels   = (/30,20,10,5/)  ; 1000m, 5000m, 10000m, and 30000m on 80 levels


;--------------------------------------------------------------- 
; plot settings
;--------------------------------------------------------------- 

 ResC = True

 ResC@gsnMaximize    = True
 ResC@vpHeightF      = 0.2
 ResC@vpWidthF       = 0.6
 ResC@vpXF           = 0.03
 ResC@vpYF           = 0.97

 ResC@tiYAxisString          = "Latitude"
 ResC@tiYAxisFontHeightF     = 0.016
 ResC@tmYLLabelFontHeightF   = 0.014
 ResC@tmYLMajorLengthF       = 0.01
 ResC@tmYLMinorLengthF       = 0.005
 ResC@tmYRMajorLengthF       = 0.01
 ResC@tmYRMinorLengthF       = 0.005
 ResC@tmYLMajorThicknessF    = 3. 
 ResC@tmYRMajorThicknessF    = 3.

 ResC@tmXBLabelFontHeightF   = 0.014
 ResC@tmXBMajorLengthF       = 0.01
 ResC@tmXBMinorLengthF       = 0.005
 ResC@tmXBMajorThicknessF    = 3. 
 ResC@tmXTMajorThicknessF    = 3. 

 ResC@cnFillOn               = True
 ResC@cnLinesOn              = False
 ResC@cnInfoLabelOn          = False

 ResC@lbLabelBarOn             = True
 ResC@lbOrientation            = "vertical"
;ResC@pmLabelBarOrthogonalPosF = -0.001
 ResC@pmLabelBarWidthF         = 0.07
 ResC@lbLabelFontHeightF       = .012
 ResC@lbLabelAutoStride        = True
 ResC@gsnSpreadColors          = True            ; use full range of colors
 ResC@gsnSpreadColorStart      = 3
 ResC@gsnSpreadColorEnd        = 101


;--------------------------------------------------------------- 
; plot
;--------------------------------------------------------------- 

 wks_type = PlotFormat        ; "ps" or "pdf"
 wks_type@wkOrientation = "landscape"
 wks = gsn_open_wks( wks_type, PlotFile )

;gsn_define_colormap( wks, "BkBlAqGrYeOrReViWh200" )
;gsn_define_colormap( wks, "testcmap" )
 gsn_merge_colormaps( wks, "WhBlGrYeRe", "testcmapshort")

 do i = 0,dimsizes(VarNames) -1
 do k = 0,dimsizes(Levels) -1

   DataFile = DataPath+DataID+"."+VarNames(i)+"_L"+Levels(k)+".zonvar.nc"
   File     = addfile( DataFile,"r" )
   var      = File->$VarNames(i)$(:,0,:,0)  ; (time,lev,lat,lon)
   var      = var * Scale(i)

   if (Model .eq. "ICOHDC") then
    ;apzero     = 101325.  ; ECHAM5 default (ncar 2008: apzero= 1.e5)
     apzero     = 1.e5
     vertc     = (File->hyam)/apzero + File->hybm  ; eta as vertical coordinate
     vertc_t   = 0.
     vertc_sfc = 1.
     height    = vertc
   end if
   if (Model .eq. "ICONAM") then
     if (VarNames(i) .eq. "W") then
      ;height  = (File->hyai) bug:file doesn't have hyai  ; half levels (assuming orography is 0)
       height  = hyai                         ; half levels (assuming orography is 0)
     else
       hyai    = (File->hyai)
       hyam    = (File->hyam)                 ; full levels (assuming orography is 0)
       height  = hyam                         ; full levels (assuming orography is 0)
     end if
   end if

   ntime    = dimsizes(var&time)
   delete(var&time)
  ;var&time = fspan(0,ntime-1,ntime)
   var&time = fspan(0,(ntime-1)*dtout/24.0,ntime)

   ResC@tiMainString        = VarNames2(i)+ "      (L" + Levels(k) + "  sigma =  " + height(Levels(k)-1) + ")"
   ResC@gsnLeftString       = evol_stat + "        mean: " + avg(var)
  ;ResC@gsnCenterString     = "Min-Mean-Max: "+min(var)+"   "+avg(var)+"   "+max(var)
   ResC@gsnCenterString     = top_title + "~C~~C~  "
   ResC@gsnRightString      = Resolution
   ResC@tiXAxisJust         = "CenterCenter"
   ResC@tiXAxisString       = "Time (days)"
   ResC@tiXAxisFontHeightF  = 0.016
   delete(height)
 
   contour = gsn_csm_lat_time( wks, var(lat|:,time|:), ResC )

 end do
 end do

 system("scp " + PlotFile + ".ps oflxs04:/fe1-daten/mkoehler/plots/icon")

end
