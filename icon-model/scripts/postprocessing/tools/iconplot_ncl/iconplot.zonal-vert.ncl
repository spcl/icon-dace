; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;---------------------------------------------------------
; Plot zonal mean cross-sections.
;
; Hui Wan ( MPI, Aug 2006 )
; Hui Wan ( MPI, Nov 2008 )
; Pilar Ripodas  ( DWD, Oct 2010)
; Martin Koehler ( DWD, Nov 2010)
;---------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;---------------------------------------------------------

begin

;---------------------------------------------------------
; input, variable setup
;---------------------------------------------------------

 PlotFormat = getenv("plot_file_format") 
 Model      = getenv("Model")
 DataPath   = getenv("DataPath")
 DataPath2  = getenv("DataPath2")
 DataID     = getenv("DataID")
 DataGrid   = getenv("trunc")
 PlotPath   = getenv("PlotPath")
 ConfigStr  = getenv("ConfigStr")
 Resolution = getenv("Resolution")
 CellType   = getenv("cell_type")
 time_stat_text = getenv("time_stat_text")
 top_title  = getenv("top_title")
 expnum     = getenv("expnum")
 expnum2    = getenv("expnum2")
 EXPNAME    = getenv("EXPNAME")

 PlotFile = PlotPath+EXPNAME+"."+Resolution+"."+expnum+"-vs-"+expnum2+".zonal_cross"
 suffix   = ".T"+DataGrid+".nc"

;zmta: zonal mean, time average (or time step or time variance)

 varNameF   = (/"CC.zmta"  ,"T.zmta"   ,"QV.zmta"     ,"Q1.zmta"     , \
                "QC.zmta"  ,"Q2.zmta"  ,"QI.zmta"     ,"Q3.zmta"     , \
                "Q4.zmta"  ,"Q5.zmta"  ,                               \
                "U.zmta"   ,"V.zmta"   ,"W.zmta"      ,                \
                "UpVp.zmta","TpVp.zmta","eddy-KE.zmta","TpTp.zmta"    /) 
 varName    = (/"CC"       ,"T"        ,"QV"          ,"Q1"          , \
                "QC"       ,"Q2"       ,"QI"          ,"Q3"          , \
                "Q4"       ,"Q5"       ,                               \
                "U"        ,"V"        ,"W"           ,                \
                "U"        ,"T"        ,"U"           ,"T"            /)

;apzero     = 101325.  ; ECHAM5 default (ncar 2008: apzero= 1.e5)
 apzero     = 1.e5

 VarNameLong = (/"Cloud Cover [1]"            , \                ;1
                 "Temperature [K]"            , \                ;2
                 "QV diagnostic [g/kg]"       , \                ;3
                 "QV grid-scale [g/kg]"       , \                ;4
                 "QC diagnostic [mg/kg]"      , \                ;5
                 "QC grid-scale [mg/kg]"      , \                ;6
                 "QI diagnostic [mg/kg]"      , \                ;7
                 "QI grid-scale [mg/kg]"      , \                ;8
                 "QR grid-scale [mg/kg]"      , \                ;9
                 "QS grid-scale [mg/kg]"      , \                ;10
                 "U-wind (ms~S~-1~N~)"        , \                ;11
                 "V-wind (ms~S~-1~N~)"        , \                ;12
                 "W-wind (ms~S~-1~N~)"        , \                ;13
                 "Eddy momentum flux (m~S~2~N~s~S~-2~N~)" , \    ;14
                 "Eddy heat flux (Kms~S~-1~N~)"           , \    ;15
                 "Eddy kinetic energy (m~S~2~N~s~S~-2~N~)", \    ;16
                 "Temperature variance (K~S~2~N~)"         /)    ;17
 nvar = dimsizes(VarNameLong)

;                 1,   2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,  9 , 10 , 11 , 12 , 13    , 14 , 15 , 16 , 17
 Scale =      (/1e2,   1., 1e3, 1e3, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6,  1.,  1.,  1.   ,   1,  1.,  1., 1./)
 Min =        (/  5, 160.,   1,   1,  10,  10,   2,   2,   2,   5, -20,  -5, -0.03 ,-200,-100, 100, 10/)
 Max =        (/ 80, 305.,  15,  15, 190, 190, 100, 100, 100, 200, 120,   5,  0.03 , 200, 100,2000,400/)
 Int =        (/ 10,    5,   2,   2,  20,  20,   5,   5,   5,  10,  10, 0.5,  0.003,  20,  10, 100, 20/)
 Max2=        (/ 30,   20, 0.5, 0.5,  50,  50,  20,  20,  20,  40,  40, 1.5,  0.005,  60,  20, 500,100/)
 Int2=        (/  3,    2, 0.1, 0.1,  10,  10,   2,   2,   2,   4,   4, 0.3,  0.001,   6,   4,  50, 10/)
					       		 		                
 colorStart = (/  0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,      0,   0,   0,   0,  0/)
 colorEnd   = (/101,  101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101,    101, 101, 101, 101,101/)

 nlev4plot  = 101
 iplot      = 0

 topheight  = 60.0     ; top in [km]

;---------------------------------------------------------
; plot setup
;---------------------------------------------------------

 wks_type = PlotFormat        ; "ps" or "pdf"
 wks_type@wkOrientation = "landscape"
 wks = gsn_open_wks( wks_type, PlotFile )

 contour = new( nvar*2, graphic )

 ResC = True

 ResC@gsnDraw  = False
 ResC@gsnFrame = False

 if (Model .eq."ICOHDC") then
   ResC@tiYAxisString          = "Hybrid coordinate ~F8~h~F~" 
 end if
 if (Model .eq."ICONAM") then
   ResC@tiYAxisString          = "Height Km" 
 end if     

 FontHeight = 0.02

 ResC@tiXAxisFontHeightF        = FontHeight
 ResC@tiYAxisFontHeightF        = FontHeight
 ResC@tmXBLabelFontHeightF      = FontHeight
 ResC@tmYLLabelFontHeightF      = FontHeight
 ResC@gsnLeftStringFontHeightF  = FontHeight + 0.003
 ResC@gsnCenterStringFontHeightF= FontHeight + 0.004
 ResC@gsnRightStringFontHeightF = FontHeight + 0.003

 ResC@lbLabelBarOn             = True
 ResC@pmLabelBarHeightF        = 0.1
 ResC@pmLabelBarWidthF         = 0.58
 ResC@lbLabelAutoStride        = True
 ResC@pmLabelBarOrthogonalPosF = 0.01
 ResC@lbLabelFontHeightF       = FontHeight 

;ResC@trYLog          = False
 ResC@tmYRMode        = "Automatic"
 if (Model .eq."ICOHDC") then     
   ResC@trYReverse    = True
   ResC@trYMaxF       = 1.
   ResC@trYMinF       = 0.
 end if
 if (Model .eq."ICONAM") then 
  ;ResC@trGridType    = "LogLin"
  ;ResC@trYAxisType   = "LogAxis"
  ;ResC@trYLog        = True    
   ResC@trYReverse    = False
   ResC@trYMaxF       = topheight
   ResC@trYMinF       = 0.
  ;ResC@trYLog        = True
  ;ResC@tmYLMinorPerMajor=8
 end if 

 ResC@cnInfoLabelOn                = False
;ResC@cnLineThicknessF             = 1.2
 ResC@cnLinesOn                    = True
 ResC@gsnContourZeroLineThicknessF = 1.5
 ResC@gsnContourNegLineDashPattern = 2
 ResC@cnFillOn                     = True
 ResC@gsnSpreadColors              = True
;ResC@cnLineLabelsOn               = True
;ResC@cnLineLabelPlacementMode     = "Computed"
;ResC@cnLineLabelInterval          = 4
;ResC@cnLineLabelAngleF            = 0.0
;ResC@cnLineLabelFontHeightF       = 0.016
;ResC@cnLineLabelFontColor         = "black"
;ResC@cnLabelMasking               = True
;ResC@cnLineLabelBackgroundColor   = "transparent"
;ResC@cnLineLabelFormat            = "0@*+^.2g"

 res_L                   = True
 res_L@gsnFrame          = False
 res_L@gsnPanelRight     = 0.5
 res_L@gsnMaximize       = True                       ; maximize plots
 res_L@txString          = top_title

 res_R                   = True
 res_R@gsnFrame          = False
 res_R@gsnPanelLeft      = 0.5
 res_R@txString          = top_title

 txres                   = True
 txres@txFontHeightF     = 0.02


;---------------------------------------------------------
; data reading and plotting
;---------------------------------------------------------

 do ivar=0,nvar-1
   File = addfile( DataPath  + DataID +"." + varNameF(ivar) +suffix, "r" )
   File2= addfile( DataPath2 + DataID +"." + varNameF(ivar) +suffix, "r" )
 
   var  = (/File->$varName(ivar)$(0,:,:,0)/)          ; (time,lev,cell,lon)
   var2 = (/File2->$varName(ivar)$(0,:,:,0)/)         ; (time,lev,cell,lon)
   var  = var * Scale(ivar)
   var2 = var2* Scale(ivar)
  
   do istat=0,1                                       ; 0: single experiment, 1:diff exp-control
     if (istat .eq. 1) then
       var = var - var2
     end if

     print( varName(ivar) + " min=" + min(var) + " max=" + max(var) )
 
     if (Model .eq. "ICOHDC") then
        vertc     = (File->hyam)/apzero + File->hybm  ; eta as vertical coordinate
        vertc_t   = 0.
        vertc_sfc = 1.
     end if
     if (Model .eq. "ICONAM") then
        if (varName(ivar) .eq. "W") then
          height  = (File->hyai)                      ; half levels (assuming orography is 0)
        else
          height  = (File->hyam)                      ; full levels (assuming orography is 0)
        end if
        vertc     =  height/1000.                     ; height in km as vertical coordinate
        vertc_t   = topheight
        vertc_sfc = 0.
     end if
 
     var!0     = "vertc"
     var&vertc = vertc
     var!1     = "lat"
     var&lat   = File->lat
 
     vertc4plot       = fspan(vertc_t,vertc_sfc,nlev4plot)
     vertc4plot!0     = "vertc"
     vertc4plot&vertc = vertc4plot
     tmp              = linint1_Wrap( vertc, var(lat|:,vertc|:), False, vertc4plot, 0 )
     var4plot         = tmp(vertc|:,lat|:)
 
     ResC@cnLevelSelectionMode  = "ManualLevels"
     if (istat .eq. 0) then
       ResC@cnMinLevelValF        =  Min(ivar)
       ResC@cnMaxLevelValF        =  Max(ivar)
       ResC@cnLevelSpacingF       =  Int(ivar)
      ;gsn_merge_colormaps( wks,"WhBlGrYeRe","testcmapshort")
      ;colormap                   = "testcmapshort"
       colormap                   = "WhBlGrYeRe"
       ResC@gsnSpreadColorStart   = colorStart(ivar)
       ResC@gsnSpreadColorEnd     = colorEnd(ivar)
       ResC@gsnCenterString       = ""
     else
       ResC@cnMinLevelValF        = -Max2(ivar)
       ResC@cnMaxLevelValF        =  Max2(ivar)
       ResC@cnLevelSpacingF       =  Int2(ivar)
      ;colormap   = "BlueGreen14"
      ;colormap   = "testcmap"
       colormap   = "ViBlGrWhYeOrRe"
       ResC@gsnSpreadColorStart   = 2
       ResC@gsnSpreadColorEnd     = 250
       ResC@gsnCenterString       = "vs. "+expnum2
     end if
     gsn_define_colormap(wks,colormap) 

     ResC@gsnLeftString         = time_stat_text
     ResC@gsnRightString        = Resolution
     ResC@tiMainString          = VarNameLong(ivar)
     ResC@tiXAxisString         = "Min-Mean-Max: "+min(var)+"   "+avg(var)+"   "+max(var)

     delete(height)
     delete(vertc)

     contour(iplot) = gsn_csm_contour( wks,var4plot,ResC ) 
     if (istat .eq. 0) then
       gsn_panel (wks, contour(iplot), (/1,1/), res_L)
     else
       gsn_panel (wks, contour(iplot), (/1,1/), res_R)
     end if
     iplot          = iplot + 1
   end do
   frame(wks)
   delete(var)                ; necessary because U/W on full/half levels
   delete(var2)               ; necessary because U/W on full/half levels
 end do


 system("scp " + PlotFile + ".ps oflxs04:/fe1-daten/mkoehler/plots/icon")

end
