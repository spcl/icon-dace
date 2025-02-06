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
; iconplot.maps.ncl
;
; Versions:
; Initial version by Kristina Froehlich (10-2010)
; Update by Martin Koehler (11-2010, 2-2011)
;
; Concepts illustrated:
;   - Plotting ICON model data
;   - Contouring one-dimensional X, Y, Z data
;   - Using triangular meshes to create contours
;   - overlaying the refined grid onto the coarser one
;
;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "~kfroehli/NCL/lm_user.ncl"
;ATTENTION: colormaps:  /e/uhome/dreinert/NCL/ncl_5.2.1/lib/ncarg/colormaps
;---------------------------------------------------------------

begin

 testcase   = "APE"
 exp1       = "exp12"       ;experiment
 exp2       = "exp12"       ;control experiment for comparison
 Resolution1= "R2B04L31"
 Model      = "ICOHDC"      ;ICONAM or ICOHDC
 Model1     = "ICOHAM: " + exp1 
          ; +      "~C~" \
          ; + " diffusion + convection + cloud      ~C~" \
          ; + " microphysics + new satad~C~" \
          ; + " radiation"
 Model2     = "ICOHAM: " + exp1 + " - " + exp2 
          ; + "~C~" \
          ; + " diffusion + convection + cloud                ~C~" \
          ; + " microphysics + new satad          ~C~" \
          ; + " radiation"
 day        = "0040"        ;file number (0001, 0002, ...)
 timestep   = 20-1          ;0 = first step !!
 dtout      = 12            ;output frequency in [h]
 timebefore = 390*24        ;hours before current file (for plot title time)

 Diri       = "/e/uwork/mkoehler/icon/experiments/"+exp1+"/"
 Diri2      = "/e/uwork/mkoehler/icon/experiments/"+exp2+"/"
 FileName   = testcase+"_"+Resolution1+"_"+day+".nc"
 Diro       = Diri+"plots/"
 PlotFile   = Diro+testcase+"."+Resolution1+"."+exp1+"-"+exp2+".maps"
 result     = systemfunc("mkdir -p "+Diro) 


;---------------------------------------------------------------
; read in the meteorological field and grid information
;---------------------------------------------------------------

 File1 = addfile(Diri +FileName, "r" )
 File2 = addfile(Diri2+FileName, "r" )
 
 rad2deg = 45./atan(1.)        ; radians to degrees 
;coordinates of global coarse grid
 x = File1->clon *rad2deg      ; cell center, lon
 y = File1->clat *rad2deg      ; cell center, lat
;coordinates of refined grid  
;x_rf = File2->clon *rad2deg   ; cell center, lon
;y_rf = File2->clat *rad2deg   ; cell center, lat
 
 LeftString   = Resolution1
 RightString  = testcase+"~C~"+Resolution1
 FontHeight   = 0.02

 wks = gsn_open_wks("ps",PlotFile)


;---------------------------------------------------------------
; defaults
;---------------------------------------------------------------

;Set up resources for contour/map plot
 ResC                      = True

 ResC@mpProjection         = "CylindricalEquidistant"
;ResC@mpProjection         = "Satellite"
;ResC@mpPerimOn            = False        ; turn off box around plot
 ResC@mpCenterLonF         = 180.         ; choose center lon
;ResC@mpCenterLatF         = 20.          ; choose center lat
 ResC@mpLimitMode          = "LatLon"     ; required 
 ResC@mpFillOn             = False
;ResC@mpCenterLonF         = 180.
 ResC@mpMinLatF            = -90.
 ResC@mpMaxLatF            = 90.
;ResC@mpGreatCircleLinesOn = True

;ResC@gsnPolar             = "NH"  
 ResC@gsnDraw              = False
 ResC@gsnFrame             = False
 ResC@gsnMaximize          = False
 ResC@gsnSpreadColors      = True
 ResC@gsnStringFontHeightF = FontHeight*0.8
 ResC@gsnLeftString        = LeftString 
 ResC@gsnRightString       = RightString  
 ResC@cnFillOn             = True
 ResC@cnLinesOn            = False
 ResC@cnInfoLabelOn        = False
 ResC@cnFillMode           = "Areafill"
 ResC@cnRasterSmoothingOn  = True

 ResC@sfXArray             = x  ; These are 1D arrays, so a triangular
 ResC@sfYArray             = y  ; mesh will be created internally.

 ResC@tiMainFont           = 4
 ResC@tiMainFontHeightF    = FontHeight   
 ResC@tiMainFontThicknessF = 4
;ResC@tiMainOffsetYF       = 0.02
 ResC@tiXAxisFontHeightF   = 0.015;FontHeight
 ResC@tiXAxisOffsetYF      = 0.01
   
 ResC@pmLabelBarHeightF        = 0.07
 ResC@pmLabelBarWidthF         = 0.7
 ResC@pmLabelBarOrthogonalPosF = 0.1
 ResC@lbLabelBarOn             = True 
 ResC@lbLabelAutoStride        = True
 ResC@lbTitlePosition          = "Bottom"
 ResC@lbTitleFontHeightF       = 0.015
 
 ResP                            = True
 ResP@gsnMaximize                = True
   ;ResP@gsnPaperMargin             = 0.2
 ResP@gsnPanelXWhiteSpacePercent = 5
 ResP@gsnPanelYWhiteSpacePercent = 5
 ResP@txFontHeightF              = 0.01
   ;ResP@gsnPanelBottom             = 0.0   ;prevent landscape

 Res1                   = True
 Res1@gsnFrame          = False
 Res1@gsnPanelBottom    = 0.5
 Res1@gsnMaximize       = True
 Res1@gsnPanelYWhiteSpacePercent = 5

 Res2                   = True
 Res2@gsnFrame          = False
 Res2@gsnPanelTop       = 0.5
 Res2@gsnMaximize       = True
 Res2@gsnPanelYWhiteSpacePercent = 5


;---------------------------------------------------------------
; 3D variables: labels and variable dependent settings
;---------------------------------------------------------------

;VarNames = (/"OMEGA"/)   ; only for hydrostatic model
 VarNames = (/"T", "ACLC", "Qv", "Qw", "Qi", "Q4", "Q5", \
              "DIV", "U", "V", "OMEGA"/)
 VarNameLong = (/"Temperature [K]"            , \
                 "Cloud Cover [%]"            , \
                 "QV diagnostic [g/kg]"       , \
                 "QC diagnostic [g/kg]"       , \
                 "QI diagnostic [g/kg]"       , \
                 "QV grid-scale [g/kg]"       , \
                 "QC grid-scale [g/kg]"       , \
                 "QI grid-scale [g/kg]"       , \
                 "QR grid-scale [g/kg]"       , \
                 "QS grid-scale [g/kg]"       , \
                 "Divergence [1/h]"           , \
                 "U-wind (ms~S~-1~N~)"        , \
                 "V-wind (ms~S~-1~N~)"        , \
                 "OMEGA (Pas~S~-1~N~)"         /)
;Levels   = (/80, 69, 52, 38, 11/)  ; 1100m and 10100m on 60 levels
 Levels   = (/30,20,10,5/)  
 greyscale = (/"white" ,"black" ,"gray100","gray95","gray90","gray85","gray80","gray75","gray70", \
               "gray65","gray60","gray55" ,"gray50","gray45","gray40","gray35","gray30","gray25","gray20","gray15"/)

 do i = 0, dimsizes(VarNames)-1
  do k = 0, dimsizes(Levels)-1

   if (Model .eq. "ICOHDC") then
    ;apzero    = 101325.  ; ECHAM5 default (ncar 2008: apzero= 1.e5)
     apzero    = 1.e5
     vertc     = (File1->hyam)/apzero + File1->hybm  ; eta as vertical coordinate
     vertc_t   = 0.
     vertc_sfc = 1.
     height    = vertc
   end if
   if (Model .eq. "ICONAM") then
     if (VarNames(i) .eq. "W") then
      ;height  = (File1->hyai) bug:file doesn't have hyai  ; half levels (assuming orography is 0)
       height  = hyai                         ; half levels (assuming orography is 0)
     else
       hyai    = (File1->hyai)
       hyam    = (File1->hyam)                ; full levels (assuming orography is 0)
       height  = hyam                         ; full levels (assuming orography is 0)
     end if
   end if
 
   print("Varname="+VarNames(i)+"   Level="+Levels(k)+"   Height="+height(Levels(k)-1))
   var  = File1->$VarNames(i)$(timestep,Levels(k)-1,:)         ; dims: (time,lev,cell)
   var2 = File2->$VarNames(i)$(timestep,Levels(k)-1,:)         ; dims: (time,lev,cell)
  ;printVarSummary(var)

   if (VarNames(i).eq."T") then  
     scale      = 1.
     colormap   = "rainbow"
     varMin     = (/ 220./)
     varMax     = (/ 330./)
     varInt     = (/   3./)
   end if

   if ( VarNames(i) .eq."ACLC") then
     scale      = 100.0
     colormap   = greyscale
     varMin     = (/  10./)
     varMax     = (/ 100./)
     varInt     = (/  10./)
   end if 

   if ( VarNames(i) .eq."Qv") then
     scale      = 1000.0
    ;colormap   = "WhiteBlueGreenYellowRed"
     colormap   = "WhBlGrYeRe"
     varMin     = (/  2./)
     varMax     = (/ 22./)
     varInt     = (/  2./)
   end if 

   if ( VarNames(i) .eq."Qw") then
     scale      = 1000.0
     colormap   = "WhBlGrYeRe"
     varMin     = (/ 0./)
     varMax     = (/ 2./)
     varInt     = (/ 0.1/)
   end if

   if (VarNames(i) .eq."Qi") then
     scale      = 1000.0
     colormap   = "WhBlGrYeRe"
     varMin     = (/ 0./)
     varMax     = (/ 0.05/)
     varInt     = (/ 0.0025/)
   end if

   if ( VarNames(i) .eq."Q1") then
     scale      = 1000.0
     colormap   = "WhBlGrYeRe"
     varMin     = (/  2./)
     varMax     = (/ 22./)
     varInt     = (/  2./)
   end if

   if ( VarNames(i) .eq."Q2") then
     scale      = 1000.0
     colormap   = "WhBlGrYeRe"
     varMin     = (/ 0./)
     varMax     = (/ 2./)
     varInt     = (/ 0.1/)
   end if

   if (VarNames(i) .eq."Q3") then
     scale      = 1000.0
     colormap   = "WhBlGrYeRe"
     varMin     = (/ 0./)
     varMax     = (/ 0.05/)
     varInt     = (/ 0.0025/)
   end if

   if (VarNames(i) .eq."Q4") then
     scale      = 1000.0
     colormap   = "WhBlGrYeRe"
     varMin     = (/ 0./)
     varMax     = (/ 0.05/)
     varInt     = (/ 0.0025/)
   end if

   if (VarNames(i) .eq."Q5") then
     scale      = 1000.0
     colormap   = "WhBlGrYeRe"
     varMin     = (/ 0./)
     varMax     = (/ 0.05/)
     varInt     = (/ 0.0025/)
   end if

   if (VarNames(i) .eq."DIV") then
     scale      = 1.e3 * 3.6 
     colormap   = "ViBlGrWhYeOrRe"
     varMin     = (/-0.2/)
     varMax     = (/ 0.2/)
     varInt     = (/ 0.02/)
     ResS       = ResC
     symMinMaxPlt (var,20,False,ResS)
   end if

   if (VarNames(i) .eq."OMEGA") then
     scale      = 1.
     var@units  = "Pa/s"
     colormap   = "ViBlGrWhYeOrRe"
     varMin     = (/-2.  /)
     varMax     = (/ 2.  /)
     varInt     = (/ 0.02/)
     ResS       = ResC
     symMinMaxPlt (var,20,False,ResS)
   end if

   if (VarNames(i).eq."U") then
     scale      = 1.
     colormap   = "ViBlGrWhYeOrRe"
     varMin     = (/ -40./)
     varMax     = (/  40./)
     varInt     = (/   3./)
     ResS       = ResC
     symMinMaxPlt (var,20,False,ResS)
   end if

   if (VarNames(i).eq."V") then
     scale      = 1.
     colormap   = "ViBlGrWhYeOrRe"
     varMin     = (/ -40./)
     varMax     = (/  40./)
     varInt     = (/   3./)
     ResS       = ResC
     symMinMaxPlt (var,20,False,ResS)
   end if

   if (VarNames(i).eq."W") then
     scale      = 1.
     colormap   = "ViBlGrWhYeOrRe"
     varMin     = (/-0.05 /)
     varMax     = (/ 0.05 /)
     varInt     = (/ 0.005/)
     ResS       = ResC
     symMinMaxPlt (var,20,False,ResS)
   end if 

   var@long_name=VarNameLong(i)
  ;var@units    = VarUnit(i)

;---------------------------------------------------------------
; 3D variables: plot
;---------------------------------------------------------------

   var     = var *scale
   var2    = var2*scale
   vardiff = var - var2

   dt_forecast = timestep*dtout + timebefore

   TitleString = var@long_name + "  [L" + Levels(k) + " = " + height(Levels(k)-1) + "m]"

   ResC@cnMinLevelValF  = varMin
   ResC@cnMaxLevelValF  = varMax 
   ResC@cnLevelSpacingF = varInt 
   ResC@gsnLeftString   = dt_forecast/24.0 + " days"
   ResC@tiMainString    = TitleString
  ;ResC@lbTitleString   = var@long_name+" in "+var@units

   ResS = ResC

   plot = new( 2, graphic )

   gsn_define_colormap(wks,colormap)
   ResS@tiXAxisString   = "Min "+min(var)+"   Mean "+avg(var)+"   Maximum "+max(var)
   ResS@gsnCenterString = Model1
   plot(0) = gsn_csm_contour_map(wks,var,ResS)
   gsn_panel (wks, plot(0), (/1,1/), Res1)

   gsn_define_colormap(wks,"ViBlGrWhYeOrRe") 
   ResS@tiXAxisString   = "Min "+min(vardiff)+"   Mean "+avg(vardiff)+"   Maximum "+max(vardiff)
   ResS@gsnCenterString = Model2
   symMinMaxPlt (vardiff,20,False,ResS)
   plot(1) = gsn_csm_contour_map(wks,vardiff,ResS) 
   gsn_panel (wks, plot(1), (/1,1/), Res2)

   frame(wks)

   delete(height)
   delete(colormap)
   delete(ResS)

  end do
 end do





;---------------------------------------------------------------
; surface variables: plot
;---------------------------------------------------------------

 VarNames2D = (/ \
   "PS"   , "SKT"  , "PHIS" , "cosmu0" , "flxdwswtoa" , "APRL"  , "APRC" , \
   "APRS" , "RSFL" , "RSFC" , "SSFL"   , "SSFC"       , "ACLCOV" /)
;  "PS"          , "T_G"         , "QV_S"        ,             \
;  "flxdwswtoa"  , "swflxsfc"    , "lwflxsfc"    , "cosmu0"  , \
;  "TOT_PREC"    , "RAIN_GSP"    , \
;  "SNOW_GSP"    , "RAIN_CON"    , "SNOW_CON"    , \
;  "PRR_GSP"     , "PRS_GSP"     , \
;  "swflxsfc_avg", "swflxtoa_avg", "lwflxsfc_avg", "lwflxtoa_avg", \
;  "swflxsfc"    , "flxdwswtoa"  , "swflxtoa"    , "lwflxsfc"    , "lwflxtoa" /)

 do i = 0,dimsizes(VarNames2D) -1

   var  = File1->$VarNames2D(i)$(timestep,:)   ; dims: (time,cell)
   var2 = File2->$VarNames2D(i)$(timestep,:)   ; dims: (time,cell)

   if (VarNames2D(i) .eq."PS") then
     scale      = 1.e-2    
     var@units  = "hPa"
     colormap   = "ViBlGrWhYeOrRe"
     varMin     = (/ 940./)
     varMax     = (/ 1040./)
     varInt     = (/  5./)
   end if

   if (VarNames2D(i).eq."T_G" .or. VarNames2D(i).eq."SKT") then
     scale      = 1.
     var@long_name = "Surface (skin) Temperature"
     var@units  = "K"
     colormap   = "rainbow"
     varMin     = (/ 270./)
     varMax     = (/ 300./)
     varInt     = (/   2./)
   end if

   if (VarNames2D(i).eq."QV_S") then
     scale      = 1000.
     var@long_name = "Surface Specific Humidity"
     var@units  = "g/kg"
     colormap   = "rainbow"
     varMin     = (/ 25./)
     varMax     = (/  1./)
     varInt     = (/  2./)
   end if

   if (VarNames2D(i).eq."flxdwswtoa"   .or. \
       VarNames2D(i).eq."swflxsfc"     .or. VarNames2D(i).eq."swflxtoa") then
     scale      = 1.
     var@units  = "Wm2"
     colormap   = "prcp_2"
     varMin     = (/-100./)
     varMax     = (/ 100./)
     varInt     = (/ 50./)
   end if

   if (VarNames2D(i).eq."lwflxtoa"     .or. VarNames2D(i).eq."lwflxsfc") then
     scale      = 1.
     var@units  = "Wm2"
     colormap   = "prcp_2"
     varMin     = (/-140./)
     varMax     = (/  20./)
     varInt     = (/  20./)
   end if

   if (VarNames2D(i).eq."swflxsfc_avg" .or. VarNames2D(i).eq."swflxtoa_avg") then
     scale      = 1.
     var@units  = "Wm2"
     colormap   = "prcp_2"
     varMin     = (/   0./)
     varMax     = (/ 260./)
     varInt     = (/  20./)
   end if

   if (VarNames2D(i).eq."lwflxsfc_avg") then
     scale      = 1.
     var@units  = "Wm2"
     colormap   = "prcp_2"
     varMin     = (/-60./)
     varMax     = (/  0./)
     varInt     = (/  5./)
   end if

   if (VarNames2D(i).eq."lwflxtoa_avg") then
     scale      = 1.
     var@units  = "Wm2"
     colormap   = "prcp_2"
     varMin     = (/-300./)
     varMax     = (/-100./)
     varInt     = (/  20./)
   end if

   if (VarNames2D(i).eq."cosmu0") then
     scale      = 1.
     var@units  = "1"
     colormap   = "prcp_2"
     varMin     = (/ 0./)
     varMax     = (/ 1./)
     varInt     = (/ 0.1/)
   end if

   if (VarNames2D(i).eq."TOT_PREC") then
     scale      = 1./(dt_forecast/24.0)   ;/day
     var@long_name = "Total Precipitation (GS+conv)"
     var@units  = "mm/day"
     colormap   = "WhBlGrYeRe"
     varMin     = (/  2./)
     varMax     = (/ 40./)
     varInt     = (/  2./)
   end if

   if (VarNames2D(i).eq."RAIN_GSP" .or. VarNames2D(i).eq."SNOW_GSP" .or. \
       VarNames2D(i).eq."RAIN_CON" .or. VarNames2D(i).eq."SNOW_CON") then
     scale      = 1./(dt_forecast/24.0)   ; total to /day
     var@units  = "mm/day"
     colormap   = "WhBlGrYeRe"
     varMin     = (/  2./)
     varMax     = (/ 40./)
     varInt     = (/  2./)
   end if

   if (VarNames2D(i).eq."PRR_GSP"  .or. VarNames2D(i).eq."PRS_GSP") then
     scale      = 86400.0                 ;/s to /day
     var@units  = "mm/day"
     colormap   = "WhBlGrYeRe"
     varMin     = (/  2./)
     varMax     = (/ 40./)
     varInt     = (/  2./)
   end if

   print("Varname="+VarNames2D(i))

   var    = var *scale
   var2   = var2*scale
   vardiff= var-var2

   ResC@cnMinLevelValF  = varMin
   ResC@cnMaxLevelValF  = varMax 
   ResC@cnLevelSpacingF = varInt
   ResC@tiMainString    = var@long_name+" ["+var@units+"]"
  ;ResC@lbTitleString   = var@long_name+" in "+var@units

   ResS = ResC

   plot = new( 2, graphic )

   gsn_define_colormap(wks,colormap)
   ResS@tiXAxisString   = "Min "+min(var)+"   Mean "+avg(var)+"   Maximum "+max(var)
   ResS@gsnCenterString = Model1
   plot(0) = gsn_csm_contour_map(wks,var,ResS)
   gsn_panel (wks, plot(0), (/1,1/), Res1)
 
   gsn_define_colormap(wks,"ViBlGrWhYeOrRe") 
   ResS@tiXAxisString   = "Min "+min(vardiff)+"   Mean "+avg(vardiff)+"   Maximum "+max(vardiff)
   ResS@gsnCenterString = Model2
   symMinMaxPlt (vardiff,20,False,ResS)
   plot(1) = gsn_csm_contour_map(wks,vardiff,ResS)
   gsn_panel (wks, plot(1), (/1,1/), Res2)

   frame(wks)

   delete(ResS)

 end do


 system("scp " + PlotFile + ".ps oflxs04:/fe1-daten/mkoehler/plots/icon")
;system("ps2pdf " + PlotFile + ".ps")
;system("gv "+PlotFile +" &")

end
