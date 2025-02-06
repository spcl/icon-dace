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
; Hui Wan ( MPI-M, Aug 2006 )
; Hui Wan ( MPI-M, Nov 2008 )
; Pilar Ripodas ( DWD, Oct 2010)
; Marco Giorgetta (MPI-M, 2017-04-26)
;---------------------------
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;---------------------------------------------------------
 begin

    PlotFormat = getenv("plot_file_format") 
    DataPath   = getenv("DataPath")
    DataID     = getenv("DataID")
    DataGrid   = getenv("trunc")
    PlotPath   = getenv("PlotPath")
    ConfigStr  = getenv("ConfigStr")
    PlotTop    = getenv("PlotTop")
;------

    suffix     = "_t"+DataGrid+"grid.nc"

    varNameF   = (/"ta_ta","ta_tv","tapvap_ta","ua_ta","ua_tv","eke_ta","va_ta","va_tv","uapvap_ta","wa_ta","wa_tv","uapwap_ta"/) 
    varName    = (/"ta"   ,"ta"   ,"ta"       ,"ua"   ,"ua"   ,"ua"    ,"va"   ,"va"   ,"ua"       ,"wa"   ,"wa"   ,"ua"       /)

;------

 LeftString = (/"(a) temperature T (K)",                                                 \
                "(b) temperature variance <T'~S~2~N~> (K~S~2~N~)",                       \
                "(c) merid. eddy heat flux <T'v'> (Kms~S~-1~N~)",                        \
                "(d) zonal wind u (ms~S~-1~N~)",                                         \
                "(e) zonal wind variance <u'~S~2~N~> (m~S~2~N~s~S~-2~N~)",               \
                "(f) eddy kinetic energy <(u'~S~2~N~+v'~S~2~N~)/2> (m~S~2~N~s~S~-2~N~)", \
                "(g) merid. wind v (ms~S~-1~N~)",                                        \
                "(h) merid. wind variance <v'~S~2~N~> (m~S~2~N~s~S~-2~N~)",              \
                "(i) merid. eddy momentum flux <u'v'> (m~S~2~N~s~S~-2~N~)",              \
                "(j) vertical wind w (1e-2 ms~S~-1~N~)",                                 \
                "(k) vertical wind variance <w'~S~2~N~> (1e-4 m~S~2~N~s~S~-2~N~)",       \
                "(l) vertical eddy momentum flux <u'w'> (1e-2 m~S~2~N~s~S~-2~N~)"        /)

 nvar = dimsizes(LeftString)
 
 Min = (/ 190,   4, -24, -32,  40,  40, -0.80,  40, -64, -0.25, 0.4, -2.5 /)
 Max = (/ 305,  52,  24,  32, 480, 480,  0.80, 480,  64,  0.25, 4.0,  2.5 /)
 Int = (/   5,   4,   3,   4,  40,  40,  0.08,  40,   8,  0.05, 0.4,  0.5 /)

 colorStart = (/  2,   2, 102, 102,   2,   2, 102,   2, 102, 102,   2, 102/)
 colorEnd   = (/101, 101, 201, 201, 101, 101, 201, 101, 201, 201, 101, 201/)

 nlev4plot = 101

;---------------------------------------------------------
; make plots
;---------------------------------------------------------

    wks = gsn_open_wks( PlotFormat,PlotPath+DataID+"t"+DataGrid+"grid_climate"+PlotTop+"km")
    gsn_merge_colormaps( wks,"WhBlGrYeRe","ViBlGrWhYeOrRe")
    ;gsn_define_colormap( wks,Col(2))


    contour = new( nvar,graphic )

     ResC = True

     ResC@gsnDraw  = False
     ResC@gsnFrame = False

     ResC@tiYAxisString          = "Height Km" 

     FontHeight = 0.02

     ResC@tiXAxisFontHeightF        = FontHeight
     ResC@tiYAxisFontHeightF        = FontHeight
     ResC@tmXBLabelFontHeightF      = FontHeight
     ResC@tmYLLabelFontHeightF      = FontHeight
     ResC@gsnLeftStringFontHeightF  = FontHeight + 0.003
    ;ResC@gsnCenterStringFontHeightF= FontHeight + 0.004
    ;ResC@gsnRightStringFontHeightF = FontHeight + 0.002

     ResC@lbLabelBarOn             = True
     ResC@pmLabelBarHeightF        = 0.1
     ResC@pmLabelBarWidthF         = 0.58
     ResC@lbLabelAutoStride        = True
     ResC@pmLabelBarOrthogonalPosF = 0.01
     ResC@lbLabelFontHeightF       = FontHeight 

     ResC@tmYRMode      = "Automatic"
     ResC@trYReverse    = False
     ResC@trYMaxF       = todouble(PlotTop)
     ResC@trYMinF       = 0.

     ResC@cnInfoLabelOn                = False
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



;---------------------------------------------------------
; data reading and plotting
;---------------------------------------------------------

 do ivar=0,nvar-1
    File = addfile( DataPath + DataID + varNameF(ivar) +suffix, "r" )
    var = (/File->$varName(ivar)$(0,:,:,0)/)  ;(time,lev,cell,lon)

    if (ivar.eq.0) then
      height = (File->alt)       ; altitude, orography ist 0. in this test case
      vertc =  height/1000.      ; height in km as vertical coordinate
      vertc_t = todouble(PlotTop)
      vertc_sfc = 0.
    end if

    var!0   = "vertc"
    var&vertc = vertc
    var!1   = "lat"
    var&lat = File->lat

    vertc4plot     = fspan(vertc_t,vertc_sfc,nlev4plot)
    vertc4plot!0   = "vertc"
    vertc4plot&vertc = vertc4plot

    tmp      = linint1_Wrap( vertc, var(lat|:,vertc|:), False, vertc4plot, 0 )
    var4plot = tmp(vertc|:,lat|:)

     ResC@cnLevelSelectionMode  = "ManualLevels"
     ResC@cnMinLevelValF        =  Min(ivar)
     ResC@cnMaxLevelValF        =  Max(ivar)
     ResC@cnLevelSpacingF       =  Int(ivar)

     ResC@gsnLeftString  = LeftString(ivar)
     ;ResC@gsnRightString = ConfigStr

     ResC@gsnSpreadColorStart = colorStart(ivar)
     ResC@gsnSpreadColorEnd   = colorEnd(ivar)

     contour(ivar) = gsn_csm_contour( wks,var4plot,ResC )
  end do

     ResP                            = True
     ResP@gsnMaximize                = True
   ; ResP@gsnPaperMargin             = 0.2
     ResP@gsnPanelXWhiteSpacePercent = 5
     ResP@gsnPanelYWhiteSpacePercent = 5
     ResP@txString = "Held-Suarez experiment"
     ResP@txFontHeightF  = 0.01

     gsn_panel (wks, contour, (/4,3/),ResP)

 end
