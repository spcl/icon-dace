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
; Hui Wan ( MPI, Aug 2006 )
; Hui Wan ( MPI, Nov 2008 )
; Pilar Ripodas  ( DWD, Oct 2010)
; Martin Koehler ( DWD, Nov 2010)
; Hui Wan (MPI-M, Apr 2011)
;---------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./hgt_to_pres.ncl"
;---------------------------------------------------------
begin

;---------------------------------------------------------
; input, variable setup
;---------------------------------------------------------

 YAxis      = getenv("YAxis")

 PlotFormat = getenv("plot_file_format") 
 PlotFile   = getenv("PlotFile")
 wkOrientation = getenv("wkOrientation")

 DataFile   = getenv("DataFile")
 DataFile2  = getenv("DataFile2")

 TopHeight  = stringtofloat(getenv("TopHeight"))
 varName    = getenv("varname")
 LongName   = getenv("LongName")
 Scale      = stringtofloat(getenv("Scale"))

 Min = (/stringtofloat(getenv("Min")),   \
        -stringtofloat(getenv("DiffMax"))/)
 Max = (/stringtofloat(getenv("Max")),  \
         stringtofloat(getenv("DiffMax"))/)
 Int = (/stringtofloat(getenv("Int")),  \
         stringtofloat(getenv("DiffInt"))/)

; For the difference plot

 ColorMap   = (/getenv("ColorMap"),getenv("DiffColorMap")/)
 ColorStart = stringtoint((/getenv("ColorStart"),getenv("DiffColorStart")/))
 ColorEnd   = stringtoint((/getenv("ColorEnd"),  getenv("DiffColorEnd")  /))

 LeftString  = "Days "+(/getenv("timerange"),getenv("ref_timerange")/)
 RightString = (/getenv("resolution"),getenv("ref_resolution")/)
 PanelTitle  = (/getenv("exp_config"),"vs. "+getenv("ref_config")/)
 MainTitle   = LongName

 nlev4plot  = 101

;---------------------------------------------------------
; plot setup
;---------------------------------------------------------

 contour = new( 2, graphic )

 ; Set the workstation orientation to "landscape" if the user wants so 
 if ( .not.ismissing(wkOrientation) ) then
 
    ; check user's specification of workstation orientation
    if (wkOrientation.eq."landscape") then
       PlotFormat@wkOrientation  = "landscape"   ;NCL default is "portrait"
       PlotFormat@wkDeviceLowerX = -60           ;36
       PlotFormat@wkDeviceLowerY = 40            ;126
       PlotFormat@wkDeviceUpperX = 710           ;-36
       PlotFormat@wkDeviceUpperY = 770           ;126
    end if
 end if

 wks = gsn_open_wks( PlotFormat, PlotFile )

 ResC = True

 ResC@gsnDraw  = False
 ResC@gsnFrame = False

 if (YAxis.eq."height") then
   ResC@tiYAxisString = "Height (km)" 
 else
   ResC@tiYAxisString = "Pressure (hPa)" 
 end if     

 FontHeight = 0.02

 ResC@tiMainFontHeightF         = FontHeight + 0.003
 ResC@tiXAxisFontHeightF        = FontHeight
 ResC@tiYAxisFontHeightF        = FontHeight
 ResC@tmXBLabelFontHeightF      = FontHeight
 ResC@tmYLLabelFontHeightF      = FontHeight
 ResC@tmYRLabelFontHeightF      = FontHeight
 ResC@gsnLeftStringFontHeightF  = FontHeight + 0.003
 ResC@gsnRightStringFontHeightF = FontHeight + 0.003

 ResC@lbLabelBarOn             = True
 ResC@pmLabelBarHeightF        = 0.1
 ResC@pmLabelBarWidthF         = 0.58
 ResC@lbLabelAutoStride        = True
 ResC@pmLabelBarOrthogonalPosF = 0.01
 ResC@lbLabelFontHeightF       = FontHeight 

 if (YAxis.eq."pressure") then     
   ResC@trYMinF          = hgt_to_pres(TopHeight) 
   ResC@tmYRTickSpacingF = 2
 else
   ResC@trYReverse    = False
   ResC@trYMaxF       = TopHeight
   ResC@trYMinF       = 0.
 end if 

 ResC@cnInfoLabelOn                = False
 ResC@cnLinesOn                    = True
 ResC@gsnContourZeroLineThicknessF = 1.5
 ResC@gsnContourNegLineDashPattern = 2
 ResC@cnFillOn                     = True
 ResC@gsnSpreadColors              = True

;---------------------------------------------------------
; data reading and plotting
;---------------------------------------------------------

   File = addfile( DataFile,  "r" )
   File2= addfile( DataFile2, "r" )
 
   var  = (/File->$varName$(0,:,:,0)/)          ; (time,lev,lat,lon)
   var2 = (/File2->$varName$(0,:,:,0)/)         ; (time,lev,lat,lon)
   var  = var * Scale
   var2 = var2* Scale

  
   do istat=0,1             ; 0: single experiment, 1:diff exp-control
     if (istat .eq. 1) then
       var = var - var2
     end if
     print( varName + " min=" + min(var) + " max=" + max(var) )

     if (YAxis.eq."pressure") then

        pres_lev = File->lev
        if (pres_lev@units.eq."Pa") then
            pres_lev@units = "hPa"
            pres_lev@long_name = "Pressure (hPa)"
            pres_lev = pres_lev*0.01  ; convert to hPa
        end if
        var4plot     = (/var/)
        var4plot!0   = "p"
        var4plot&p   = pres_lev
        var4plot!1   = "lat"
        var4plot&lat = File->lat

     else
        if (varName .eq. "W") then
          height  = (File->hyai)       ; half levels (assuming orography is 0)
        else
          height  = (File->hyam)       ; full levels (assuming orography is 0)
        end if
        vertc     =  height/1000.      ; height in km as vertical coordinate
        vertc_t   = TopHeight
        vertc_sfc = 0.
 
        var!0     = "vertc"
        var&vertc = vertc
        var!1     = "lat"
        var&lat   = File->lat
 

        vertc4plot       = fspan(vertc_t,vertc_sfc,nlev4plot)
        vertc4plot!0     = "vertc"
        vertc4plot&vertc = vertc4plot

        tmp      = linint1_Wrap(vertc,var(lat|:,vertc|:),False,vertc4plot,0)
        var4plot = tmp(vertc|:,lat|:)

        delete(vertc)
        delete(tmp)
        delete(height)
     end if
 
     ResC@cnLevelSelectionMode  = "ManualLevels"
     ResC@cnMinLevelValF        = Min(istat)
     ResC@cnMaxLevelValF        = Max(istat)
     ResC@cnLevelSpacingF       = Int(istat)

     ResC@gsnSpreadColorStart   = ColorStart(istat)
     ResC@gsnSpreadColorEnd     = ColorEnd(istat)

     gsn_define_colormap(wks,ColorMap(istat)) 

     ResC@gsnLeftString   = LeftString(istat) 
     ResC@gsnRightString  = RightString(istat)
     ResC@gsnCenterString = ""
     ResC@tiMainString    = PanelTitle(istat) 
     ResC@tiXAxisString   = "Min-Mean-Max: " \
                          + min(var)+"   "+avg(var)+"   "+max(var)

     if (YAxis.eq."pressure") then
       contour(istat) = gsn_csm_pres_hgt( wks,var4plot,ResC ) 
     else
       contour(istat) = gsn_csm_contour(  wks,var4plot,ResC ) 
     end if

     ResP = True
     ResP@gsnFrame      = False
     ResP@gsnPanelLeft  = istat*0.52
     ResP@gsnPanelRight = 0.48 + ResP@gsnPanelLeft
     ResP@txString      = LongName
     ResP@txFontHeightF = FontHeight + 0.003 

     gsn_panel (wks, contour(istat), (/1,1/), ResP)

   end do
   frame(wks)

end
