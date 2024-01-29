; ICON
;
; ---------------------------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ---------------------------------------------------------------

;*************************************************
; NAO.ncl
;
; cdo -f nc -r -selyear,1979/2008 -selvar,tas,psl -setpartabc,/mnt/lustre01/pf/zmaw/m214091/ClimVarDiaGit/echam6cmip5 
;          /pool/data/MOP/AMIP/ECHAM6.3.02/T63L47/DATA/BOT_mbe0711_1976-2008  /work/mh0081/m214091/GitTest/NAO/BOT_mbe0711_1979-2008
;
;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
undef("plot_cont")
 function plot_cont  ( wks,\
                        var,x,y,cArea,res_in)
;************************************************
;  contour plots for ICON
;************************************************
local plot, res

begin

;***resource list for color contour levels
    res=res_in
    res@gsnDraw      = False
    res@gsnFrame     = False 
    res@cnLineLabelsOn = False
    res@cnFillOn        = True
    res@cnLinesOn       = False
    res@lbLabelBarOn    = False
    res@gsnPolar = "NH"     
    res@mpMinLatF    = 20.

    res@mpGridLatSpacingF =  45.            ; change latitude  line spacing
    res@mpGridLonSpacingF =  45.            ; change longitude  line spacing
    res@mpCenterLonF         = 0         ; center of map
    res@gsnSpreadColors      = True     ; use full range of colors


    res@mpGeophysicalLineColor = "gray42"
  
    res@mpGeophysicalLineThicknessF = 2.     
    res@mpGridAndLimbOn = True    ; Grid line
    res@mpFillOn = False
    res@mpOutlineOn = True  
    res@mpDataSetName = "Earth..4"
    res@sfXArray = x      ; These are 1D arrays, so a triangular mesh
    res@sfYArray = y      ; will be used to create the contours.

    contour_object = gsn_csm_contour_map(wks,var,res)


  return(contour_object)

end
;***end function
undef("overlay_plot_cont")
 function overlay_plot_cont  ( wks,\
                        var,x,y,cArea,res_in)
;************************************************
;  contour plots for ICON
;************************************************
local plot, res

begin

;***resource list for color contour levels
    res=res_in
    res@gsnDraw = False
    res@gsnFrame = False
    res@cnFillOn       = False
    res@cnLinesOn      = True
 

    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@gsnCenterString = ""

    res@gsnSpreadColors      = True     ; use full range of colors

    res@sfXArray = x      ; These are 1D arrays, so a triangular mesh
    res@sfYArray = y      ; will be used to create the contours.

    contour_object = gsn_csm_contour(wks,var,res)


  return(contour_object)

end
;***end function

;*********************************************************
;           MAIN
;*********************************************************
 begin

;***adjustment
; ==============================================================
; User defined parameters that specify region of globe and
; ==============================================================
; Traute: NAO boundaries nach Hurrell
; Traute: https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-station-based
;  latS   =  20.,  latN   =  80.,  lonL   = -90.,  lonR   =  40. define in prepare_NAO_withEOF_ECHAM


  values    = asciiread("var.txt",7,"string")
  run = values(0)
  atm_RES = values(1)
  diri = values(2)
  latS = stringtofloat(values(3))
  latN = stringtofloat(values(4))
  lonL = stringtofloat(values(5))
  lonR = stringtofloat(values(6))

;************************************************
;  Read grid information
;  define the x-, y-values and the polygon points
;************************************************

  GridDir="/pool/data/ICON/post/"+atm_RES+"_amip"
  GridInfoFileName=GridDir+"/"+atm_RES+"_amip.nc"
 
  print(GridInfoFileName)

  GridInfoFile= addfile( GridInfoFileName, "r" )

  rad2deg = 45./atan(1.)            ; radians to degrees

  x = GridInfoFile->clon *rad2deg   ; cell center, lon
  y = GridInfoFile->clat *rad2deg   ; cell center, lat
  cArea = GridInfoFile->cell_area   ; area of grid cell

; note: clon and clat are longitude and latitude of triangle centers.
;       Locations of the cell corners are given by 
;       clon_vertices and clat_vertices in the nc file.


  x!0     = "lon"
  x@units = "degrees_east"
  y!0     = "lat"
  y@units = "degrees_north"

  vlon = GridInfoFile->clon_vertices * rad2deg    ; cell longitude vertices
  vlon = where(vlon.lt.0, vlon + 360, vlon)       ; longitude: 0-360
  vlat = GridInfoFile->clat_vertices * rad2deg    ; cell lattitude vertices

;***define plot

   pltType = "pdf"                        ;-- set to pdf 
   pltName  = diri+"/NAO_"+run
   wks  = gsn_open_wks(pltType,pltName)       ; specifies a plot for NAO pattern
   gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map   

   pltName_eof  = diri+"/NAO_eof_"+run
   wkseof  = gsn_open_wks(pltType,pltName_eof)       ; specifies a plot for eof1 pattern
   gsn_define_colormap(wkseof,"amwg_blueyellowred") ; Define color map

   pltName_eof_ts  = diri+"/NAO_eof_ts_"+run   
   wksts  = gsn_open_wks(pltType,pltName_eof_ts)       ; specifies a plot eof1 timesieres

   plotxy  = new(4,graphic)
   ploteof = new(4,graphic)
   plot_tas = new(4,graphic)
   plot_psl = new(4,graphic)
   noPlot=-1

  seas  = (/"DJF","MAM","JJA","SON"/)  
 
  do i=0,dimsizes(seas)-1

   noPlot=noPlot+1

;************************************************
; Read from netCDF File: psl
;************************************************
   Unit = "[hPa]"
   fili   = "/psl_"+run+"_"+seas(i)+".nc"
   f      = addfile(diri+fili,"r")   
   psl      = f->psl(:,:)                        

   x_psl_new=psl(ncells|:,time|:)
   printVarSummary(psl)
   printMinMax(psl, True)
                            
   printVarSummary(x_psl_new)
   fili   = "/eof1_psl_"+run+"_"+seas(i)+"_fullgrid.nc"
   f      = addfile(diri+fili,"r")   
   eof_new      = f->psl(0,:) 
   printVarSummary (eof_new )

   fili   = "/eof1_ts_psl_"+run+"_"+seas(i)+".nc"
   f      = addfile(diri+fili,"r")   
   eof_ts      = f->psl(:,0,0)   
  printVarSummary (eof_ts )
        
;exit
  yyyy = (cd_calendar(eof_ts&time,-1) /100)    

  nyrs  = dimsizes(yyyy)
  print ("nyrs: "+nyrs)  ;no of years

  values    = asciiread("var_"+run+"_"+seas(i)+".txt",1,"string")
  pc_eof1 = stringtofloat(values(0))  

   

;************************************************
; Calculate the regression coefficients psl   
;************************************************

   rc_psl           = regCoef(eof_ts(time|:),x_psl_new(ncells|:,time|:))            
   
   rc_psl!0   = "ncells"    ; name dimensions
 

   printVarSummary(rc_psl)
   printMinMax (rc_psl, True)
;exit
;*************************************************
; plot eof-timeseries
;*************************************************

   rts           = True 
   rts@gsnDraw      = False
   rts@gsnFrame     = False
   rts@gsnYRefLine           = 0.              ; reference line   
   rts@gsnXYBarChart         = True            ; create bar chart 
   rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
   rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
   
   rts@tmYROn = False            ; Turn off right tickmarks.
   rts@tmXTOn = False            ; Turn off top tickmarks.   
   rts@trYMinF = -3			; set minimum Y-axis value
   rts@trYMaxF = 3			; set maximum Y-axis value
   rts@tiYAxisString = " "

   rts@gsnRightString = "  EOF1: "+sprintf("%5.1f",  pc_eof1 ) +"%"
   rts@gsnLeftString =  seas(i)+"("+yyyy(0)+"-"+yyyy(nyrs-1) +")"

   plotxy(i) = gsn_csm_xy (wksts,yyyy,eof_ts(time|:),rts)
;exit
;*************************************************
; plot eof-pattern
;*************************************************
 
   reseof = True
   reseof@cnLevelSelectionMode = "ExplicitLevels"
   reseof@cnLevels = (/-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7./)
   reseof@gsnRightString = "EOF1: "+sprintf("%5.1f",  pc_eof1 ) +"%"
   reseof@gsnLeftString = seas(i)+" ("+yyyy(0)+"-"+yyyy(nyrs-1)+")"


    eof_new=eof_new/100
    ploteof(i) =     plot_cont ( wkseof,eof_new,x,y,cArea,reseof)

;exit
;************************************************
; Read from netCDF File: tas
;************************************************
  
   fili   = "/tas_"+run+"_"+seas(i)+".nc"
   f      = addfile(diri+fili,"r") 

   x_tas      = f->tas(time|:,ncells|:)
   x_tas_new  = x_tas(ncells|:,time|:)
   
   printVarSummary(x_tas_new)                            ; [time| 720]x[lat| 91]x[lon| 180]

;exit



;************************************************
; Calculate the regression coefficients tas   
;************************************************
   rc_tas           = regCoef(eof_ts(time|:),x_tas_new(ncells|:,time|:))             
     
   rc_tas!0   = "ncells"    ; name dimensions

   printVarSummary(rc_tas)
;exit
  
;************************************************
; plotting parameters 
;************************************************

 
   res                       = True     
   res@gsnMaximize           = False           ; make large

   res@cnLevelSelectionMode = "ExplicitLevels"
   res@cnLevels = fspan(-1.05,1.05,15)

   res@gsnCenterString       = seas(i)+" ("+yyyy(0)+"-"+yyyy(nyrs-1)+")"
   

     contourMax=7
     contourMin=-7
     contourInc=1
     InfoContour="Contour from "+ contourMin+" to "+contourMax+" by "+contourInc+" hPa"
    
     res3 = True          ; PSL resources
     res3@cnLineColor = "black"
     res3@cnLineLabelsOn = False
 
     res3@gsnContourZeroLineThicknessF = 2.       ; thickness of zero contour
     res3@tiMainOn = False      
     res3@gsnContourNegLineDashPattern  = 1   
     res3@cnLineDashSegLenF = 0.08
     if (noPlot.eq.3)then
        res3@cnInfoLabelOn = True
        res3@cnInfoLabelString = InfoContour 
     else
        res3@cnInfoLabelOn = False
     end if
     res3@cnLevelSelectionMode = "ExplicitLevels"
     res3@cnLevels = ispan(contourMin,contourMax,contourInc)

     plot_tas(i) =     plot_cont ( wks,rc_tas,x,y,cArea,res)

     rc_psl=rc_psl/100
     plot_psl(i) =     overlay_plot_cont ( wks,rc_psl,x,y,cArea,res3) 

     overlay(plot_tas(i),plot_psl(i))


delete([/psl,eof_ts,x_tas,x_tas_new,x_psl_new,rc_tas,rc_psl,yyyy/])
   end do


  pres          =  True
  pres@gsnPanelLabelBar    = True                ; add common colorbar
  pres@txString           = run+" - NAO - regression: tas(color), psl "
  gsn_panel(wks, (/plot_tas(0),plot_tas(1),\
                 plot_tas(2),plot_tas(3) \
                /), (/2,2/), pres)

;*** add test at the bottom of the frame
    txres = True
    txres@txFontHeightF = 0.012
    txres@txFontColor   = "black" 
    gsn_text_ndc(wkseof,Unit,0.75,0.06,txres)
  preseof          =  True
  preseof@gsnPanelLabelBar    = True                ; add common colorbar
  preseof@txString           = run+" - NAO - EOF1 (20N-80N, 90W-40E): psl " 
  gsn_panel(wkseof, (/ploteof(0),ploteof(1),\
                    ploteof(2),ploteof(3) \
                   /), (/2,2/), preseof)

  prests          =  True 
  prests@txString           = run+" NAO index "
  gsn_panel(wksts, (/plotxy(0),plotxy(1),\
                   plotxy(2),plotxy(3) \
                   /), (/2,2/), prests)


  end
