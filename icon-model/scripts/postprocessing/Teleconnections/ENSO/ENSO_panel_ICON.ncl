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
; ENSO with ICON-grid
;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("plot_cont")
 function plot_cont  ( wks,var,x,y,cArea,res_in )
;************************************************
;  contour plots for ICON
;************************************************
local plot, res

begin

;***resource list for color contour levels
    res=res_in

    res@gsnFrame             = False             ; don't advance frame
    res@gsnDraw              = False             ; don't draw plot yet

    res@mpProjection = "WinkelTripel"
    res@mpGeophysicalLineColor = "gray42"
     
    res@mpPerimOn    = False
    res@mpGridLatSpacingF =  45.            ; change latitude  line spacing
    res@mpGridLonSpacingF =  90.           ; change longitude line spacing
;    res@mpGridLineColor   = "grey"  ; trick ncl into drawing perimeter
    res@mpGridAndLimbOn   = True           ; turn on lat/lon lines  
    res@mpFillOn = False
    res@mpCenterLonF = 210.
    res@mpOutlineOn = True  
    res@mpDataSetName = "Earth..4"


    res@cnFillOn              = True             ; turn on color
    res@cnLinesOn             = False            ; turn off contour lines
    res@cnLineLabelsOn        = False            ; turn off contour line labels
    res@lbLabelBarOn         = False   
    res@gsnLeftString        = " " 
    res@gsnRightString       = " " 

    res@gsnSpreadColors      = True     ; use full range of colors

    res@pmLabelBarOrthogonalPosF = 0.2           ; move farther to plot

    res@sfXArray = x      ; These are 1D arrays, so a triangular mesh
    res@sfYArray = y      ; will be used to create the contours.

    contour_object = gsn_csm_contour_map(wks,var,res)

    return(contour_object)

end
;***end function 

;*********************************************************
;           MAIN
;*********************************************************
begin
 


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
   pltName  = diri+"/ENSO_"+run
   print("pltName:  "+pltName)
   wks  = gsn_open_wks(pltType,pltName)       ; specifies a plot
   gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map

   listplots = new(4,graphic)
   listplots_tas = new(4,graphic)
   listplots_psl = new(4,graphic)
   noPlot=-1
 
  seas  = (/"DJF","MAM","JJA","SON"/)
  
   do i=0,dimsizes(seas)-1
print("seas "+seas(i))

    noPlot=noPlot+1
;************************************************
; Read from netCDF File: psl
;************************************************
   
   fili   = "psl_"+run+"_"+seas(i)+"_"+atm_RES+"_detrend.nc"
   f      = addfile(diri+fili,"r") 

   x_psl      = f->psl(time|:,ncells|:)   ; dims: (time,ncells)                       
   x_psl_new  = x_psl(ncells|:,time|:)
   
   printVarSummary(x_psl_new)                            

;exit
;************************************************
; Read from netCDF File: tas
;************************************************
  
   fili   = "tas_"+run+"_"+seas(i)+"_"+atm_RES+"_detrend.nc"
   f      = addfile(diri+fili,"r") 

   x_tas      = f->tas(time|:,ncells|:)   ; dims: (time,ncells)        
   x_tas_new  = x_tas(ncells|:,time|:)
   
   printVarSummary(x_tas_new)                            

   yyyymm = cd_calendar(x_tas&time, -1)
   yyyy   = yyyymm/100

   nyrs  = dimsizes(yyyy)
   print ("nyrs: "+nyrs)  ;no of years

   print("Time: "+yyyy(0)+"-"+yyyy(nyrs-1) )
;exit

;****************************
; Read from netCDF file: ts nino3.4
;***********************************************
  
   fili   = "ts_"+run+"_"+seas(i)+"_"+atm_RES+"_nino34_detrend.nc"
   f      = addfile(diri+fili,"r") 
   printVarSummary(f) 
    print(f)          
   ts_nino34      = f->ts(time|:,lat|0,lon|0)

;   ts_n      = f->ts(time|:,lat|0,lon|0)

   printVarSummary(ts_nino34)                           
;  printVarSummary(ts_n) 


;************************************************
; Calculate the regression coefficients tas   
;************************************************

  rc_tas           = regCoef(ts_nino34,x_tas_new(ncells|:,time|:))            
   
;   printVarSummary(rc_tas)  
     
  rc_tas!0   = "ncells"    ; name dimensions

   printVarSummary(rc_tas)
   printMinMax (rc_tas, True)

  filo = "rc_tas_"+run+"_"+seas(i)+".nc"
  system("/bin/rm -f "+diri+filo)   ; remove any pre-existing file
  ncdf = addfile(diri+filo,"c")     ; open output netCDF file
                                    ; make time an UNLIMITED dimension


  ncdf->rc_tas = rc_tas

;;cdo -setgrid,/pool/data/ICON/post/r2b4_amip/r2b4_amip.nc rc_tasmbe0780_SON.nc rc_tasmbe0780_SONgrid.nc

;exit
  


;************************************************
; Calculate the regression coefficients psl   
;************************************************
   rc_psl           = regCoef(ts_nino34,x_psl_new(ncells|:,time|:))            
   
   rc_psl!0   = "ncells"    ; name dimensions

   printVarSummary(rc_psl)
   printMinMax (rc_psl, True)
;exit

;************************************************
; plotting parameters 
;************************************************

   
   res                       = True     
;   res@gsnMaximize           = True             ; make large

   res@cnLevelSelectionMode = "ExplicitLevels"
   res@cnLevels = (/-4,-3,-2,-1.5,-1,-.5,-.25,0,.25,.5,1,1.5,2,3.,4/)

   res@gsnCenterString       = seas(i)+" ("+yyyy(0)+"-"+yyyy(nyrs-1) +")"

     contourMax=8
     contourMin=-8
     contourInc=1

     res3 = True                        ; PSL resources
     res3@cnLineColor = "black"
     res3@cnLineLabelsOn = False
     
     if (noPlot.eq.3)then
        res3@cnInfoLabelOn = True
     else
        res3@cnInfoLabelOn = False
     end if
     res3@cnInfoLabelString =  "Contour from "+ contourMin+" to "+contourMax+" by "+contourInc+" hPa"
     res3@gsnContourZeroLineThicknessF = 2.       ; thickness of zero contour
     res3@tiMainOn = False      
     res3@gsnContourNegLineDashPattern  = 1   
     res3@cnLineDashSegLenF = 0.08
     res3@gsnDraw = False
     res3@gsnFrame = False
     res3@gsnLeftString = ""
     res3@gsnRightString = ""
     res3@gsnCenterString = ""
     res3@cnLevelSelectionMode = "ExplicitLevels"
     res3@cnLevels = ispan(contourMin,contourMax,contourInc)

     res3@sfXArray = x      ; These are 1D arrays, so a triangular mesh
     res3@sfYArray = y      ; will be used to create the contours.

     listplots_tas(noPlot) =     plot_cont ( wks,rc_tas,x,y,cArea,res)

     rc_psl=rc_psl/100
     listplots_psl(noPlot) = gsn_csm_contour(wks,rc_psl,res3) 

     overlay(listplots_tas(noPlot),listplots_psl(noPlot))

     delete([/ts_nino34,x_tas,x_psl,x_tas_new,x_psl_new,rc_tas,rc_psl,yyyymm,yyyy/])
   end do
  
   
  pres          =  True
  pres@txString           = run+ " -  Nino3.4 regression: tas(color), psl    "
  pres@gsnPanelLabelBar    = True                ; add common colorbar
  pres@lbLabelFontHeightF  = 0.007               ; make labels smaller
  gsn_panel(wks, (/listplots_tas(0),listplots_tas(1),\
                    listplots_tas(2),listplots_tas(3) \
                   /), (/2,2/), pres)

  end
