; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;***********************************************************
; ICON - MJO Analysis Suite
; --------------------------
; Skripts are based on NCL MJO Climate Variability Tools
; https://www.ncl.ucar.edu/Applications/mjoclivar.shtml
; 
; These were developed following the recommendations 
; of US-CLIVAR MJO Working Group
; 
; Adapted for use at DWD by Julia Keller, FE 13, August 2015
; Works with NCL 6.3.0
; 
; Program:
;   mjo_analysis_clivar_eofs.ncl  
; Purpose:
;   Compute and plot EOFs for particular variable
; Origin:
;   This script is based on mjoclivar_12.ncl
;-----------------------------------------------------------
; Required Input:
;   Daily Anomalies of variable to be processed
;   Created by mjo_analysis_clivar_anomalies.ncl
;   
;   
; Program flow: 
;   - Read daily anomalies
;   - Apply bandpass filter
;   - Compute EOFs
;   - Plot EOFs
;
;***********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "diagnostics_cam.ncl"

begin

;###########################################################
;
; Preparation  
;
;###########################################################

;------------------------------------------------
;Requires to be initialized with the following variables
;
;------ File and Data Info------
;      var      - variable
;      infile   - name of input file
;      plotdir  - name of plot directory
;      datainfo - info on dataset
;      ymdStrt  - start date of time series
;      ymdLast  - end date of time series
;------ Definition for Bandpass ----
;      bpmin    - bandpass filter minimum
;      bpmax    - bandpass filter maximum
;      bpwgt    - bandpass filter weight
;      spd      - interval of input data in days
;      latmin     - general latitude band to be considered
;      latmax   
;      neof     - number of EOFs to be considered


 
; Define plot settings  
   wksType = "png"
   if (wksType.eq."png") then
     wksTypeLocal = "eps"
   else
     wksTypeLocal = wksType
   end if
   

   yrStrt  = ymdStrt/10000
   yrLast  = ymdLast/10000

   monsea  = (/ (/ 5, 6, 7, 8, 9,10/) \      ; May-Oct [summer]
              , (/ 1, 2, 3, 4,11,12/) /)     ; Nov-Apr [winter]


;------------------------------------------------
; Define band pass
;------------------------------------------------
   ihp      = 2                             ; bpf=>band pass filter
   nWgt     = bpwgt
   sigma    = 1.0                           ; Lanczos sigma
   fca      = 1./bpmax
   fcb      = 1./bpmin
   wgt      = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
   bpf      = (/bpmin,bpmax,bpwgt/)  


;###########################################################
;
; Get data  
;
;###########################################################


;------------------------------------------------
; Read daily anomalies
;------------------------------------------------
   varname     = var+"_anom"                  ; replace if smoothed var is wanted
   f       = addfile (infile, "r")                         
   TIME    = f->time                          ; days since ...
   YMD     = cd_calendar(TIME, -2)            ; entire (time,6)

   iStrt   = ind(YMD.eq.ymdStrt)              ; index start
   iLast   = ind(YMD.eq.ymdLast)              ; index last 
   delete(TIME)
   delete(YMD )

   if (getfilevartypes(f,varname) .eq. "short") then
      x = short2flt( f->$varname$(iStrt:iLast,{latmin:latmax},:) ) 
   else
      x = f->$varname$(iStrt:iLast,{latmin:latmax},:) 
   end if
   dimx    = dimsizes( x )
   ntim    = dimx(0)
   nlat    = dimx(1)
   mlon    = dimx(2)

   printVarSummary(x)
   printMinMax(x, True)

   time = f->time(iStrt:iLast)                ; days since ...
   ymdhms         = cd_calendar(time, 0 )     ; (ntim,2)
   mm             = floattoint( ymdhms(:,1) ) ; just months
   date  = cd_calendar( time , -2 )   ; yyyymmdd

   iSummer = ind( mm.ge.5 .and. mm.le.10)
   iWinter = ind((mm.ge.1 .and. mm.le. 4) .or.  \
                 (mm.ge.11.and. mm.le.12) )

   lat = f->lat
   lon = f->lon
   delete (varname)

;------------------------------------------------
; Apply the band pass filter to the original anomalies (lat,lon,time)
;------------------------------------------------

   xBPF = wgt_runave_Wrap (x(lat|:, lon|:, time|:), wgt, 0)


;------------------------------------------------
; Compute EOFs; no need to areal weight (15S to 15N)
;------------------------------------------------
   eof    = eofunc_Wrap(xBPF, neof, False)
   eof_ts = eofunc_ts_Wrap(xBPF,eof,False)

   printVarSummary( eof )
   printVarSummary( eof_ts )
;------------------------------------------------------------
; PLOTS
;------------------------------------------------------------
   print("-----plot EOFs-----")
   wksNameBase = plotdir+datainfo+"_"+var+"_EOF"
   wks = gsn_open_wks(wksTypeLocal,wksNameBase)
   gsn_define_colormap(wks,"BlwhRe") 
   plot = new(neof,graphic)                ; create graphic array
                                           ; only needed if paneling
; EOF patterns

   res                      = True         
   res@gsnDraw              = False        ; don't draw yet
   res@gsnFrame             = False        ; don't advance frame yet
   res@gsnSpreadColors      = True         ; spread out color table
   res@gsnStringFontHeightF = 0.015        ; make larger than default

   res@mpFillOn             = False        ; turn off map fill
   res@mpMinLatF            = latmin       ; zoom in on map
   res@mpMaxLatF            = latmax
   res@mpCenterLonF         = 180.
   res@cnFillOn             = True         ; turn on color fill
   res@cnLinesOn            = False        ; True is default
   res@cnLineLabelsOn       = False        ; True is default
   res@lbLabelBarOn         = False        ; turn off individual lb's
                                           ; set symmetric plot min/max
   symMinMaxPlt(eof,  8, False, res)       ; contributed.ncl

;   panel plot only resources
   resP                     = True         ; modify the panel plot
   resP@gsnMaximize         = True         ; large format
   resP@gsnPanelLabelBar    = True         ; add common colorbar
   resP@lbLabelAutoStride   = True         ; auto stride on labels
   resP@lbLabelFontHeightF  = 0.01
   txString            = var+": "+yrStrt+"-"+yrLast
   if (isvar("plev")) then
     txString = txString+": "+plev+" hPa"
   end if

;*******************************************
; first plot
;*******************************************
   do n=0,neof-1
     res@gsnLeftString   = "EOF "+(n+1)
     if (n.eq.0) then
       res@gsnCenterString = txString
     else
       res@gsnCenterString = ""
     end if
     res@gsnRightString  = sprintf("%5.1f", eof@pcvar(n)) +"%"
     plot(n)=gsn_csm_contour_map_ce(wks,eof(n,:,:),res)
  end do
  gsn_panel(wks,plot,(/neof,1/),resP)     ; now draw as one plot
;*******************************************
; second plot
;*******************************************
  rts           = True
  rts@gsnDraw   = False       ; don't draw yet
  rts@gsnFrame  = False       ; don't advance frame yet
  rts@gsnScale  = True        ; force text scaling               

; these four resources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

  rts@vpHeightF = 0.40        ; Changes the aspect ratio
  rts@vpWidthF  = 0.85
  rts@vpXF      = 0.10        ; change start locations
  rts@vpYF      = 0.75        ; the plot

  rts@gsnYRefLine           = 0.              ; reference line   
 ;rts@gsnXYBarChart         = True            ; create bar chart 
  rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
  rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue

; panel plot only resources
  rtsP                     = True             ; modify the panel plot
  rtsP@gsnMaximize         = True             ; large format
  rtsP@txString            = txString

; create individual plots

  yyyymmdd= cd_calendar(time, -2)
  yrfrac  = yyyymmdd_to_yyyyfrac(yyyymmdd, 0.0)
  delete(yrfrac@long_name)

  nGrd    = nlat*mlon 
  eof_ts  = eof_ts/nGrd
  rts@tiYAxisString = x@units  

  do n=0,neof-1
     rts@gsnLeftString  = "EOF "+(n+1)
     rts@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
     work = wgt_runave_Wrap(eof_ts(n,:),121, 1)
     plot(n) = gsn_csm_xy (wks,yrfrac,   work    ,rts)
     ;plot(n) = gsn_csm_xy (wks,yrfrac,eof_ts(n,:),rts)
  end do
  gsn_panel(wks,plot,(/neof,1/),rtsP)     ; now draw as one plot

  if (wksType.eq."png") then
    delete(wks)
    system("convert "+wksNameBase+".eps "+wksNameBase+".png")
    system("/bin/rm -f "+wksNameBase+".eps")
  end if
end



