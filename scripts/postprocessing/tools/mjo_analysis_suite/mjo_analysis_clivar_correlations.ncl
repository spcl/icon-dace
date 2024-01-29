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
;   mjo_analysis_clivar_correlations.ncl  
; Purpose:
;   Compute various sorts of correlation among TOT_PREC,
;   OLR and Wind daily anomalies.
; Origin:
;   This script is a merger of NCL scripts mjoclivar_9.ncl
;   mjoclivar_11.ncl, mjoclivar_14.ncl and mjoclivar_16.ncl
;-----------------------------------------------------------
; Required Input:
;   Daily Anomalies of variable to be processed
;   Created by mjo_analysis_clivar_anomalies.ncl
;   
;   
; Program flow: 
;   - Read daily anomalies for all relevant variables
;   - Check for their temporal and spatial consistency
;   - Compute and plot lagged correlation for entire dataset
;     (mjoclivar_9.ncl)
;   - Compute and plot cross-spectra (mjoclivar_11.ncl)
;   - Compute and plot multivariate EOFs (mjoclivar_14.ncl)
;   - Write out PC1 and PC2 for RMM Index
;   - Create and plot MJO Lifecycle (mjoclivar_16)
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
;      (*name  - name of variable)
;      (*file  - name of file)	 
;      precname 
;      olrname
;      u850name            			  
;      u200name 
;      v850name 
;      precfile
;      olrfile 
;      u850file
;      u200file  
;      v850file
;      outpath  - path to output folder
;      plotdir  - name of plot directory
;      datainfo  - info on dataset
;      ymdStrt  - start date of time series
;      ymdLast  - end date of time series
;      dayint   - interval of input data in days
;------ Definition for Bandpass ----
;      bpmin    - bandpass filter minimum
;      bpmax    - bandpass filter maximum
;      bpwgt    - bandpass filter weight
;------ General Area to be considered ----
;      latmin   - general latitude band to be considered
;      latmax   
;------ Area to be considered for Indian Ocean----
;      latS_IO
;      latN_IO     
;      lonW_IO    
;      lonE_IO     
;------ Area to be considered for lag-corr----
;      latS_bnd
;      latN_bnd    
;      lonW_bnd    
;      lonE_bnd 
;------ Switches for various plots (on/off)------
;      neof      - number of EOFs to be considered
;      lagcorr   - Lagged Correlation plot
;      xspectra  - Cross Spectrum Plot
;--------------------------------------------------   
; Specify seasons 
   nameSeason = (/"winter","summer","annual"/)
   nameRegion = "IO"        ; Indian Ocean base region



; Define plot settings 
   wksType = "png"
   if (wksType.eq."png") then
     wksTypeLocal = "eps"
   else
     wksTypeLocal = wksType
   end if
   

   yrStrt  = ymdStrt/10000
   yrLast  = ymdLast/10000

 
;------------------------------------------------
; Define band pass
;------------------------------------------------
   ihp      = 2                             ; bpf=>band pass filter
   nWgt     = bpwgt
   sigma    = 1.0                           ; Lanczos sigma
   fca      = dayint/bpmax
   fcb      = dayint/bpmin
   print(fca)
   print(fcb)
   wgt      = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
   bpf      = (/bpmin,bpmax,bpwgt/)  
   spd      = 1/dayint



;###########################################################
;
; Get data  
;
;###########################################################

  print("-------Get data and prepare-------")


;------------------------------------------------
; Read daily anomalies for all five variables
; check for consistency in dimensions
;------------------------------------------------
   olr_var    = olrname+"_anom"                  ; replace if smoothed var is wanted
   folr         = addfile (olrfile, "r")                         
   prec_var   = precname+"_anom"                  ; replace if smoothed var is wanted
   fprec         = addfile (precfile, "r")                         
   u850_var   = u850name+"_anom"                  ; replace if smoothed var is wanted
   fu850      = addfile (u850file, "r")                         
   v850_var   = v850name+"_anom"                  ; replace if smoothed var is wanted
   fv850      = addfile (v850file, "r")                         
   u200_var   = u200name+"_anom"                  ; replace if smoothed var is wanted
   fu200      = addfile (u200file, "r")                         
   TIME    = folr->time                          ; days since ...
   YMD     = cd_calendar(TIME, -2)            ; entire (time,6)

   iStrt   = ind(YMD.eq.ymdStrt)              ; index start
   iLast   = ind(YMD.eq.ymdLast)              ; index last 
   lat = fprec->lat
   lon = fprec->lon
  
   delete(TIME)
   delete(YMD)

   if (getfilevartypes(folr,olr_var) .eq. "short") then
      olr = short2flt( folr->$olr_var$(iStrt:iLast,{latmin:latmax},:) ) 
   else
      olr = folr->$olr_var$(iStrt:iLast,{latmin:latmax},:) 
   end if
;   printVarSummary(olr)
;   printMinMax(olr, True)
   delete(olr_var)

   if (getfilevartypes(fprec,prec_var) .eq. "short") then
      prec = short2flt( fprec->$prec_var$(iStrt:iLast,{latmin:latmax},:) ) 
   else
      prec = fprec->$prec_var$(iStrt:iLast,{latmin:latmax},:) 
   end if
;   printVarSummary(prec)
;   printMinMax(prec, True)
   delete(prec_var)

   if (getfilevartypes(fu850,u850_var) .eq. "short") then
      u850 = short2flt( fu850->$u850_var$(iStrt:iLast,{latmin:latmax},:) ) 
   else
      u850 = fu850->$u850_var$(iStrt:iLast,{latmin:latmax},:) 
   end if
;   printVarSummary(u850)
;   printMinMax(u850, True)
   delete(u850_var)



   if (getfilevartypes(fu200,u200_var) .eq. "short") then
      u200 = short2flt( fu200->$u200_var$(iStrt:iLast,{latmin:latmax},:) ) 
   else
      u200 = fu200->$u200_var$(iStrt:iLast,{latmin:latmax},:) 
   end if
;   printVarSummary(u200)
;   printMinMax(u200, True)
   delete(u200_var)

   if (getfilevartypes(fv850,v850_var) .eq. "short") then
      v850 = short2flt( fv850->$v850_var$(iStrt:iLast,{latmin:latmax},:) ) 
   else
      v850 = fv850->$v850_var$(iStrt:iLast,{latmin:latmax},:) 
   end if
;   printVarSummary(v850)
;   printMinMax(v850, True)
   delete(v850_var)

   time1 = folr->time(iStrt:iLast)                ; days since ...
   ymdhms1         = cd_calendar(time1, 0 )     ; (ntim,2)
   mm1             = floattoint( ymdhms1(:,1) ) ; just months
   date_1  = cd_calendar( time1 , -2 )   ; yyyymmdd
   time2 = fprec->time(iStrt:iLast)                ; days since ...
   ymdhms2         = cd_calendar(time2, 0 )     ; (ntim,2)
   mm2             = floattoint( ymdhms2(:,1) ) ; just months
   date_2  = cd_calendar( time2 , -2 )   ; yyyymmdd
   time3 = fu850->time(iStrt:iLast)                ; days since ...
   ymdhms3         = cd_calendar(time3, 0 )     ; (ntim,3)
   mm3             = floattoint( ymdhms3(:,1) ) ; just months
   date_3  = cd_calendar( time3 , -2 )   ; yyyymmdd
   time4 = fu200->time(iStrt:iLast)                ; days since ...
   ymdhms4         = cd_calendar(time4, 0 )     ; (ntim,3)
   mm4             = floattoint( ymdhms4(:,1) ) ; just months
   date_4  = cd_calendar( time4 , -2 )   ; yyyymmdd
   time5 = fv850->time(iStrt:iLast)                ; days since ...
   ymdhms5         = cd_calendar(time5, 0 )     ; (ntim,3)
   mm5             = floattoint( ymdhms5(:,1) ) ; just months
   date_5  = cd_calendar( time5 , -2 )   ; yyyymmdd
   imon    = floattoint( ymdhms1(:,1) )   ; convenience
   iday    = floattoint( ymdhms1(:,2) )   ; subscripts must be integer

 ;------------------------------------------------
 ; check for consistency in dimensions
 ;------------------------------------------------

   if (.not.all(date_1.eq.date_2)) then
      print("date mismatch: exit")
      exit
   end if
   if (.not.all(date_2.eq.date_3)) then
      print("date mismatch: exit")
      exit
   end if
   if (.not.all(date_2.eq.date_4)) then
      print("date mismatch: exit")
      exit
   end if
   if (.not.all(date_2.eq.date_5)) then
      print("date mismatch: exit")
      exit
   end if

   if (.not.all(dimsizes(olr).eq.dimsizes(prec))) then
      print("****size mismatch: exit****")
      print("dimsizes(olr)="+dimsizes(olr))
      print("dimsizes(prec)="+dimsizes(prec))
      exit  
   end if
   if (.not.all(dimsizes(prec).eq.dimsizes(u850))) then
      print("****size mismatch: exit****")
      print("dimsizes(prec)="+dimsizes(prec))
      print("dimsizes(u850)="+dimsizes(u850))
      exit
   end if
   if (.not.all(dimsizes(u850).eq.dimsizes(v850))) then
      print("****size mismatch: exit****")
      print("dimsizes(u850)="+dimsizes(u850))
      print("dimsizes(v850)="+dimsizes(v850))
      exit
   end if
   if (.not.all(dimsizes(v850).eq.dimsizes(u200))) then
      print("****size mismatch: exit****")
      print("dimsizes(v850)="+dimsizes(v850))
      print("dimsizes(u200)="+dimsizes(u200))
      exit
   end if

   if (.not.all(olr&lon.eq.prec&lon)) then
      print("****longitude mismatch: exit****")
      exit
   end if
   if (.not.all(prec&lon.eq.u850&lon)) then
      print("****longitude mismatch: exit****")
      exit
   end if
   if (.not.all(u850&lon.eq.v850&lon)) then
      print("****longitude mismatch: exit****")
      exit
   end if
   if (.not.all(v850&lon.eq.u200&lon)) then
      print("****longitude mismatch: exit****")
      exit
   end if
   time   = time1 
   ymdhms = ymdhms1        
   mm     = mm1             
   date   = date_1 
   delete([/date_1,date_2,date_3,date_4,date_5, \
            time1,time2,time3,time4,time5, \
            ymdhms1,ymdhms1,ymdhms1,ymdhms1,ymdhms1,mm1,mm2,mm3,mm4,mm5 /])
  
;###########################################################
;
; Create lagged correlation plots
;
;###########################################################
   if (lagcorr .eq. "true") then
   print("-------Create lagged correlation plots-------")

   ;------------------------------------------------
   ; Further specify parameters
   ;------------------------------------------------

     wyIO   = fprec->lat({latS_IO:latN_IO})  
     wyIO   = cos(0.017459*wyIO)         ; spatial wgt
     wyU    = fu850->lat({latS_IO:latN_IO})
     wyU    = cos(0.017459*wyU)          ; MJO Clivar says cos(lat) wgting  

   ;------------------------------------------------
   ; Create wgted area average of the base IO precip series  (time)
   ; Really, no need to area weight here .... area is very small.
   ;------------------------------------------------

    precIO    = wgt_areaave_Wrap(prec(:,{latS_IO:latN_IO},{lonW_IO:lonE_IO}), wyIO, 1., 0)
    precIO    = dtrend (precIO, False)                                   ; rmv overall trend
    precIO    = wgt_runave_leftdim( precIO, wgt, 0 )                     ; apply filter
  
   ;------------------------------------------------
   ; Create LAT average of the global Precip and U series  (time,lon)
   ; Really, no need to area weight here
   ;------------------------------------------------

    prec_timeLon = dim_avg_Wrap(prec( time|:,lon|:,{lat|latS_bnd:latN_bnd}) )  ; (time,lon)
    prec_timeLon = dtrend_leftdim (prec_timeLon, False)                ; rmv overall trend
    prec_timeLon = wgt_runave_leftdim(prec_timeLon, wgt, 0 )          ; apply filter

    u850_timeLon = dim_avg_Wrap(u850( time|:,lon|:,{lat|latS_bnd:latN_bnd}) )  ; (time,lon)
    u850_timeLon = dtrend_leftdim (u850_timeLon, False)                ; rmv overall trend
    u850_timeLon = wgt_runave_leftdim(u850_timeLon, wgt, 0 )

   ;------------------------------------------------
   ; Create LON average of the global Precip and U series  (time,lat)
   ; Really, no need to area weight here
   ;------------------------------------------------

    prec_timeLat   = dim_avg_Wrap(prec( time|:,lat|:,{lon|lonW_bnd:lonE_bnd}) )  ; (time,lat)
    prec_timeLat   = dtrend_leftdim(prec_timeLat, False)                ; rmv overall trend
    prec_timeLat   = wgt_runave_leftdim(prec_timeLat, wgt, 0 )          ; apply filter

    u850_timeLat   = dim_avg_Wrap(u850( time|:,lat|:,{lon|lonW_bnd:lonE_bnd}) )  ; (time,lat)
    u850_timeLat   = dtrend_leftdim(u850_timeLat, False)                ; rmv overall trend
    u850_timeLat   = wgt_runave_leftdim(u850_timeLat, wgt, 0 )
   ;------------------------------------------------
   ; Calculate/Plot the mean seasonal cross-correlations at +/- 'mxlag' lags
   ;------------------------------------------------
    wksNameBase = datainfo+"_Lag_Correlation_"
    optXcor   = False
    mxlag     = floattoint(26*spd)
    nSeason   = dimsizes(nameSeason)

    optPlot                  = True
    optPlot@gsnLeftString    = precname + "(color)"    
    optPlot@gsnRightString   = u850name + "(lines)"  
    optPlot@smth9            = 0.25                ; local spatial smoothing

    timePeriod               =  ymdStrt+"-"+ymdLast  ; panel title 

    do ns=0,nSeason-1                              ; loop over each season
      prec_timelon = mjo_xcor_lag (precIO, prec_timeLon, date, mxlag, nameSeason(ns), spd,  optXcor)
      u850_timelon = mjo_xcor_lag (precIO, u850_timeLon, date, mxlag, nameSeason(ns),spd, optXcor)

      prec_timelat = mjo_xcor_lag (precIO, prec_timeLat, date, mxlag, nameSeason(ns),spd, optXcor)
      u850_timelat = mjo_xcor_lag (precIO, u850_timeLat, date, mxlag, nameSeason(ns),spd, optXcor)
 
      ;------------------------------
	  ; adapt lag steps to actual days
      ;------------------------------
      prec_timelon&lag=ispan(floattoint(-mxlag/spd),floattoint(mxlag/spd),floattoint(1/spd))
      prec_timelat&lag=ispan(floattoint(-mxlag/spd),floattoint(mxlag/spd),floattoint(1/spd))
      u850_timelon&lag=ispan(floattoint(-mxlag/spd),floattoint(mxlag/spd),floattoint(1/spd))
      u850_timelat&lag=ispan(floattoint(-mxlag/spd),floattoint(mxlag/spd),floattoint(1/spd))

      optPlot@txString      = datainfo +" "+nameSeason(ns)+": "+timePeriod
      wksNameSeason = wksNameBase+"_"+precname+"_"+u850name+"_"+nameSeason(ns)
      mjo_xcor_lag_ovly_panel(prec_timelon, u850_timelon, prec_timelat, u850_timelat \
                            ,wksTypeLocal, plotdir, wksNameSeason, optPlot)
      print("convert "+plotdir+wksNameSeason+".eps "+plotdir+wksNameSeason+".png")
      system("convert "+plotdir+wksNameSeason+".eps "+plotdir+wksNameSeason+".png")
      system("/bin/rm -f "+plotdir+wksNameSeason+".eps")
     end do
  
   ;------------------------------------------------
   ; delete unnecessary variables
   ;------------------------------------------------
     delete([/wyIO,wyU, precIO, prec_timeLon, u850_timeLon, \
			  prec_timeLat, u850_timeLat,  \ 
              prec_timelon, u850_timelon,prec_timelat,u850_timelat/])
   end if
;###########################################################
;
; Create cross-spectra plots
;
;###########################################################

   if (xspectra .eq. "true") then
    print("-------Create Cross Spectra Plots-------")
    if (isvar("wksNameBase")) then
       delete(wksNameBase)
     end if
     if (olr&lon(0) .lt. 0) then
       olr = lonFlip( olr )       ; force 0 to 360 to match Y
     end if
     if (u850&lon(0) .lt. 0) then
       u850 = lonFlip( u850 )       ; force 0 to 360 to match Y
     end if
   ;------------------------------------------------
   ; Calculate the cross-spectra
   ;------------------------------------------------

     segLen = floattoint(256/dayint)
     segOverLap = floattoint(-50/2)
     STC    =  mjo_cross (olr(:,{-15:15},:), u850(:,{-15:15},:), segLen, segOverLap, False)

     printVarSummary(STC)
     print("===") 
     do n=0,15
       print(STC@varName(n)+":  min="+min(STC(n,:,:)) \
                         +":  max="+max(STC(n,:,:)) )
     end do
     print( STC@prob+"   "+STC@prob_coh2 )

   ;------------------------------------------------
   ; plot
   ;------------------------------------------------
     wksNameBase=datainfo+"_Cross_Spectra_"+olrname+"_"+u850name
     ; default
     opt = True
     opt@txString = datainfo + " Default: -15S-15N: "+yrStrt+"-"+yrLast
     mjo_cross_plot(STC ,plotdir, wksType ,wksNameBase+"all", opt)

     ;   sample of plotting just values a where coherence^2 >= 0.925
     opt@pltProb  = 0.925        
     np  = ind(STC@prob .eq. opt@pltProb)
     opt@txString = "prob="+opt@pltProb+" ==> "+STC@prob_coh2(np)
     mjo_cross_plot(STC ,plotdir, wksType ,wksNameBase+"prob", opt)
     delete(opt@pltProb)

     ; sample ... do not plot the phases
     opt@pltPhase = False        ; no phases ploted       
     opt@txString = "No Phases Plotted"
     mjo_cross_plot(STC ,plotdir,wksType ,wksNameBase+"nophase", opt)

   ;------------------------------------------------
   ; delete unnecessary variables
   ;------------------------------------------------
    delete(STC)
   end if

;###########################################################
;
; Create univariate EOFs
;
;###########################################################

   if ( unieofs .eq. "true") then
    print("-------Create univariate EOFs-------")
    if (isvar("wksNameBase")) then
       delete(wksNameBase)
     end if
     ueofolr      = dim_avg_Wrap( olr(time|:,lon|:,{lat|-15:15}))                      ; (time,lon)
     printVarSummary(ueofolr)
     ueofu850     = dim_avg_Wrap(u850(time|:,lon|:,{lat|-15:15}))                     ; (time,lon)
     printVarSummary(ueofu850)
     ueofu200     = dim_avg_Wrap(u200(time|:,lon|:,{lat|-15:15}))                     ; (time,lon)
     printVarSummary(ueofu200)
     printVarSummary(wgt)
     dimw    = dimsizes( olr )
     ntim    = dimw(0)
     nlat    = dimw(1)
     mlon    = dimw(2)
   ;------------------------------------------------
   ; Apply the band pass filter to the original anomalies
   ;------------------------------------------------
     ueofolrbp   = wgt_runave_Wrap ( ueofolr(lon|:, time|:), wgt, 0) ; (time,lon)
     ueofu850bp  = wgt_runave_Wrap (ueofu200(lon|:, time|:), wgt, 0)
     ueofu200bp  = wgt_runave_Wrap (ueofu850(lon|:, time|:), wgt, 0)
   ;------------------------------------------------
   ; remove temporal means of band pass series: *not* necessary 
   ;------------------------------------------------
     ueofolrbp   = dim_rmvmean(ueofolrbp)              ; (time,lon)
     ueofu850bp  = dim_rmvmean(ueofu850bp)
     ueofu200bp  = dim_rmvmean(ueofu200bp)

   ;------------------------------------------------
   ; Compute the temporal variance 
   ;------------------------------------------------
     vari_ueofolr  = dim_variance_Wrap(ueofolrbp)     ; (lon)
     vari_ueofu850 = dim_variance_Wrap(ueofu850bp)
     vari_ueofu200 = dim_variance_Wrap(ueofu200bp)

   ;------------------------------------------------
   ; Compute the zonal mean of the temporal variance
   ;------------------------------------------------
     zavg_vari_ueofolr  = dim_avg_Wrap(vari_ueofolr)     
     zavg_vari_ueofu850 = dim_avg_Wrap(vari_ueofu850)
     zavg_vari_ueofu200 = dim_avg_Wrap(vari_ueofu200)
   ;------------------------------------------------
   ; Normalize by sqrt(avg_var*)
   ;------------------------------------------------
     ueofolrbp   = ueofolrbp/sqrt(zavg_vari_ueofolr)          ; (time,lon)
     ueofu850bp  = ueofu850bp/sqrt(zavg_vari_ueofu850)
     ueofu200bp  = ueofu200bp/sqrt(zavg_vari_ueofu200)
;************************************************
; Compute EOFs
;************************************************
     eof_olr     = eofunc_Wrap(ueofolrbp    , neof, False)      ; (neof,lon)
     eof_ts_olr  = eofunc_ts_Wrap(ueofolrbp , eof_olr ,False)   ; (neof,time)
     print("==============")
     printVarSummary(eof_olr )
     printMinMax(eof_olr , True)

     eof_u850    = eofunc_Wrap(ueofu850bp   , neof, False)      ; (neof,lon)
     eof_ts_u850 = eofunc_ts_Wrap(ueofu850bp, eof_u850,False)   ; (neof,time)
     print("==============")
     printVarSummary(eof_u850)
     printMinMax(eof_u850, True)

     eof_u200    = eofunc_Wrap(ueofu200bp   , neof, False)      ; (neof,lon)
     eof_ts_u200 = eofunc_ts_Wrap(ueofu200bp,eof_u200,False)    ; (neof,time)
     print("==============")
     printVarSummary(eof_u200)
     printMinMax(eof_u200, True)

   ;************************************************
   ; Compute cross correlation of each EOF time series at zero-lag
   ;************************************************
     mxlag     = floattoint(26*spd)

     rlag_ou8  = esccr(eof_ts_olr , eof_ts_u850, mxlag)          ; (N,mxlag+1)
     rlag_ou2  = esccr(eof_ts_olr , eof_ts_u200, mxlag)          ; (N,mxlag+1)
     rlag_u8u2 = esccr(eof_ts_u850, eof_ts_u200, mxlag)          ; (N,mxlag+1)

     rlag_u8o  = esccr(eof_ts_u850, eof_ts_olr , mxlag)          ; (N,mxlag+1)
     rlag_u2o  = esccr(eof_ts_u200, eof_ts_olr , mxlag)          ; (N,mxlag+1)
     rlag_u2u8 = esccr(eof_ts_u200, eof_ts_u850, mxlag)          ; (N,mxlag+1)

     ccr_ou8   = new ( (/neof,2*mxlag+1/), float)    
     ccr_ou8(:,0:mxlag-1) = rlag_ou8(:,1:mxlag:-1)  ; "negative lag", -1 reverses order
     ccr_ou8(:,mxlag:)    = rlag_u8o(:,0:mxlag)     ; "positive lag"

     ccr_ou2   = new ( (/neof,2*mxlag+1/), float)    
     ccr_ou2(:,0:mxlag-1) = rlag_ou2(:,1:mxlag:-1)  ; "negative lag", -1 reverses order
     ccr_ou2(:,mxlag:)    = rlag_u2o(:,0:mxlag)     ; "positive lag"

     ccr_u8u2  = new ( (/neof,2*mxlag+1/), float)    
     ccr_u8u2(:,0:mxlag-1) = rlag_u2u8(:,1:mxlag:-1); "negative lag", -1 reverses order
     ccr_u8u2(:,mxlag:)    = rlag_u8u2(:,0:mxlag)   ; "positive lag"

     ccr_ou8  = where( ccr_ou8.gt. 1.0, 1.0, ccr_ou8)   
     ccr_ou8  = where( ccr_ou8.lt.-1.0,-1.0, ccr_ou8)     
     ccr_ou2  = where( ccr_ou2.gt. 1.0, 1.0, ccr_ou2)   
     ccr_ou2  = where( ccr_ou2.lt.-1.0,-1.0, ccr_ou2)     
     ccr_u8u2 = where(ccr_u8u2.gt. 1.0, 1.0, ccr_u8u2)   
     ccr_u8u2 = where(ccr_u8u2.lt.-1.0,-1.0, ccr_u8u2)     

;------------------------------------------------------------
; PLOTS
;------------------------------------------------------------
     wksNameBase=plotdir+datainfo
     yyyymmdd = cd_calendar(time, -2)
     yrfrac   = yyyymmdd_to_yyyyfrac(yyyymmdd, 0.0)

     day      = ispan(-mxlag, mxlag, 1)
  
     wks1 = gsn_open_wks(wksType,wksNameBase+"_uni_eof")
     wks2 = gsn_open_wks(wksType,wksNameBase+"_lag_uni_eof")
     plot = new(3,graphic)                ; create graphic array
                                          ; only needed if paneling
     rts           = True
     rts@gsnDraw   = False       ; don't draw yet
     rts@gsnFrame  = False       ; don't advance frame yet
     rts@gsnScale  = True        ; force text scaling               

     rts@vpHeightF = 0.40        ; Changes the aspect ratio
     rts@vpWidthF  = 0.85
     rts@vpXF      = 0.10        ; change start locations
     rts@vpYF      = 0.75        ; the plot
     rts@xyLineThicknesses = (/2, 2, 2/)
     rts@gsnYRefLine            = 0.                  ; reference line   
     rts@gsnXRefLine            = 0.                  ; reference line   

     rts@pmLegendDisplayMode    = "Always"            ; turn on legend
     rts@pmLegendSide           = "Top"               ; Change location of 
     rts@pmLegendParallelPosF   = 0.86                ; move units right
     rts@pmLegendOrthogonalPosF = -0.50               ; move units down
     rts@pmLegendWidthF         = 0.15                ; Change width and
     rts@pmLegendHeightF        = 0.15                ; height of legend.
     rts@lgLabelFontHeightF     = 0.0175

     rtsP                       = True                ; modify the panel plot
     rtsP@gsnMaximize           = True                ; large format

     z_eof = new ( (/3,mlon/), typeof(ueofolr), getFillValue(ueofolr) )
     z_ts  = new ( (/3,ntim/), typeof(ueofolr), getFillValue(ueofolr) )

;------------------------------------------------------------
; EOF spatial patterns   [neof,time]
;------------------------------------------------------------
     rtsP@txString              = "Univariate EOF: 15S-15N: "+yrStrt+"-"+yrLast 
  
     do n=0,neof-1
       z_eof(0,:)  = eof_olr (n,:)
       z_eof(1,:)  = eof_u850(n,:)
       z_eof(2,:)  = eof_u200(n,:)
       delete(z_eof@long_name)
       rts@xyExplicitLegendLabels = \ 
          (/"OLR:"+sprintf("%5.1f", eof_olr@pcvar(n))  +"%" \
           ,"U850:"+sprintf("%5.1f",eof_u850@pcvar(n)) +"%" \
           ,"U200:"+sprintf("%5.1f",eof_u200@pcvar(n)) +"%" /) 
  
       rts@gsnLeftString  = "EOF "+(n+1)
       plot(n) = gsn_csm_xy (wks1,lon,z_eof,rts)
     end do
     gsn_panel(wks1,plot(0:1),(/2,1/),rtsP)     ; now draw as one plot

    delete(rts@xyExplicitLegendLabels) 


   ;------------------------------------------------------------
   ; plot cross correlations   [neof,time]
   ;------------------------------------------------------------
     rts@xyExplicitLegendLabels = (/ "EOF 1", "EOF 2"/) 
     rts@pmLegendOrthogonalPosF = -1.135              ; move units down
     rts@pmLegendParallelPosF   =  0.095              ; move units right
     rts@pmLegendWidthF         =  0.15               ; Change width and
     rts@pmLegendHeightF        =  0.10               ; height of legend.

     delete(rts@gsnLeftString)
     rts@gsnCenterString = "r(OLR,U850)"
     plot(0) = gsn_csm_xy (wks2,day,ccr_ou8 ,rts)
     rts@gsnCenterString = "r(OLR,U200)"
     plot(1) = gsn_csm_xy (wks2,day,ccr_ou2 ,rts)
     rts@gsnCenterString = "r(U850,U200)"
     plot(2) = gsn_csm_xy (wks2,day,ccr_u8u2,rts)
     rtsP@txString = "Lag(r): Univariate EOF: 15S-15N: "+yrStrt+"-"+yrLast 
     gsn_panel(wks2,plot,(/3,1/),rtsP)     ; now draw as one plot
  
     delete([/day,rts@gsnCenterString,rtsP@txString,rts@xyExplicitLegendLabels,rts@gsnLeftString/])
   end if



;###########################################################
;
; Create multivariate EOFs
;
;###########################################################

   if ( multieofs .eq. "true") then
    print("-------Create multivariate EOFs-------")
    if (isvar("wksNameBase")) then
       delete(wksNameBase)
     end if
     eofvar1     = dim_avg_n_Wrap(olr(:,{-15:15},:), 1)                         ; (time,lon)
     printVarSummary(eofvar1)
     eofvar2     = dim_avg_n_Wrap(u850(:,{-15:15},:), 1)                         ; (time,lon)
     printVarSummary(eofvar2)
     eofvar3     = dim_avg_n_Wrap(u200(:,{-15:15},:), 1)                         ; (time,lon)
     printVarSummary(eofvar3)

     dimw    = dimsizes( olr )
     ntim    = dimw(0)
     nlat    = dimw(1)
     mlon    = dimw(2)
  
   ;------------------------------------------------
   ; Apply the band pass filter to the original anomalies
   ;------------------------------------------------
     eofvar1bp  = wgt_runave_n_Wrap (eofvar1, wgt, 0, 0) ; (time,lon)
     eofvar2bp  = wgt_runave_n_Wrap (eofvar2, wgt, 0, 0)
     eofvar3bp  = wgt_runave_n_Wrap (eofvar3, wgt, 0, 0)

   ;------------------------------------------------
   ; remove temporal means of band pass series: *not* necessary 
   ;------------------------------------------------
     eofvar1bp  = dim_rmvmean_n(eofvar1bp, 0)              ; (time,lon)
     eofvar2bp  = dim_rmvmean_n(eofvar2bp, 0)
     eofvar3bp  = dim_rmvmean_n(eofvar3bp, 0)

   ;------------------------------------------------
   ; Compute the temporal variance at each lon
   ;------------------------------------------------
     vari_eofvar1 = dim_variance_n_Wrap(eofvar1bp, 0)     ; (lon)
     vari_eofvar2 = dim_variance_n_Wrap(eofvar2bp, 0)
     vari_eofvar3 = dim_variance_n_Wrap(eofvar3bp, 0)

   ;------------------------------------------------
   ; Compute the zonal mean of the temporal variance
   ;------------------------------------------------
     zavg_vari_eofvar1 = dim_avg_n_Wrap( vari_eofvar1 , 0)     
     zavg_vari_eofvar2 = dim_avg_n_Wrap( vari_eofvar2, 0)
     zavg_vari_eofvar3 = dim_avg_n_Wrap( vari_eofvar3, 0)

   ;------------------------------------------------
   ; Normalize by sqrt(avg_var*)
   ;------------------------------------------------
     eofvar1bp  = eofvar1bp/sqrt(zavg_vari_eofvar1)          ; (time,lon)
     eofvar2bp  = eofvar2bp/sqrt(zavg_vari_eofvar2)
     eofvar3bp  = eofvar3bp/sqrt(zavg_vari_eofvar3)

   ;------------------------------------------------
   ; Combine the normalized data into one variable
   ;------------------------------------------------
     cdata     = new ( (/3*mlon,ntim/), typeof(eofvar1bp), getFillValue(eofvar1bp))
     do ml=0,mlon-1
       cdata(ml       ,:) = (/ eofvar1bp (:,ml) /)
       cdata(ml+  mlon,:) = (/ eofvar2bp (:,ml) /)
       cdata(ml+2*mlon,:) = (/ eofvar3bp(:,ml) /)
     end do

   ;------------------------------------------------
   ; Compute **combined** EOF; Sign of EOF is arbitrary
   ;------------------------------------------------
     eof_cdata    = eofunc(cdata   , neof, False)      ; (neof,3*mlon)
    print("==============")
    printVarSummary(eof_cdata)
    printMinMax(eof_cdata, True)

    eof_ts_cdata = eofunc_ts(cdata,eof_cdata,False)   ; (neof,3*ntim)
    print("==============")                                  
    printVarSummary(eof_ts_cdata)
    printMinMax(eof_ts_cdata, True)

   ;------------------------------------------------
   ; For clarity, explicitly extract each variable. Create time series 
   ;------------------------------------------------
 
    nvar = 3  ; "var1", "var2", "var3"
    ceof = new( (/nvar,neof,mlon/), typeof(cdata), getFillValue(cdata))

    do n=0,neof-1
      ceof(0,n,:) = eof_cdata(n,0:mlon-1)      ; var1
      ceof(1,n,:) = eof_cdata(n,mlon:2*mlon-1) ; var2
      ceof(2,n,:) = eof_cdata(n,2*mlon:)       ; var3
    end do

    ceof!0   = "var"
    ceof!1   = "eof"
    ceof!2   = "lon"   
    ceof&lon =  olr&lon

    ceof_ts        = new( (/nvar,neof,ntim/), typeof(cdata), getFillValue(cdata))
    ceof_ts(0,:,:) = eofunc_ts_Wrap(eofvar1bp(lon|:,time|:),ceof(0,:,:),False)   ; (0,neof,ntim)
    ceof_ts(1,:,:) = eofunc_ts_Wrap(eofvar2bp(lon|:,time|:),ceof(1,:,:),False)   ; (1,neof,ntim)
    ceof_ts(2,:,:) = eofunc_ts_Wrap(eofvar3bp(lon|:,time|:),ceof(2,:,:),False)   ; (2,neof,ntim)
    printVarSummary(ceof_ts)       ; (nvar,neof,time)

   ;------------------------------------------------
   ; Add code contributed by Marcus N. Morgan, Florida Institute of Technology; Feb 2015
   ; Calculate % variance (pcv_ )accounted for by OLR, U850 and U200
   ;------------------------------------------------

    pcv_eof_var1 = new(neof,typeof(ceof))
    pcv_eof_var2 = new(neof,typeof(ceof))
    pcv_eof_var3 = new(neof,typeof(ceof))
      
    do n=0,neof-1
       pcv_eof_var1(n) = avg((ceof(0,n,:)*sqrt(ceof@eval(n)))^2)*100
       pcv_eof_var2(n) = avg((ceof(1,n,:)*sqrt(ceof@eval(n)))^2)*100
       pcv_eof_var3(n) = avg((ceof(2,n,:)*sqrt(ceof@eval(n)))^2)*100
       print("pcv: neof="+(n+1)+":  "+pcv_eof_var1(n)+"  "+pcv_eof_var2(n)+"  "+pcv_eof_var3(n))
    end do

   ;------------------------------------------------
   ; Change sign? Ad hoc approach [ Eun-Pa Lim]
   ;------------------------------------------------
                                          ; index (subscript)
    imax_var2   = maxind(ceof(1,0,:))     ; max u850 eof 1 
    imin_var2   = minind(ceof(1,0,:))     ; min     

    imax_var3   = maxind(ceof(1,1,:))     ; max u200 eof 2
    imin_var3   = minind(ceof(1,1,:))     ; min

    lonmax_var2 = ceof&lon(imax_var2)      ; longitude of max value
    lonmax_var3 = ceof&lon(imax_var3)

    if (lonmax_var2.lt.50 .or. lonmax_var2.gt.120) then
      ceof(:,0,:)       = -ceof(:,0,:)         
      ceof_ts(:,0,:)    = -ceof_ts(:,0,:)  
      eof_cdata(0,:)    = -eof_cdata(0,:)
      eof_ts_cdata(0,:) = -eof_ts_cdata(0,:)
    end if

    if (lonmax_var3.lt.90 .or. lonmax_var3.gt.150) then
      ceof(:,1,:)       = -ceof(:,1,:)         
      ceof_ts(:,1,:)    = -ceof_ts(:,1,:)  
      eof_cdata(1,:)    = -eof_cdata(1,:)
      eof_ts_cdata(1,:) = -eof_ts_cdata(1,:)
    end if

    print("==============")
    printVarSummary(eof_cdata)
    printMinMax(eof_cdata, True)

   ;------------------------------------------------
   ; Compute cross correlation of each variable's EOF time series at zero-lag
   ;------------------------------------------------
    r_var1_var2  = escorc(ceof_ts(0,:,:) , ceof_ts(1,:,:) )   ; (neof)
    r_var1_var3  = escorc(ceof_ts(0,:,:) , ceof_ts(2,:,:) )
    r_var2_var3 = escorc(ceof_ts(1,:,:) , ceof_ts(2,:,:) )

    print("==============")
    do n=0,neof-1
      print("neof="+n \
          +"  r_var1_var2=" +sprintf("%4.3f",r_var1_var2(n))  \
          +"  r_var1_var3=" +sprintf("%4.3f",r_var1_var3(n))  \
          +"  r_var2_var3=" +sprintf("%4.3f",r_var2_var3(n)) ) 
    end do
    print("==============")

   ;------------------------------------------------
   ; Compute cross correlation of the multivariate EOF; EOF 1 vs EOF 2
   ;------------------------------------------------

    mxlag     = floattoint(26*spd)
    rlag_01   = esccr(eof_ts_cdata(0,:),eof_ts_cdata(1,:), mxlag)   ; (N,mxlag+1)
    rlag_10   = esccr(eof_ts_cdata(1,:),eof_ts_cdata(0,:), mxlag)   ; (N,mxlag+1)
    ccr_12    = new ( (/2*mxlag+1/), float)    

    ccr_12(mxlag:)    = rlag_10(0:mxlag)   
    ccr_12(0:mxlag)   = rlag_01(::-1)       ; reverse order

   ;------------------------------------------------
   ; Normalize the multivariate EOF 1&2 component time series
   ; Compute (PC1^2+PC2^2): values > 1 indicate "strong" periods
   ;------------------------------------------------
    eof_ts_cdata(0,:) = eof_ts_cdata(0,:)/stddev(eof_ts_cdata(0,:))
    eof_ts_cdata(1,:) = eof_ts_cdata(1,:)/stddev(eof_ts_cdata(1,:))

    mjo_ts_index      = eof_ts_cdata(0,:)^2 + eof_ts_cdata(1,:)^2 
    mjo_ts_index_smt  = runave(mjo_ts_index, floattoint(91*spd), 0) ; 91-day running mean

    nGood   = num(.not.ismissing(mjo_ts_index))     ; # non-missing
    nStrong = num(mjo_ts_index .ge. 1.0)
    print("nGood="+nGood+"   nStrong="+nStrong+"   nOther="+(nGood-nStrong))

   ;------------------------------------------------
   ; Write PC results to netCDF for use in another example.
   ;------------------------------------------------
    mjo_ts_index!0    = "time"
    mjo_ts_index&time = time 
    mjo_ts_index@long_name = "MJO PC INDEX" 
    mjo_ts_index@info      = "(PC1^2 + PC2^2)" 

    PC1  = eof_ts_cdata(0,:)
    PC1!0= "time"
    PC1&time      =  time
    PC1@long_name = "PC1"
    PC1@info      = "PC1/stddev(PC1)"

    PC2  = eof_ts_cdata(1,:)
    PC2!0= "time"
    PC2&time      =  time
    PC2@long_name = "PC2"
    PC2@info      = "PC2/stddev(PC2)"

    
    outfile = datainfo+"_MJO_PC_INDEX.nc"
    system("/bin/rm -f "+outpath+outfile)   ; remove any pre-existing file
    ncdf = addfile(outpath+outfile,"c")     ; open output netCDF file
                                            ; make time an UNLIMITED dimension 
    filedimdef(ncdf,"time",-1,True)         ; recommended  for most applications
                                            ; output variables directly
    ncdf->MJO_INDEX = mjo_ts_index    
    ncdf->PC1       = PC1     
    ncdf->PC2       = PC2     

   ;------------------------------------------------
   ; PLOTS
   ;------------------------------------------------
    wksNameBase=plotdir+datainfo
    yyyymmdd = cd_calendar(time, -2)
    yrfrac   = yyyymmdd_to_yyyyfrac(yyyymmdd, 0.0)

    day      = ispan(-mxlag, mxlag, 1)

    wks1 = gsn_open_wks(wksType,wksNameBase+"_multi_eof")
    wks2 = gsn_open_wks(wksType,wksNameBase+"_ccr_eof")
    wks3 = gsn_open_wks(wksType,wksNameBase+"_mjo_index")
    print(wks1)
    plot = new(3,graphic)                

   ;------------------------------------------------
   ; Multivariate EOF plots
   ;------------------------------------------------
                                      
    rts           = True
    rts@gsnDraw   = False       ; don't draw yet
    rts@gsnFrame  = False       ; don't advance frame yet
    rts@gsnScale  = True        ; force text scaling               

    rts@vpHeightF = 0.40        ; Changes the aspect ratio
    rts@vpWidthF  = 0.85
    rts@vpXF      = 0.10        ; change start locations
    rts@vpYF      = 0.75        ; the plot
    rts@xyLineThicknesses = (/2, 2, 2/)
    rts@xyLineColors      = (/"black","red","green"/)
    rts@gsnYRefLine            = 0.                  ; reference line   

    rts@pmLegendDisplayMode    = "Always"            ; turn on legend
    rts@pmLegendSide           = "Top"               ; Change location of 
    rts@pmLegendParallelPosF   = 0.86                ; move units right
    rts@pmLegendOrthogonalPosF = -0.50               ; move units down
    rts@pmLegendWidthF         = 0.15                ; Change width and
    rts@pmLegendHeightF        = 0.15                ; height of legend.
    rts@lgLabelFontHeightF     = 0.0175

    rtsP                       = True                ; modify the panel plot
    rtsP@gsnMaximize           = True                ; large format
    rtsP@txString              = datainfo + " Multivariate EOF: 15S-15N: "+yrStrt+"-"+yrLast 
  
    do n=0,neof-1
      rts@xyExplicitLegendLabels = (/u200name+": "+sprintf("%4.1f", pcv_eof_var3(n)) +"%" \
                                    ,u850name+": "+sprintf("%4.1f", pcv_eof_var2(n))+"%" \
                                    , olrname+": "+sprintf("%4.1f", pcv_eof_var1(n))+"%" /)
      rts@gsnLeftString  = "EOF "+(n+1)
      rts@gsnRightString = sprintf("%3.1f",ceof@pcvar(n))  +"%"
      plot(n) = gsn_csm_xy (wks1,lon,ceof(:,n,:),rts)
    end do
    gsn_panel(wks1,plot(0:1),(/2,1/),rtsP)     ; now draw as one plot

   ;------------------------------------------------
   ; cross correlation plots; delete unneeded resources 
   ;------------------------------------------------
    delete([/rts@xyExplicitLegendLabels, rts@pmLegendDisplayMode, rts@xyLineThicknesses \
            ,rts@gsnLeftString, rts@gsnRightString, rts@xyLineColors/] )
  
    lag              = ispan(-mxlag,mxlag,1)
    lag@long_name    = "lag (days)"

    plot(0)          = gsn_csm_xy (wks2, lag ,ccr_12,rts)
    rtsP@txString    = datainfo + " Cross Correlation: Multivariate EOF: 15S-15N: " \
                   +  yrStrt+"-"+yrLast 
    rtsP@gsnPaperOrientation = "portrait"        ; force portrait
    gsn_panel(wks2,plot(0),(/1,1/),rtsP)     ; now draw as one plot

   ;------------------------------------------------
   ; MJO "strong" index 
   ;------------------------------------------------
    rts@gsnYRefLine        = 1.0
    rts@gsnYRefLineColor   = "black"
    rts@xyMonoDashPattern  = True
    rts@xyLineColors       = (/"black", "blue"/)
    rts@xyLineThicknesses  = (/1, 2/)
    rts@pmLegendDisplayMode    = "Always"            ; turn on legend
    rts@pmLegendWidthF         = 0.12                ; Change width and
    rts@pmLegendHeightF        = 0.10                ; height of legend.
    rts@pmLegendParallelPosF   = 0.86                ; move units right
    rts@pmLegendOrthogonalPosF = -0.40               ; move units down
    rts@xyExplicitLegendLabels = (/"daily", "91-day runavg" /)

    mjo_ind_plt = new ( (/2,ntim/), typeof(mjo_ts_index))
    mjo_ind_plt(0,:) = mjo_ts_index
    mjo_ind_plt(1,:) = (/ mjo_ts_index_smt /)
    plot(0) = gsn_csm_xy(wks3, yrfrac,mjo_ind_plt,rts)

    rtsP@txString  = datainfo + " MJO Index: (PC1^2+ PC2^2) : 15S-15N: "+yrStrt+"-"+yrLast 
    gsn_panel(wks3,plot(0),(/1,1/),rtsP)     ; now draw as one plot

   ;------------------------------------------------
   ; delete unnecessary variables
   ;------------------------------------------------
    delete([/eofvar1,eofvar2,eofvar3,eofvar1bp,eofvar2bp,eofvar3bp, \
            vari_eofvar1,vari_eofvar2,vari_eofvar3, \
            zavg_vari_eofvar1,zavg_vari_eofvar2,zavg_vari_eofvar3, \
            cdata,eof_cdata,eof_ts_cdata,ceof,ceof_ts, \
            pcv_eof_var1, pcv_eof_var2 , pcv_eof_var3, \
            r_var1_var2, r_var2_var3, r_var1_var3, \
            rlag_01, rlag_10/])
 

;###########################################################
;
; Create MJO Composite Life Cycles
;
;###########################################################
    print("-------Create MJO Composite Life Cycles-------")
    if (isvar("wksNameBase")) then
      delete(wksNameBase)
    end if
    if (isvar("plot")) then
	  delete(plot)
    end if
    wksNameBase=datainfo+"mjo_composite_lifecycle"
    
   ;------------------------------------------------
   ; Place each array into an appropriate array
   ;------------------------------------------------

    phase       = (/ 202.5, 247.5, 292.5, 337.5 \   ; phi: Matt Wheeler's f77 code
                   ,  22.5,  67.5, 112.5, 157.5 /) 
    nPhase      = dimsizes(phase)
    phaBnd      = 180.0/nPhase
    angBnd      = new( (/2,nPhase/), "float")
    angBnd(0,:) = phase-phaBnd 
    angBnd(1,:) = phase+phaBnd
  
    r2d         = 180./(4.*atan(1.0))
    ang         = atan2(PC2,PC1)*r2d     ; phase space 
    nn          = ind(ang.lt.0)
    ang(nn)     = ang(nn) + 360          ; make 0 to 360


    nDays       = new (nPhase, "integer")
    pLabel      = "P"+ispan(1,nPhase,1)+": "



   ;------------------------------------------------------------
   ; PLOTS
   ;------------------------------------------------------------
    pltPath = plotdir+wksNameBase
    plotSubTitle = "Anomalous: OLR, U850, V850"
    wks  = gsn_open_wks(wksType,pltPath)
    gsn_define_colormap(wks,"ViBlGrWhYeOrRe") 
    plot = new(nPhase,graphic)              ; create graphic array

    res                      = True         
    res@gsnDraw              = False        ; don't draw yet
    res@gsnFrame             = False        ; don't advance frame yet
    res@gsnSpreadColors      = True         ; spread out color table
 
    res@mpFillOn             = False        ; turn off map fill
    res@mpMinLatF            = -20          ; zoom in on map
    res@mpMaxLatF            = 20
    res@mpCenterLonF         = 210.
    res@cnFillOn             = True         ; turn on color fill
    res@cnLinesOn            = False        ; True is default
    res@cnLineLabelsOn       = False        ; True is default
    res@lbLabelBarOn         = False        ; turn off individual lb's
    res@gsnScalarContour     = True         ; contour 3rd array 
    res@gsnMajorLatSpacing   = 15
    res@gsnMajorLonSpacing   = 60
    res@tmXBLabelFontHeightF = 0.01
    res@tmYLLabelFontHeightF = 0.01

    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF       =  -40         ; -100; mnmxint(0)
    res@cnMaxLevelValF       =   40         ;   80; mnmxint(1)
    res@cnLevelSpacingF      =    5         ;   20; mnmxint(2)

    res@vcMinDistanceF            = 0.01            ; thin the vector density
    res@vcRefMagnitudeF           = 2.0             ; define vector ref mag
    res@vcRefLengthF              = 0.025           ; define length of vec ref
    res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
    res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
    res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

                                            ; panel plot only resources
    resP                     = True         ; modify the panel plot
    resP@gsnMaximize         = True         ; large format
    resP@gsnPanelLabelBar    = True         ; add common colorbar
    resP@lbLabelAutoStride   = True         ; auto stride on labels
    resP@lbLabelFontHeightF  = 0.01
    resP@gsnPanelBottom      = 0.05         ; add some space at bottom
    resP@pmLabelBarWidthF    = 0.8          ; label bar width
    resP@pmLabelBarHeightF   = 0.05
    resP@gsnPanelFigureStringsFontHeightF = 0.0125  ; bit larger than default
 
    txres                  = True             
    txres@txFontHeightF    = 0.01
    txid = gsn_create_text(wks, plotSubTitle, txres)

    amres                  = True 
    amres@amOrthogonalPosF =  0.75  
    amres@amJust           = "CenterCenter"

   ;------------------------------------------------------------
   ; Loop over each phase
   ;------------------------------------------------------------

    res@gsnLeftString   = ""
    res@gsnRightString  = ""
    do nSeason=1,2
      if (nSeason.eq.1) then
        resP@txString  = datainfo + " " +yrStrt+"-"+yrLast+": May to Oct"
      else
        resP@txString  = datainfo + " " +yrStrt+"-"+yrLast+": Nov to Apr"
      end if
      do n=0,nPhase-1
     
       ; All times where the following criteria are met
        if (nSeason.eq.1) then
          nt = ind(mjo_ts_index.gt.1.0    .and.                     \
             (imon.ge.5          .and. imon.le.10).and.    \
             ang.ge.angBnd(0,n) .and. ang.lt.angBnd(1,n))
        else
          nt = ind(mjo_ts_index.gt.1.0    .and.                     \
             (imon.ge.11         .or.  imon.le. 4).and.    \
             ang.ge.angBnd(0,n) .and. ang.lt.angBnd(1,n))
        end if
        if (.not.all(ismissing(nt))) then
          olrAvg  = dim_avg_n_Wrap(  olr(nt,:,:), 0)           ; better
          u850Avg = dim_avg_n_Wrap( u850(nt,:,:), 0)
          v850Avg = dim_avg_n_Wrap( v850(nt,:,:), 0)

          nDays(n) = dimsizes(nt)
          res@tmXBLabelsOn         = False    ; do not draw lon labels
          res@tmXBOn               = False    ;             lon tickmarks
          if (n.eq.(nPhase-1)) then           ; 
            res@tmXBLabelsOn     = True       ; draw lon labels
            res@tmXBOn           = True       ;          tickmarks
          end if

          plot(n) = gsn_csm_vector_scalar_map_ce(wks,u850Avg,v850Avg,olrAvg,res)
        end if
        delete(nt)                  ; will change next iteration
      end do
    resP@gsnPanelFigureStrings= pLabel+nDays
    gsn_panel(wks,plot,(/nPhase,1/),resP)     ; now draw as one plot
    end do

  end if
end



