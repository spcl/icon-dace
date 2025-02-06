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
;   mjo_analysis_clivar_filtering.ncl  
; Purpose:
;   Applies bandpass filter to daily anomalies
; Origin:
;   This script is a merger of NCL Scripts mjoclivar_3.ncl
;   mjoclivar_4.ncl mjoclivar_5.ncl, mjoclivar6.ncl, 
;   mjoclivar_7.ncl and mjoclivar_10.ncl
;-----------------------------------------------------------
; Required Input:
;   Daily Anomalies of variable to be processed
;   Created by mjo_analysis_clivar_anomalies.ncl
;   
; Program flow: 
;   - Read daily anomalies
;   - Compute unfiltered and filtered seasonal variances 
;     and their relation, and plot them (mjoclivar_3)
;   - Apply bandpass filter and plot results as
;     * Area-Time Plot (mjoclivar_4)
;     * Hovmöller Plot (mjoclivar_5)
;     * Lat-lon map    (mjoclivar_6)
;   - Compute and plot spectra over region (mjoclivar_7)
;   - Wavenumber Frequency Spectra over season (mjoclivar_10) 
;     [This one is not tested for 48h data! Ok for 24h and less]
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
;      dayint   - interval input data in days
;------ Definition for Bandpass ----
;      bpmin    - bandpass filter minimum
;      bpmax    - bandpass filter maximum
;      bpwgt    - bandpass filter weight
;      yperhov  - no of years plotted in hov
;------ Area to be considered ----
;      latmin
;      latmax     
;      lonmin    
;      lonmax  
;------ Switches for various plots (on/off)------
;      vari     - create variance plot
;      bptime   - create bp-filtered time series
;      bphov    - create bp-filtered hovmoeller diagram
;      bparea   - create bp-filtered map
;      spectra  - create spectres 
;--------------------------------------------------    


; Define band pass
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

;###########################################################
;
; Get data  
;
;###########################################################

print("-------Get data and prepare-------")
  


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
      x = short2flt( f->$varname$(iStrt:iLast,:,:) ) 
   else
      x = f->$varname$(iStrt:iLast,:,:) 
   end if
; extract limited area for Hovmöller Plots 
   x_area = x(:,{latmin:latmax},{lonmin:lonmax})  ; explicitly extract  

   printVarSummary(x)
   printMinMax(x, True)

;  process date information and create start and end date for hovmoeller diagrams
   time = f->time(iStrt:iLast)                ; days since ...
   ymdhms         = cd_calendar(time, 0 )     ; (ntim,2)
   mm             = floattoint( ymdhms(:,1) ) ; just months
   NoOfTimeSteps  = dimsizes( ymdhms(:,0) )-1
   NoHovplots=floattoint(NoOfTimeSteps/(yperhov*365/dayint))
   
   HovInterval=floattoint(NoOfTimeSteps/NoHovplots)
   hStrt= new(NoHovplots,integer)
   hLast= new(NoHovplots,integer)
   hovStrt= new(NoHovplots,integer)
   hovLast= new(NoHovplots,integer)
   do i = 0,NoHovplots-1
     hStrt(i)= i*HovInterval
     hLast(i)= hStrt(i)+HovInterval
     hovStrt(i)=cd_calendar(time(hStrt(i)),-2)
     hovLast(i)=cd_calendar(time(hLast(i)),-2)
   end do
   
  

   
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
   
;################################################
;
; Create Variance Plots
;
;################################################

  if ( vari .eq. "true" ) then
    print("-------Create Variance Plots-------") 
    if (isvar("wksNameBase")) then
      delete(wksNameBase)
    end if

  ;------------------------------------------------
  ; Compute unfiltered seasonal variances (lat,lon)
  ;------------------------------------------------

    print("---unfiltered seasonal variances---")
    xVarWinter = dim_variance_Wrap( x(lat|:,lon|:,time|iWinter) ) 
    xVarSummer = dim_variance_Wrap( x(lat|:,lon|:,time|iSummer) ) 
    xVarWinter_nozero=where(xVarWinter.ne.0,xVarWinter,xVarWinter@_FillValue)
    xVarSummer_nozero=where(xVarSummer.ne.0,xVarSummer,xVarSummer@_FillValue)
    printMinMax(xVarWinter,True)
    printMinMax(xVarSummer,True)
    ;printMinMax(xVarWinter_nozero,True)
    ;printMinMax(xVarSummer_nozero,True)

 
  ;------------------------------------------------
  ; Compute filtered seasonal variances (lat,lon)
  ;------------------------------------------------

    print("---filtered seasonal variances---")
    xVarWinterBPF = dim_variance_Wrap( xBPF(:,:,iWinter) ) 
    xVarSummerBPF = dim_variance_Wrap( xBPF(:,:,iSummer) ) 
    ;printMinMax(xVarWinterBPF,True)
    ;printMinMax(xVarSummerBPF,True)

  ;------------------------------------------------
  ; Compute (filtered_variance)/(unfiltered_variance)
  ;------------------------------------------------
    print("---fraction filt/unfilt variances---")
    xRatioVarWinter = (xVarWinterBPF/xVarWinter_nozero)*100
    xRatioVarSummer = (xVarSummerBPF/xVarSummer_nozero)*100
    ;xRatioVarWinter = where(ismissing(xVarWinter_nozero),0,xVarWinter_nozero)
    ;xRatioVarSummer = where(ismissing(xVarSummer_nozero),0,xVarSummer_nozero)
 

    copy_VarCoords( xVarWinter, xRatioVarWinter)
    copy_VarCoords( xVarSummer, xRatioVarSummer)


  ;------------------------------------------------
  ; Plots
  ;------------------------------------------------

    if (wksType.eq."png") then
      wksTypeLocal = "eps"
    else
      wksTypeLocal = wksType
    end if
  
    wksNameBase = plotdir+datainfo+"_"+var+"_variance"
    plot = new ( 2, "graphic")

  ;------------------------------------------------
  ; set minimum and maximum contour interval
  ; in dependence on the variable to plot
  ; values are based on observation/reanalysis data
  ;------------------------------------------------
    if (var .eq. "U200") then
      conmin =  20
      conmax = 300
      conint =  20
    end if
    if (var .eq. "U850") then
      conmin =  10
      conmax = 130
      conint =  10
    end if
    if (var .eq. "V850") then
      conmin =  5
      conmax = 70
      conint =  5
    end if
    if ((var .eq. "OLR").or.(var .eq. "TOT_PREC")) then
      conmin =  200
      conmax = 1800
      conint =  200
    end if
 


  ;------------------------------------------------
  ; resource list for first (color) data array
  ;------------------------------------------------
    res1                      = True
    res1@gsnDraw              = False              ; don't draw
    res1@gsnFrame             = False              ; don't advance frame
    res1@gsnStringFontHeightF = 0.0125             ; make larger than default
    res1@lbLabelBarOn         = False              ; turn off individual cb's
    res1@gsnSpreadColors      = True               ; use full range of colors
    res1@mpCenterLonF         = 180.               ; center plot at 180
    res1@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
    res1@cnMinLevelValF       =  10.
    res1@cnMaxLevelValF       =  50.
    res1@cnLevelSpacingF      =   5. 

    res1@mpMinLonF            = lonmin               ; select a subregion
    res1@mpMaxLonF            = lonmax
    res1@mpMinLatF            = latmin 
    res1@mpMaxLatF            = latmax

   ;res1@mpLandFillColor      = "background"       ; color of land
   ;res1@cnFillDrawOrder      = "Predraw"         

  ;------------------------------------------------
  ; resource list for second (contour only) data array
  ;------------------------------------------------
    res2                      = True
    res2@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
    res2@cnMinLevelValF       =  conmin
    res2@cnMaxLevelValF       =  conmax
    res2@cnLevelSpacingF      =  conint
    ;res2@cnInfoLabelOn        = False             ; do not want
    res2@cnInfoLabelOrthogonalPosF = -0.11        ; move the label inside th plot
    res2@cnInfoLabelFontHeightF    = 0.008
    res2@cnInfoLabelString    =  "$CMN$ TO $CMX$ BY $CIU$"  
    res2@cnLineLabelsOn       = True
    res2@cnLabelMasking       = True		
    res2@cnLineLabelFontHeightF    = 0.008


  ;------------------------------------------------
  ; resource list for panel
  ;------------------------------------------------
    resP                     = True                ; modify the panel plot
    resP@gsnMaximize         = True                ; make large
    resP@gsnPanelLabelBar    = True                ; add common colorbar
    resP@lbLabelFontHeightF  = 0.0125

  ;------------------------------------------------
  ; plot data
  ;------------------------------------------------

    wks     = gsn_open_wks(wksTypeLocal, wksNameBase)        
    gsn_define_colormap(wks,"amwg")

    delete(xVarWinter@long_name)
    res1@gsnLeftString        = "Winter (N-A)"
    plot(0)   = gsn_csm_contour_map_overlay(wks,xRatioVarWinter,xVarWinter,res1,res2)

    delete(xVarSummer@long_name)
    res1@gsnLeftString        = "Summer (M-O)"
    plot(1)   = gsn_csm_contour_map_overlay(wks,xRatioVarSummer,xVarSummer,res1,res2)
    resP@txString  = datainfo + " Variance: "+x@long_name +" "+ var + " "+ yrStrt+"-"+yrLast

    gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot


    if (wksType.eq."png") then
      delete(wks)
      print("convert "+wksNameBase+".eps "+wksNameBase+".png")
      system("convert -rotate 270 "+wksNameBase+".eps "+wksNameBase+".png")
      system("/bin/rm -f "+wksNameBase+".eps")
    end if
  end if


;################################################
;
; Create Bandpass Area-Time Result
;
;################################################
 
  if (bptime .eq. "true") then
    print("-------Create BPF Area Time Plots-------") 
    if (isvar("wksNameBase")) then
      delete(wksNameBase)
    end if


  ;------------------------------------------------
  ; Extract Area averages and filtered series 
  ; Use NCL's coordinate subscripting
  ;------------------------------------------------

    wgty    = latRegWgt(x_area&lat, "double", 0)
    optArea = False
    xAreaT  = band_pass_area_time(x_area, spd, bpf, wgty, optArea)

    printVarSummary(xAreaT)              ; (time)
    printMinMax(xAreaT, True)
 
  ;------------------------------------------------
  ; Plot  band_pass_area_time result 
  ;------------------------------------------------
  
    wksNameBase = datainfo+"_"+var+"_bp_timeseries"
    optArea = True
    optArea@pltTitle = datainfo + " "+ xAreaT@long_name+": Area ["+latmin+":"+latmax+" , "+lonmin+":"+lonmax+"]"
    band_pass_area_time_plot(xAreaT, time, plotdir, wksTypeLocal ,wksNameBase, optArea)

    if (wksType.eq."png") then
      if (isvar("wks")) then
          delete(wks)
      end if
      print("convert "+plotdir+wksNameBase+".eps "+plotdir+wksNameBase+".png")
      system("convert -rotate 270 "+plotdir+wksNameBase+".eps "+plotdir+wksNameBase+".png")
      system("/bin/rm -f "+plotdir+wksNameBase+".eps")
    end if

  end if

;################################################
;
; Create Bandpass Hovmöller Result
;
;################################################

  if (bphov .eq. "true") then
    print("-------Create BPF Hovmoeller Plots-------") 
    i=0
    do i=0,NoHovplots-1
      if (isvar("wksNameBase")) then
        delete(wksNameBase)
      end if
      wksNameBase = datainfo+"_"+var+"_bp_hovmoeller"+hovStrt(i)+"-"+hovLast(i)

      wgty    = latRegWgt(x_area&lat, "double", 0)

      optHov  = False
      hov     = band_pass_hovmueller (x_area, spd, bpf, wgty, optHov) 
      printVarSummary(hov)               ; (time,lon)
      printMinMax(hov, True)

  ;------------------------------------------------
  ; Create generic Hovmueller plot using the
  ; output from band_pass_hovmueller.
  ; This is for the specific time period of interest 
  ;------------------------------------------------
      
      if (var.eq."OLR") then
        optHovP = True
        optHovP@cnLevelSpacingF =   4
        optHovP@cnMinLevelValF  = -24
        optHovP@cnMaxLevelValF  =  24
      end if
      if (var.eq."TOT_PREC") then
        optHovP = True
        optHovP@cnLevelSpacingF =   0.5
        optHovP@cnMinLevelValF  = -5
        optHovP@cnMaxLevelValF  =  5
      end if
      if (var.eq. "U200" .or. var.eq. "V850" .or.  var.eq. "U850" ) then
        optHovP = True
        optHovP@cnLevelSpacingF =   1
        optHovP@cnMinLevelValF  = -8
        optHovP@cnMaxLevelValF  =  8
      end if
      band_pass_hovmueller_plot(hov(hStrt(i):hLast(i),:) \
                           ,plotdir, wksTypeLocal ,wksNameBase, datainfo, optHovP)  

      if (wksType.eq."png") then
        if (isvar ("wks")) then
           delete(wks)
        end if
        system("convert "+plotdir+wksNameBase+".eps "+plotdir+wksNameBase+".png")
        system("/bin/rm -f "+plotdir+wksNameBase+".eps")
      end if
    end do
  end if

;################################################
;
; Create Bandpass Lat-Lon-Map Result
;
;################################################

  if (bparea .eq. "true") then
    print("-------Create BPF Map Plots-------") 
    if (isvar("wksNameBase")) then
      delete(wksNameBase)
    end if
  ;------------------------------------------------
  ; Extract filtered series at each grid point
  ;------------------------------------------------

    optLatLon  = False
    mjo_latlon = band_pass_latlon_time(x_area, spd, bpf, optLatLon)
    printVarSummary(mjo_latlon)
    printMinMax(mjo_latlon, True)

  ;------------------------------------------------
  ; create spatial plot at one time
  ;------------------------------------------------
    plot = new ( 2, "graphic")

    res                       = True     ; plot mods desired
    res@gsnDraw               = False    ; don't draw
    res@gsnFrame              = False    ; don't advance frame

    res@cnFillOn              = True     ; turn on color fill
    res@cnLinesOn             = False    ; turn of contour lines
    res@cnLineLabelsOn        = False    ; turn of contour line labels
    res@gsnSpreadColors       = True     ; use full range of color map
   ;res@lbLabelBarOn          = False    ; turn off individual cb's
    res@lbLabelAutoStride     = True     ; let NCL figure spacing
    res@lbOrientation         = "horizontal" ; vertical label bar

    res@mpMinLatF             = latmin     ; range to zoom in on
    res@mpMaxLatF             = latmax
    res@mpMinLonF             = lonmin
    res@mpMaxLonF             = lonmax
    res@mpCenterLonF          = 180      ; DL is middle
   ;res@pmTickMarkDisplayMode = "Always" ; use default lat/lon labels

                                       ; nice symmetric min&max for region
    rSym    = True
    symMinMaxPlt (mjo_latlon,20,False,rSym)
    mnmxint = nice_mnmxintvl( min(x_area) , max(x_area), 16, False)

    rOrig   = True
    rOrig@cnLevelSelectionMode = "ManualLevels"
    rOrig@cnMinLevelValF       = mnmxint(0)
    rOrig@cnMaxLevelValF       = mnmxint(1)
    rOrig@cnLevelSpacingF      = mnmxint(2)

    resP                       = True           ; modify the panel plot
    resP@gsnMaximize           = True
    resP@gsnPaperOrientation   = "portrait"

    res@gsnAddCyclic           = False          ; {lonL:lonR}
 
    do nt=hStrt,hLast,2                         ; all times
      wksNameBase = plotdir+datainfo+"_"+var+"_bp_latlonmap"+date(nt)
      wks  = gsn_open_wks(wksTypeLocal,wksNameBase)
      gsn_merge_colormaps(wks,"amwg_blueyellowred","BlAqGrYeOrReVi200")
      res@gsnSpreadColorStart  = 18
      res@gsnSpreadColorEnd    = -1
      res                      = rOrig
      print("Regular Anomalies")

      plot(0) = gsn_csm_contour_map_ce(wks,x_area(nt,:,:), res)

      res@gsnSpreadColorStart  =  2
      res@gsnSpreadColorEnd    = 17
      res                      = rSym
      print("Bandpass Anomalies")
      plot(1) = gsn_csm_contour_map_ce(wks,mjo_latlon(nt,:,:), res)

      resP@txString            = datainfo +" "+date(nt)
      gsn_panel(wks,plot,(/2,1/),resP)

      if (wksType.eq."png") then
         delete(wks)
         system("convert "+wksNameBase+".eps "+wksNameBase+".png")
         system("/bin/rm -f "+wksNameBase+".eps")
    end if

    end do
  end if

;################################################
;
; Create Spectra over Regions
;
;################################################
  
  if (spectra .eq. "true") then
    print("-------Create Spectra over Regions-------") 
    if (isvar("wksNameBase")) then
      delete(wksNameBase)
    end if

  
        if ((var .eq. "OLR") .or. (var .eq. "TOT_PREC")) then
          regionName       = (/ "IO-Winter" , "WP-Winter" , "MC-Winter" ,"IO-Summer" , "WP-Summer" , "BB-Summer" /)
          regionLatS       = (/ -10. , -20. , -18. , -10. ,  10. ,  10./)    
          regionLatN       = (/   5. ,  -5. ,  -3. ,   5. ,  25. ,  20./)
          regionLonW       = (/  75. , 160. , 115. ,  75. , 115. ,  80./)
          regionLonE       = (/ 100. , 185. , 145. , 100. , 140. , 100./)
        end if

        if ((var .eq. "U850") .or. (var .eq. "V850")) then
          regionName       = (/ "IO-Winter" , "WP-Winter" ,"IO-Summer" , "WP-Summer" , "EP-Summer" /)
          regionLatS       = (/ -16. , -14. ,  4. ,   4. ,   6./)    
          regionLatN       = (/  -1. ,   1. , 21. ,  21. ,  16./)
          regionLonW       = (/  69. , 164. , 69. , 119. , 241./)
          regionLonE       = (/  96. , 191. , 96. , 146. , 266./)
        end if

        if (var .eq. "U200") then
          regionName       = (/ "IO-Winter" , "WP-Winter" , "EP-Winter" ,  "IO-Summer" , "WP-Summer" , "EP-Summer" /)
          regionLatS       = (/   4. ,   4. , -16. , -16. ,   4. , -16./)    
          regionLatN       = (/  21. ,  21. ,   1. ,   1. ,  21. ,   1./)
          regionLonW       = (/  56. , 124. , 256. ,  44. , 123. , 239./)
          regionLonE       = (/  79. , 151. , 279. ,  71. , 151. , 266./)
        end if

      nRegion = dimsizes(regionName)
      nBounds = dimsizes(regionLatS)

      if (nRegion.ne.nBounds) then
        print("mjo_spectra: nRegion and nBounds do not match")
        print("             nRegion="+nRegion )
        print("             nBounds="+nBounds )
        exit
      end if



  ;------------------------------------------------
  ; MJO CLivar specify cos(lat)
  ;------------------------------------------------

      weighty  = cos(0.01745329*lat)       
      copy_VarCoords(lat, weighty)

  ;------------------------------------------------
  ; driver to calculate and plot the spectra
  ;------------------------------------------------
      wksNameBase = datainfo+"_"+var+"_time_spectra"


      optMJO= True
      mjo_spectra (x, date, weighty, regionLatS, regionLatN, regionLonW, regionLonE, regionName \
              ,plotdir, wksTypeLocal, wksNameBase, optMJO,spd)
      if (wksType.eq."png") then
         plots=dimsizes(regionName)
         do j=0,plots-1
           system("convert "+plotdir+wksNameBase+"_"+regionName(j)+".eps "+plotdir+wksNameBase+"_"+regionName(j)+".png")
           system("/bin/rm -f "+plotdir+wksNameBase+"_"+regionName(j)+".eps")
         end do
      end if
  end if
;################################################
;
; Create Wavenumber-Frequency Spectra per Season
;
;################################################

  if (wavfreq .eq. "true") then
   print("-------Create Wav-Freq Spectra per Season-------") 
   if (isvar("wksNameBase")) then
      delete(wksNameBase)
    end if
  
    X = dim_avg_n_Wrap( x(:,{-10:10},:), 1 )                   ; (time,lon) ; v5.1.1 onward
  
  ;------------------------------------------------
  ; calculate and plot the OLR spectra
  ;------------------------------------------------

    optWavFrq     = False
    optPlot       = True
    optPlot@cnLinesOn    = False

    nameSeason = (/"winter" , "summer"/)

    do ns=0,dimsizes(nameSeason)-1
      wksNameBase = datainfo+"_"+var+"_"+nameSeason(ns)+"_wn_freq_spectra" 
      wf      = mjo_wavenum_freq_season (X, date,nameSeason(ns), spd,optWavFrq)  
      optPlot@tiMainString = nameSeason(ns)+": "+ymdStrt+"-"+ymdLast
      mjo_wavenum_freq_season_plot (wf, nameSeason(ns), plotdir, wksType, wksNameBase, optPlot) 
    end do
  end if

end



