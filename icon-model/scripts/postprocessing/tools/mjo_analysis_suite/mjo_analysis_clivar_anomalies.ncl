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
;   mjo_analysis_clivar_anomalies.ncl
; Purpose:
;   Creates daily anomalies from daily means (wind) or daily sums (precip, olr)
;   using both raw and smoothed climatology
; Origin:
;   This script is based on mjoclivar_2.ncl from the 
;   NCL MJO Climate Variability Tools
;-----------------------------------------------------------
; Required Input:
;   Time series of daily means, e.g. constructed from 
;   6 hourly output via cdo daymean
;   
; Program flow: 
;   - Read daily means of variable 
;     (preprocessed from 6 hourly to daily means via cdo)
;   - Compute daily climatology for specified period (2001-2005)
;     either raw or smoothed (retaining first two harmonics)
;   - Compute daily anomalies from raw or smoothed climatology
;   - Write daily anomalies to NetCDF output file
;   - Plot raw and smoothed daily climatology at various places
;   - Plot series of daily anomalies from raw or smoothed 
;     climatology at one particular place  
;    
;***********************************************************


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
print("started")
begin

;***********************************************************
; Preparation  
;***********************************************************

;------------------------------------------------
; Requires to receive the following variables
; from shell
;
;------ File and Data Info------
;      var        - variable
;      infile     - name of input file
;      outfile    - name of output file
;      plotfile   - name of plot file
;      datainfo   - experiment/analysis id
;      ymdStrt    - start date of time series
;      ymdLast    - end date of time series
;      lats       - latitudes to be
;                   considered
;      lons       - longitudes to be 
;                   considered
;      times      - number of time steps
;--------------------------------------------------    

 

; Define plot settings  
   wksType = "png"
   wksNameBase = plotfile+datainfo+var

; Specify number of fourier comp 
   nhar    = 4
   yrStrt  = ymdStrt/10000
   yrLast  = ymdLast/10000

;***********************************************************
; Read data and assign variables accordingly 
;-----------------------------------------------------------
; To allow flexibility in variable and coordinate names
; in input data, the script does not retrieve them via their
; names, but in dependence on their dimensions (assuming 
; that lat/lon/time differ in their dimension)
; Further, the name of the dimensions is unified to allow
; data from different origins to be used.
;***********************************************************

   fin = addfile(infile,"r")
   vars  = getfilevarnames(fin)                  ; write all variable names in array
   
   nvars = dimsizes(vars)
   timecount=0
   do i = 0,nvars-1
     dummy      = fin->$vars(i)$                 ; read variable based on array entry 
     dummydims  = getvardims(dummy)              ; get dimensions of variable
     ndims   = dimsizes(dummydims)               ; get number of dimensions
    
                                                 ; separate between data and coordinates
     if ( ndims .eq. 1) then                     ; data = several dimensions  
                                                 ; coordinates =  one dimension       
       print("--------Get info on dimensions--------")       
                                                 
       if (dimsizes(dummy) .eq. lats) then       ; identify lat, lon and time coords
          print("---get lats---")
          lat     = fin->$vars(i)$               ; set dimension names to standard 
          lat!0   = "lat"
          lat@lat = "latitude"
       end if
        
       if (dimsizes(dummy) .eq. lons) then
          print("---get lons---")
          lon     = fin->$vars(i)$
          lon!0   = "lon"
          lon@lon = "longitude"
       end if

       if ((dimsizes(dummy) .eq. times) .and. (timecount .eq. 0)) then
          timecount=1                            ; just read one time axis (if there are several)
          print("---get time---")
          time      = fin->$vars(i)$
          time!0    = "time"
          time@time = "time"
          time@units= "days since 2001-01-01 00:00:0.0"
       end if
     end if    
                                                 ; end of process for coordinates
     delete(dummy)
     if (ndims .gt. 1) then                      ; identify data and its characteristics
       print("---get data---")
       if (getfilevartypes(fin,vars(i)) .eq. "short") then
          datain = short2flt(fin->$vars(i)$ )
       else if (getfilevartypes(fin,vars(i)) .eq. "integer") then
          datain = int2flt(fin->$vars(i)$ )
       else 
          datain = fin->$vars(i)$
       end if
       end if                                     ; determine size of dimensions and 
       subdims = dimsizes(datain)                 ; rename dimensions accordingly
       do j = 0, ndims-1
           if (subdims(j) .eq. times) then
              datain!(j)  = "time"
              datain@time = "time"
           end if
           if (subdims(j) .eq. lats) then
              datain!(j)  = "lat"
              datain@lat  = "latitude"
           end if
           if (subdims(j) .eq. lons) then
              datain!(j)  = "lon"
              datain@lon  = "longitude"
           end if
       end do
     end if  
     delete(dummydims)
     delete(ndims)
   end do

; further specify timing
   ymd     = cd_calendar(time, -2)            ; YYYYMMDD integer output
   iStrt   = ind(ymd.eq.ymdStrt)              ; index start
   iLast   = ind(ymd.eq.ymdLast)              ; index last 
   delete(ymd)
   print(iStrt)
   print(iLast)
; get variable

   timeselection = time(iStrt:iLast)
   TIME          = cd_calendar(timeselection, 0)          ; type float 
   year          = floattointeger( TIME(:,0) )
   month         = floattointeger( TIME(:,1) )
   day           = floattointeger( TIME(:,2) ) 
   hour          = floattointeger( TIME(:,3) )
   ddd           = day_of_year(year, month, day) 
   yyyyddd       = year*1000 + ddd                        ; needed for input
   
   x       =  datain(iStrt:iLast,:,:)                     ; read data to x 
   printVarSummary( x ) 

;***********************************************************
; Compute raw daily climatology  
;***********************************************************
   xClmDay = clmDayTLL(x, yyyyddd)                        ; daily climatology 
                                                          ; at each grid point
   printVarSummary(xClmDay)            

;***********************************************************
; Compute smoothed daily climatology using 'nhar' harmonics
;***********************************************************
   xClmDay_sm = smthClmDayTLL(xClmDay, nhar)              ; daily smothed climatology
                                                          ; at each grid point
   printVarSummary(xClmDay_sm)

;***********************************************************
; Compute daily anomalies using raw and smoothed climatologies
;***********************************************************
    xAnom      = calcDayAnomTLL (x, yyyyddd, xClmDay)     
    printVarSummary(xAnom)
    printMinMax(xAnom, True)

    xAnom_sm   = calcDayAnomTLL (x, yyyyddd, xClmDay_sm)     
    xAnom_sm@long_name = "Anomalies from Smooth Daily Climatology"
    printVarSummary(xAnom_sm)
    printMinMax(xAnom_sm, True)

    delete( x )                                          ; no longer needed

;***********************************************************
; Create output file: daily anomalies
;***********************************************************

    dimx   = dimsizes(xAnom)
    ntim   = dimx(0)
    nlat   = dimx(1)
    nlon   = dimx(2)

    system("/bin/rm -f "+outfile)                        ; rm any pre-exist file, if any
    fnc    = addfile (outfile,"c")
    filAtt = 0
    filAtt@title         = var +": Daily Anomalies: "+yrStrt+"-"+yrLast+" "+datainfo  
    filAtt@source_file   = fin
    filAtt@creation_date = systemfunc("date")
    setfileoption(fnc,"DefineMode",True)
        
    varout      = var+"_anom"
    varout_sm   = var+"_anom_sm"

    dimNames = (/"time", "lat", "lon"/)  
	dimSizes = (/ -1 ,  nlat,  nlon/) 
	dimUnlim = (/ True , False, False/)   
    filedimdef(fnc,dimNames,dimSizes,dimUnlim)

    filevardef(fnc, "time"     ,typeof(time)  , getvardims(time)) 
    filevardef(fnc, "lat"      ,typeof(lat)   , getvardims(lat)) 
    filevardef(fnc, "lon"      ,typeof(lon)   , getvardims(lon))
    filevardef(fnc, varout     ,typeof(xAnom) , getvardims(xAnom))    
    filevardef(fnc, varout_sm  ,typeof(xAnom) , getvardims(xAnom))    
      
            
    filevarattdef(fnc, "time"    ,time)                      ; copy time attributes
    filevarattdef(fnc, "lat"     ,lat)                       ; copy lat attributes
    filevarattdef(fnc, "lon"     ,lon)                       ; copy lon attributes
    filevarattdef(fnc, varout    ,xAnom)                
    filevarattdef(fnc, varout_sm ,xAnom_sm)                

    fnc->time          = (/timeselection/)     
    fnc->lat           = (/lat/)
    fnc->lon           = (/lon/)
    fnc->$varout$      = (/xAnom   /)
    fnc->$varout_sm$   = (/xAnom_sm/)
   
;************************************************
; Plot results
;************************************************
    LAT   = (/ 50, 45,  5,  -5, -45, -50 /)
    LON   = (/270, 30, 90,  90, 180,  0 /)
    nPts  = dimsizes( LAT )

    plot  = new ( nPts, graphic)
    data  = new ( (/2,366/), typeof(xClmDay), getFillValue(xClmDay))

    wks   = gsn_open_wks (wksType,wksNameBase+"_daily_climatology")    

    res                   = True                      ; plot mods desired
    res@gsnDraw           = False
    res@gsnFrame          = False
    res@trXMinF           =   1
    res@trXMaxF           = 365
       
    res@xyLineThicknesses = (/1.0, 2.0/)              ; make 2nd lines thicker
    res@xyLineColors      = (/"blue","red"/)          ; change line color
    res@xyMonoDashPattern = True                      ; all solid

    do np=0,nPts-1 
       data(0,:) = xClmDay(:,{LAT(np)},{LON(np)})
       data(1,:) = xClmDay_sm(:,{LAT(np)},{LON(np)})
       res@gsnCenterString = "lat="+LAT(np)+"  lon="+LON(np)
       plot(np)  = gsn_csm_y (wks,data,res) ; create plot
    end do
       
    resP                     = True                ; modify the panel plot
    resP@txString            = var+": Raw/Smooth Daily Climatology: "+yrStrt+"-"+yrLast
    resP@gsnMaximize         = True
    resP@gsnPaperOrientation = "portrait"
    gsn_panel(wks,plot,(/(nPts/2),2/),resP)               ; now draw as one plot

    ;==========
    ; Plot anomalies for an arbitrarily selected near equatorial location
    ; Time as specified  
    ;==========
    LATX    = 0
    LONX    = 0

    yyyymmdd  = cd_calendar(timeselection, -2)


    xAnom@long_name    = "Anomalies from Raw"   ; short labels for plot
    xAnom_sm@long_name = "Anomalies from Smooth"
    wks   = gsn_open_wks (wksType,wksNameBase+"_daily_anomalies")
    ntBegin   = ind(yyyymmdd.eq.ymdStrt)
    ntEnd     = ind(yyyymmdd.eq.ymdLast)
       
    monthLabels    = (/1,4,7,10/)
    monNam = (/"Jan","Feb","Mar","Apr","May","Jun" \
                  ,"Jul","Aug","Sep","Oct","Nov","Dec" /)
    xVal   = new(ntim, typeof(xAnom&time) , "No_FillValue") ; bigger than
    xLab   = new(ntim, "string", "No_FillValue")        ; needed
    xValm  = new(ntim, typeof(xAnom&time) , "No_FillValue") ; bigger than

    ntm            = -1
    cr             = inttochar(10)                     ; carriage return
    do nt=ntBegin,ntEnd
      if (day(nt).eq.1) then
        ntm       = ntm + 1
        xVal(ntm) = xAnom&time(nt)
        xLab(ntm) = monNam(month(nt)-1)
        if (month(nt).eq.1) then
           xLab(ntm) = xLab(ntm) + " " + sprinti("%0.4i", year(nt))
        end if
      end if
    end do
    print(ntm)
    rxy  = True      
    rxy@gsnDraw     = False
    rxy@gsnFrame    = False
    rxy@gsnYRefLine = 0.0                    ; create a reference line   
    rxy@gsnAboveYRefLineColor = "red"        ; above ref line fill red
    rxy@gsnBelowYRefLineColor = "blue"       ; below ref line fill blue
    rxy@xyLineThicknessF  = 2.0                                               
    rxy@vpHeightF  = 0.4                     ; resize            
    rxy@vpWidthF   = 0.8                  
    rxy@tmXBMode   = "Explicit"
    rxy@tmXBValues = xVal(0:ntm:6)
    rxy@tmXBLabels = xLab(0:ntm:6)
    rxy@tmXBMinorValues = xVal(0:ntm)

    plot(0)  = gsn_csm_xy (wks,timeselection(ntBegin:ntEnd) \
                           ,xAnom(ntBegin:ntEnd,{LATX},{LONX}),rxy) 
    plot(1)  = gsn_csm_xy (wks,timeselection(ntBegin:ntEnd) \
                           ,xAnom_sm(ntBegin:ntEnd,{LATX},{LONX}),rxy) 
    resP@txString  = var+": Anomalies: (0N,0E)"
    gsn_panel(wks,plot(0:1),(/2,1/),resP)   

end
quit

