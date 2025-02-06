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
;   mjo_analysis_clivar_wk_spacetime.ncl  
; Purpose:
;   Create Wheeler-Kiladis Wavenumber-Frequency Spectra plot
; Origin:
;   This script is based on the NCL Scripts wkSpaceTime_1-3.ncl
;   available at www.ncl.ucar.edu/Applications/space_time.shtml
;-----------------------------------------------------------
; Required Input:
;   Time series of daily or 6/12hourly means, e.g. constructed from 
;   6 hourly output via cdo daymean (not the processed anomalies!)
; Program Flow:
;   Read data and plot Space-Time Spectra using the NCL Function
;   wkSpaceTime
;
;   Note:
;   NCL-Function "wkSpaceTime" assigns colortables 
;   in dependence on variable names -> commit required names 
;   to get same color scheme for different experiments!
;   Variable Names are "FLUT", "OLR", "olr", "U200", "U850",
;   "OMEGA500" and "PRECT"
;
; Oct 2015, DWD, Julia Keller, FE13, julia.keller@dwd.de

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "diagnostics_cam.ncl"
;    
;***********************************************************

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

; Set plotting options
  opt = True
  opt@debug = True
  opt@pltType = "png"
  opt@cnLinesOn = True ; turn off contour lines

print("-------Get data and prepare-------")
;***********************************************************
; Load file and process information on variables etc   
;   - Get info on variables
;   - Identify which variable contains coordinates or data,
;     based on dimensions and sizes
;   - Bring dimension and coordinate names to standard                
;***********************************************************

; Read data and assign variables accordingly
   fin   = addfile (infile, "r") 

   vars  = getfilevarnames(fin)                  ; write all variable names in array
   
   nvars = dimsizes(vars)
   timecount=0
   do i = 0,nvars-1
     print(vars)
     print(nvars)
     dummy      = fin->$vars(i)$                 ; read variable 
     dummydims  = getvardims(dummy)           ; get dimensions of variable
     ndims   	= dimsizes(dummydims)              ; get number of dimensions
     delete(dummy)
                                                   ; separate data and coordinates
     if ( ndims .eq. 1) then                     ; data = several dimensions        
       print("--------Get info on dimensions--------")
                                                 ; coordinates =  one dimension 
       if (ndims .eq. lats) then   ; identify lat, lon and time coords
         print("---get lats---")
         lat     = fin->$vars(i)$                ; set dimension names to standard 
         lat!0   = "lat"
         lat@lat = "latitude"
       end if
        
       if (ndims .eq. lons) then
         print("---get lons---")
         lon     = fin->$vars(i)$
         lon!0   = "lon"
         lon@lon = "longitude"
       end if
    
       if ((ndims .eq. times) .and. (timecount .eq. 0)) then
          timecount=1
          print("---get time---")
          time      = fin->$vars(i)$
          time!0    = "time"
          time@time = "time"
          time@units="hours since 2001-01-01 00:00:0.0"
       end if
       
     end if                                ; end of process for coordinates
     
     if (ndims .ge. 3) then                        ; identify data characteristics
      print("---get data---")
       if (getfilevartypes(fin,vars(i)) .eq. "short") then
          datain = short2flt(fin->$vars(i)$ )
       else if (getfilevartypes(fin,vars(i)) .eq. "integer") then
          datain = int2flt(fin->$vars(i)$ )
       else 
          datain = fin->$vars(i)$
       end if
       end if
; determine size of dimensions and 
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
   print("ok")
   print("--------Extract subregion--------")

   if (lons .gt. 720) then
      xin = datain(:,{latmin:latmax:2},:)                       ; read data to x 
   else
      xin = datain(:,{latmin:latmax},:)                       ; read data to x 
   end if
 
   delete(datain)
   datatype=typeof(xin)
  
   if (datatype .eq. "short") then
     x=short2flt(xin)
   else
     x=xin
   end if


  

   wkSpaceTime(x,plotfile,datainfo,varname,latmax,floattoint(nSpD),nDayWin,nDaySkip,opt)


