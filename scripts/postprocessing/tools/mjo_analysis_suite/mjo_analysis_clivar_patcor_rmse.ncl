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
;   mjo_analysis_clivar_mean_state.ncl  
; Purpose:
;   Plot mean state fields
;-----------------------------------------------------------
; Required Input:
;   Boreal winter or summer Mean of 
;   * All ICON Forecast
;   * Analysis/Observations
;   Both files have to have same lat, lon and time dimensions
;
;
;
;
;    
;***********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "/e/uhome/jukeller/mjo_tools/mjo_clivar_ncl/diagnostics_cam.ncl"

begin
   print("ok")
   lats_full=181
   lons_full=360
   lats_half=361
   lons_half=720
   lats_f_tr=101
   lats_h_tr=201

;************************************************
; open file and read in data
;***********************************************
   nana=2
   analysisfile_full="/e/uwork/jukeller/icon4mjo/data/TRMM/TRMM_TOT_PREC_nh_winter_fulldeg.nc"
   analysisfile_half="/e/uwork/jukeller/icon4mjo/data/TRMM/TRMM_TOT_PREC_nh_winter_halfdeg.nc"
   nexp=3
   expfile_ICON_80km="/e/uwork/jukeller/icon4mjo/data/clim_80km/new_tuning/processed/ICON_80km_TOT_PREC_nh_winter.grb2"
   expfile_ICON_40km="/e/uwork/jukeller/icon4mjo/data/clim_40km/new_tuning/processed/ICON_40km_TOT_PREC_nh_winter.grb2"
   expfile_ICON_20km="/e/uwork/jukeller/icon4mjo/data/clim_20km/new_tuning/processed/ICON_20km_TOT_PREC_nh_winter.grb2"

; Define plot settings  
;   wksType = "png"
;   wksNameBase = plotdir+resinfo+"mean_state_"+season
;   wksPlotTitle = " Pattern Correlation and normalized RMSE ~C~        mean state in boreal "+season
; Read analysis data and assign variables accordingly
   analysis = addfiles((/analysisfile_full,analysisfile_half/),"r")
   
; Read and process analysis data
   do nf=0,nana-1
     vars  = getfilevarnames(analysis[nf])
     nvars = dimsizes(vars)
     do i = 0,nvars-1
       dummy      = analysis[nf]->$vars(i)$            ; read variable 
       dummydims  = getvardims(dummy)              ; get dimensions of variable
       ndims      = dimsizes(dummydims)            ; get number of dimensions
       print(vars)
       print(ndims)
       delete(dummydims)      
                                                   ; separate between data and coordinates
       if ( ndims .eq. 1) then                     ; data = several dimensions        
                                                   ; coordinates =  one dimension 
         if ((dimsizes(dummy) .eq. lats_full) .or. (dimsizes(dummy) .eq. lats_half)) then 
           print("---get lats---")                 ; identify lat, lon and time coords
           lat     = analysis[nf]->$vars(i)$          ; set dimension names to standard 
           lat!0   = "lat"
           lat@lat = "latitude"
         end if        
         if ((dimsizes(dummy) .eq. lons_full) .or. (dimsizes(dummy) .eq. lons_half)) then
           print("---get lons---")
           lon     = analysis[nf]->$vars(i)$
           lon!0   = "lon"
           lon@lon = "longitude"
         end if
       end if                                        ; end of process for coordinates

       delete(dummy) 
 
       if (ndims .eq. 2) then                        ; identify data and its characteristics
         print("---get data---")
         if (getfilevartypes(analysis[nf],vars(i)) .eq. "short") then
           analysisdata[nf] = short2flt(analysis[nf]->$vars(i)$ )
         else if (getfilevartypes(analysis[nf],vars(i)) .eq. "integer") then
           analysisdata[nf] = int2flt(analysis[nf]->$vars(i)$ )
         else 
           analysisdata[nf] = analysis[nf]->$vars(i)$
         end if
         end if                                        ; determine size of dimensions and 
         subdims = dimsizes(analysisdata[nf])              ; rename dimensions accordingly
         do j = 0, ndims-1
           if ((subdims(j) .eq. lats_full).or.(subdims(j) .eq. lats_half)) then
             analysisdata[nf]!(j)  = "lat"
             analysisdata[nf]@lat  = "latitude"
           end if
           if ((subdims(j) .eq. lons_full).or.(subdims(j) .eq. lons_half)) then
             analysisdata[nf]!(j)  = "lon"
             analysisdata[nf]@lon  = "longitude"
           end if
         end do
       end if  
       delete(ndims)
     end do
   end do



; Read and process ICON data
; In the same loop, do also compute RMSE and pattern correlation and store in array
   do nf=0,nexp-1
     vars  = getfilevarnames(experiments[nf])
     nvars = dimsizes(vars)
     do i = 0,nvars-1
       dummy      = experiments[nf]->$vars(i)$            ; read variable 
       dummydims  = getvardims(dummy)              ; get dimensions of variable
       ndims      = dimsizes(dummydims)            ; get number of dimensions
       print(vars)
       print(ndims)
       delete(dummydims)      
                                                   ; separate between data and coordinates
       if ( ndims .eq. 1) then                     ; data = several dimensions        
                                                   ; coordinates =  one dimension 
         if ((dimsizes(dummy) .eq. lats_full) .or. (dimsizes(dummy) .eq. lats_half)) then  
           print("---get lats---")
           lat     = experiments->$vars(i)$          ; set dimension names to standard 
           lat!0   = "lat"
           lat@lat = "latitude"
         end if        
         if ((dimsizes(dummy) .eq. lons_full) .or. (dimsizes(dummy) .eq. lons_half))then
           print("---get lons---")
           lonf     = experiments->$vars(i)$
           lon!0   = "lon"
           lon@lon = "longitude"
         end if
       end if                                        ; end of process for coordinates
       delete(dummy) 
 
       if (ndims .eq. 2) then                        ; identify data and its characteristics
         print("---get data---")
         if (getfilevartypes(experiment[nf],vars(i)) .eq. "short") then
           expdata = short2flt(experiment[nf]->$vars(i)$ )
         else if (getfilevartypes(experiment[nf],vars(i)) .eq. "integer") then
           expdata = int2flt(experiment[nf]->$vars(i)$ )
         else 
           expdata = experiment[nf]->$vars(i)$
         end if
         end if                                        ; determine size of dimensions and 
         subdims = dimsizes(expdata)              ; rename dimensions accordingly
         do j = 0, ndims-1
           if ((subdims(j) .eq. lats_full).or.(subdims(j) .eq. lats_half)) then
             expdata!(j)  = "lat"
             expdata@lat  = "latitude"
           end if
           if ((subdims(j) .eq. lons_full).or.(subdims(j) .eq. lons_half)) then
             expdata!(j)  = "lon"
             expdata@lon  = "longitude"
           end if
         end do
       end if  
       delete(ndims)
     end do
;  compute weighted area RMSE
     if (nf .lt. 2) then
       rmse[nf]=wgt_arearmse2(experiment[nf],analysisdata[0],wgt,1)
     else if (nf .eq. 2) then
       rmse[nf]=wgt_arearmse2(experiment[nf],analysisdata[1],wgt,1)
     end if
     end if

;  compute pattern correlation
     if (nf .lt. 2) then
       patcor[nf]=pattern_cor(experiment[nf],analysisdata[0],1.0,1)
     else if (nf .eq. 2) then
       patcor[nf]=pattern_cor(experiment[nf],analysisdata[1],1.0,1)
     end if
     end if


   end do
end

 
