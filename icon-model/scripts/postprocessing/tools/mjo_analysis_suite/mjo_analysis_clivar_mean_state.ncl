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
;
;-----------------------------------------------------------
; Required Input:
;   Boreal winter or summer mean of 
;   * 850hPa Zonal Wind
;   * Daily Total Precipitation
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
;      infile1    - seasonal Tot Prec file
;      infile2    - seasonal U850hPa file
;      plotdir    - directory for plots
;      datainfo   - info on data name
;      season     - winter or summer season
;      lats       - number of latitudes
;      lons       - number of longitudes
;--------------------------------------------------    

;************************************************
; open file and read in data
;***********************************************
 
  if ( season .eq. "summer" )
    months = "(May - Oct)"
  end if
  if (season .eq. "winter" )
    months = "(Nov - Apr)"
  end if

; Define plot settings  
   wksType = "png"
   wksNameBase = plotdir+datainfo+"_mean_state_"+season
   wksPlotTitle = datainfo+" 850hPa Zonal Wind & Total Precipitation ~C~        mean state in boreal "+season+months
; Read data and assign variables accordingly
   fin1 = addfile(infile1,"r")
   vars1  = getfilevarnames(fin1)                  ; write all variable names in array
   fin2 = addfile(infile2,"r")
   vars2  = getfilevarnames(fin2)                  ; write all variable names in array
  
; Process file 1
   nvars = dimsizes(vars1)
   timecount=0
   do i = 0,nvars-1
     dummy      = fin1->$vars1(i)$               ; read variable 
     dummydims  = getvardims(dummy)              ; get dimensions of variable
     ndims      = dimsizes(dummydims)            ; get number of dimensions
     delete(dummydims)      
                                      ; separate between data and coordinates
     if ( ndims .eq. 1) then                      ; data = several dimensions        
                                                  ; coordinates =  one dimension 
       if (dimsizes(dummy) .eq. lats) then        ; identify lat, lon and time coords
          lat     = fin1->$vars1(i)$              ; set dimension names to standard 
          lat!0   = "lat"
          lat@lat = "latitude"
       end if        
       if (dimsizes(dummy) .eq. lons) then
          lon     = fin1->$vars1(i)$
          lon!0   = "lon"
          lon@lon = "longitude"
       end if
 
     end if                                        ; end of process for coordinates
     delete(dummy) 
     if (ndims .eq. 2) then                        ; identify data and its characteristics
       print("---get data file 1---")
       if (getfilevartypes(fin1,vars1(i)) .eq. "short") then
          datain1 = short2flt(fin1->$vars1(i)$ )
       else if (getfilevartypes(fin1,vars1(i)) .eq. "integer") then
          datain1 = int2flt(fin1->$vars1(i)$ )
       else 
          datain1 = fin1->$vars1(i)$
       end if
       end if  
                                                   ; determine size of dimensions and 
       subdims = dimsizes(datain1)                 ; rename dimensions accordingly
       do j = 0, ndims-1
           if (subdims(j) .eq. lats) then
              datain1!(j)  = "lat"
              datain1@lat  = "latitude"
           end if
           if (subdims(j) .eq. lons) then
              datain1!(j)  = "lon"
              datain1@lon  = "longitude"
           end if
       end do
       datain1@long_name = "Daily Precipitation"
       datain1@units = "mm day~S~-1"
     end if  
     delete(ndims)
 
   end do
   delete(subdims)
; Process file 2
   nvars = dimsizes(vars2)
   timecount=0
   do i = 0,nvars-1
     dummy      = fin2->$vars2(i)$                 ; read variable 
     dummydims  = getvardims(dummy)                ; get dimensions of variable
     ndims      = dimsizes(dummydims)              ; get number of dimensions
     delete(dummy)
     delete(dummydims)
     if (ndims .eq. 2) then                        ; identify data and its characteristics
       print("---get data file 2---")
       if (getfilevartypes(fin2,vars2(i)) .eq. "short") then
          datain2 = short2flt(fin2->$vars2(i)$ )
       else if (getfilevartypes(fin2,vars2(i)) .eq. "integer") then
          datain2 = int2flt(fin2->$vars2(i)$ )
       else 
          datain2 = fin2->$vars2(i)$
       end if
       end if                                      ; determine size of dimensions and 
       subdims = dimsizes(datain2)                 ; rename dimensions accordingly
       do j = 0, ndims-1
           if (subdims(j) .eq. lats) then
              datain2!(j)  = "lat"
              datain2@lat  = "latitude"
           end if
           if (subdims(j) .eq. lons) then
              datain2!(j)  = "lon"
              datain2@lon  = "longitude"
           end if
       end do
       datain2@long_name = "850hPa Zonal Wind"
       datain2@units = "m s~S~-1"
     end if  
     delete(ndims)
   end do  
     
     printVarSummary(datain1)
 
;************************************************
; create overlay plot
;************************************************
  wks = gsn_open_wks(wksType,wksNameBase)      ; send graphics to PNG file
; define options for map
  res                             = True     ; plot mods desired
  res@mpFillOn                    = False    ; turn on color fill
  res@mpMinLatF                   = -25.     ; range to zoom in on
  res@mpMaxLatF                   =  15.
  res@mpMinLonF                   =  40.
  res@mpMaxLonF                   = 320.
  res@mpCenterLonF 				  = 180.
  res@mpOutlineOn                 = True
  res@mpGeophysicalLineColor      = "blue"
  res@mpGeophysicalLineThicknessF = 2
  res@gsnDraw                     = False
  res@gsnFrame                    = False

; define options for first plot (shaded)
  res@cnLevelSelectionMode = "ExplicitLevels" ; use explicit levels
  res@cnLevels    = (/2.,3.5,5,6.5,8,9.5,11,12.5,14/)
  res@cnLineLabelsOn       = False            ; do not use line labels
  res@cnFillOn             = True             ; color fill
  res@cnLinesOn            = False            ; do not draw contour lines
  res@cnFillColors = (/"white","(/0.192,0.404,0.141/)","(/0.278,0.545,0.200/)", \
					   "(/0.365,0.686,0.255/)","(/0.996,0.886,0.388/)", \
					   "(/0.953,0.639,0.173/)","(/0.875,0.420,0.118/)", \
					   "(/0.765,0.204,0.082/)","(/0.592,0.086,0.078/)", \
               		   "(/0.408,0.039,0.067/)" /)  
  res@gsnAddCyclic          = False          ; data already has cyclic point
  res@gsnMajorLatSpacing = 10              ; change maj lat tm spacing
  res@gsnMajorLonSpacing = 40              ; change maj lon tm spacing
  res@tmXBLabelFontHeightF = 0.01          ; change maj lat tm spacing
  res@tmYLLabelFontHeightF = 0.01          ; change maj lat tm spacing

  res@tiMainString = wksPlotTitle            ; set the main title
  res@tiMainFontHeightF = 0.01
  res@gsnRightString = ""
  res@gsnLeftString = ""
  res@lbLabelFontHeightF = 0.01
  res@pmLabelBarOrthogonalPosF = 0.4 
; define options for second plot
  sres                      = True              ; set up a second resource list
  sres@gsnDraw              = False             ; do not draw the plot
  sres@gsnFrame             = False             ; do not advance the frame
  sres@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
  sres@cnLevels             = ispan(-9,9,3)     ; set the contour levels
  sres@gsnContourZeroLineThicknessF = 2.	    ; doubles thickness of zero contour
  sres@cnLineLabelsOn		= False
  sres@gsnContourNegLineDashPattern = 2
  sres@cnInfoLabelOn        = False
  sres@gsnRightString = ""
  sres@gsnLeftString = ""
 
  plot    = gsn_csm_contour_map(wks,datain1,res)   ; create the precipitation plot
  plot_ov = gsn_csm_contour(wks,datain2,sres)      ; create the U-wind plot

  overlay(plot,plot_ov)                       ; overlay the U-wind plot on the temperature plot
  draw(plot)                                  ; draw the temperature plot (with the U-wind plot overlaid)
  frame(wks)                                  ; advance the frame
end

  
 
