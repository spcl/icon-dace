#!/bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#

# plot T,S and potential density variation to initial values from a list of ICON input files

#  Usage: script expPath - where expPath is output path of experiment to be plotted
#         scrip  

#  Updates Ralf
#   - ncl (icon_plot) with dozens of errors on blizzard, but okay on workstation
#   - cdo -r cat necessary - done

#  Updates Stephan:
#   - use all timesteps available -> removes bias at decades/eof-file
#   - use new variable names t_acc, s_acc, rhopot_acc of 'nml'-output
#   - use ncl-script developed by Stephan Lorenz and Karin Meier-Fleischer
#   - an init file is generated using initial output - currently rhopot is still zero (Oct 2013)
#   - now also for plotting of final data only - output identifier must be provided, no paths
#   - latest change - 2013-12-06


# ==============================================================================
set -ex
#  Tools and their paths
     CDO="cdo"                                   #  climate data operators
ICONPLOT=/pool/data/ICON/tools/icon_plot.ncl     #  NCL-script
 ICONLIB=/pool/data/ICON/tools

# ==============================================================================
# Input argument $1 is path of experiment
          expPath=${1:-=xxx}                     # import experiment path
[[ "$expPath" == "y" ]] && plotOnly="y"          # plot data via $1=y: outputIdent must be set accordingly, see below
# ==============================================================================
#        expIdent='Data.r13694.hupw.vppm'                                            # experiment identifier (typed)
#         expPath='/scratch/mpi/CC/mh0287/users/m211032/Icon/icon-dev.new/experiments' # experiment path (typed)
         expIdent=${expPath##*/}    # delete path on left before '/'                 # automatic identifier
      outputIdent=${expIdent#*.}    # delete smallest match on left before '.'       # automatic suffix
     fileListPath="$expPath/Output"                                                  # output data path
     fileListPath="$expPath"                                                         # output data path
# these lines must not be used if standards are matched; edit if plotOnly==y or other needs
if [[ "$plotOnly" == "y" ]]; then 
      outputIdent='r15085.qio1'                                                      # output file name appendix
fi
# fileListPattern="${expIdent}_iconR2B04-ocean_etopo40_planet_000[1-5].nc"           # 'nml' naming convention
  fileListPattern="${expIdent}_R2B04_oce_DOM01_ML_000[1-3].nc"                       # 'vlist' naming convention
   outputDataFile="ano.TSrho.$outputIdent.nc"                                        # output data file name
          Tempvar='T'                                                                # temperature variable name
           Salvar='S'                                                                # salinity variable name
           Rhovar='rhopot'                                                           # density variable name
          Tempvar='t_acc'                                                            # temperature variable name
           Salvar='s_acc'                                                            # salinity variable name
           Rhovar='rhopot_acc'                                                       # density variable name
       PlotScript='none'              #  no plots
       PlotScript='icon'              #  script icon_plot_ncl
       PlotScript='ncl'               #  ncl-script included, see below
# ==============================================================================
# these variable names must not necessarily be changed
      MaskVarName='wet_c'
          ScrDatN='tsrho'
          ScrDatN='scr_tsrho'
# ==============================================================================
#declare a fileListArray
if [[ "$plotOnly" == "y" ]]; then 
  echo " generate plot tsrho - file ano.TSrho.$outputIdent.nc must be available"
else
fileList=$(ls $fileListPath/$fileListPattern)
i=0
for file in $fileList; do
  fileListArray[$i]=$file
  i=$((i+1))
done
numberOfFiles=$i
echo "$numberOfFiles number of files"
echo $fileList
# ==============================================================================
# handling of the initial values
# get the initial values before any averaging is done
#initFile='init-phc.r2b4.20L.ts_acc.nc'     #  orig PHC T/S data
# include rhopot in PHC data:
# cdo chname,rhopoto,rhopot_acc -rhopot -adisit -chname,t_acc,tho,s_acc,sao init-phc.r2b4.20L.ts_acc.nc xx
#  - now use first data set of output file for subtracting init!
#    otherwise sequence, names, dimensions etc. are different when using PHC data
#  - caution: first set of accumulated values is zero
#initFile='init-phc.rho.r2b4.20L.ts_acc.nc'  #  contains names as in ocean output: t_acc, s_acc, rhopot_acc
initFile=fldmean_init.$outputIdent.nc        #  to be created from first input data file

# ==============================================================================
# using a mask file - contains no missing values
#maskFile='wetc.r2b4.r11xxx.nc'    # easy to generate - cdo selvar,wet_c
maskFile="wetc.$expIdent.nc"      # current mask created automatically
# ==============================================================================
# Loop over all files in serial
for file in $fileList; do
  fnum=${file##*_}  #  '_' must be last part before "0001.nc" of ocean output filename
  pattern_num=${outputIdent}_$fnum
  if [[ -f fldmean_${pattern_num} ]]; then
    echo  "fldmean_${pattern_num} already exists!"
  else
    if [[ -f ${ScrDatN}_${pattern_num} ]]; then
      echo  "${ScrDatN}_${pattern_num} already exists!"
    else
      if [[ ! -f "$maskFile" ]]; then 
        # create maskFile
        $CDO selvar,$MaskVarName $file $maskFile
      fi
      # select variables, mask out land points - takes some time
      $CDO -div -selname,$Tempvar,$Salvar,$Rhovar $file $maskFile ${ScrDatN}_${pattern_num}
    fi
    # compute the fldmean using masked values only
    $CDO fldmean ${ScrDatN}_${pattern_num} fldmean_${pattern_num}

    # special treatment of first file "*0001.nc": contains zero!
    #  - for _acc accumulated variables only
 #  if [[ $Tempvar == "t_acc" ]]; then 
 #    if [[ $fnum == "0001.nc" ]]; then 
 #      # no deldate available?
 #      # seldate,2001-01-11,2010-12-30 for 10-days output only!
 #      $CDO seldate,2001-01-11,2010-12-30 fldmean_${pattern_num} fldmean.noinit.nc
 #      mv fldmean.noinit.nc fldmean_${pattern_num}
 #    fi  # file no=1
 #  fi  # tempvar=t_acc
  fi  # create fldmean_${pattern_num}
  if [[ $fnum == "0001.nc" ]]; then 
    if [[ -f $initFile ]]; then
      echo  "$initFile already exists!"
    else
      # use first=initial timestep for hovmoeller anomaly plot:
      $CDO seltimestep,1 fldmean_${pattern_num} $initFile
    fi
  fi
done

# ==============================================================================
# Cat the files together - results would be appended to last file, delete it
[[ -f fldmean_$outputIdent.nc ]] && rm fldmean_$outputIdent.nc
$CDO -r cat fldmean_${outputIdent}_* fldmean_$outputIdent.nc

# ==============================================================================
# Subtract initial values
[[ -f $outputDataFile ]] || $CDO yearmean -sub fldmean_$outputIdent.nc $initFile $outputDataFile
  
fi  #  plotOnly
# ==============================================================================
# Plot a hovmoeller type graph using icon_plot.ncl
#  reports lots of errors, difficult to put on one page (psnup -3 if.ps of.ps)
if [[ $PlotScript == "icon" ]]; then 
  for varname in $Tempvar $Salvar $Rhovar; do 
    nclsh $ICONPLOT  -altLibDir=$ICONLIB -varName=$varname -iFile=$outputDataFile \
    -oFile=hov.ano.${outputIdent}_$varname -oType=ps -isIcon -hov=true
  done
fi

# ==============================================================================
# Plot a hovmoeller type graph using ncl-script
if [[ $PlotScript == "ncl" ]]; then 

cat >scr_plot_tsrho.ncl << EOF
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;-----------------------------------------------------------------------------
;-- NCL needs monotonous increasing time steps
;-- if time steps are missing NCL draws no time axis
;-- New arrays for time and data are created and filled with missing values
;-- Here: annual data
;--                                      Karin Meier-Fleischer, DKRZ, 18.09.13
;-----------------------------------------------------------------------------
;-- Use cdo to show the yearly time steps and get missing years or incorrect/
;-- uncomplete years.
;-- for yearly data ! mon=12=Dec
;-----------------------------------------------------------------------------
;-- more TBD:
;    - geometry of the plot (postscript) - do not use automatic scaling (gsn_panel)
;    - independent placing of string with $outputIdent - not in title
;    - shorter minor x-tick marks - do not use stride for x-axis labelling
;-----------------------------------------------------------------------------

   dates       = systemfunc("cdo showdate $outputDataFile")
   dates_count = str_fields_count(dates," ")
   dates_str   = new(dates_count,"string")
   yearlist    = new(dates_count,"integer")
   monlist     = new(dates_count,"integer")
   indmis      = new(dates_count,"integer")
   missing     = new(dates_count,"integer")

   j=0
   do i=0,dates_count-1
      dates_str = str_get_field(dates,i+1,"  ")
      yearlist(i) = tointeger(str_get_field(dates_str(i),1,"-"))
      if (i.ne.0) then
         last = last + 1
;         print(last + "     "+ yearlist(i))
         if (yearlist(i).ne.last) then
            m = j
            n = i
            do while(yearlist(i).ne.last)
               indmis(m) = n
               missing(m) = tointeger(last)
               last = last + 1
               m = m + 1
               n = n + 1
            end do
            j=m
         end if
      else
         last = yearlist(0)
      end if
      monlist(i) = tointeger(str_get_field(dates_str(i),2,"-"))
   end do

;-- get the indices of the wrong years (yearly data not complete)
   print("")
   print("Missing time values: ")
   do i=0,dimsizes(missing)-1
      if (.not.ismissing(missing(i))) then
         print(""+missing(i))
      end if
   end do

;-- get the indices of the missing years
   indlist = ind(monlist.ne.12)       ;-- for yearly data month must be 12

   do i=0,dimsizes(indlist)-1
      if (.not.ismissing(indlist(i))) then
         indlast = indlist(i)
      end if
   end do

;-- how many indices have incorrect or missing time steps
   nmisy = dimsizes(ind(.not.ismissing(missing)))
   nnoty = dimsizes(ind(.not.ismissing(indlist)))

   print("")
   print("Number of missing years:    "+nmisy)
   print("Number of incomplete years: "+nnoty)
   print("")

;-----------------------------------------------------------------------------
;-- read file original data file containing uncomplete years and/or missing years
;-----------------------------------------------------------------------------
   f     = addfile ("$outputDataFile", "r")

   temp  = f->t_acc(depth|:,time|:,lat|0,lon|0)
   s     = f->s_acc(depth|:,time|:,lat|0,lon|0)
   rho   = f->rhopot_acc(depth|:,time|:,lat|0,lon|0)

   time  = f->time
   depth = f->depth

;-- assign depth array for correct y-axis labeling
   ndepth    = dimsizes(depth)
   dlabels   = new(ndepth,string)
   dlabels   = depth

;-- set _FillValue to wrong time steps of original data
   temp(depth|:,time|indlist) = temp@_FillValue
   s(depth|:,time|indlist)    = s@_FillValue
   rho(depth|:,time|indlist)  = rho@_FillValue

;-- create new time dimension with continuing time steps
   tdim1 = dimsizes(time)
   tdim2 = tdim1  + nmisy
   tdim11 = tdim1-1
   if((tdim11).eq.indlast) then
      print("Last time step incorrect")
      tdim2 = tdim2-1                        ; --> KMF do not increase tdim2 if last time step
      lastyear = tointeger(time(dimsizes(time)-2)) + 365
      time2 = int2dble(ispan(tointeger(time(0)), lastyear, 365))
   else
      time2 = int2dble(ispan(tointeger(time(0)), tointeger(time(dimsizes(time)-1)),365))
   end if
   time2@standard_name = time@standard_name
   time2@units = time@units

;-- create new data arrays containing _FillValue at missing/wrong years
   tempdim2 = (/ndepth, tdim2/)

   temp2 = new(tempdim2, typeof(temp),temp@_FillValue)
   s2    = new(tempdim2, typeof(s),s@_FillValue)
   rho2  = new(tempdim2, typeof(rho),rho@_FillValue)

   tdim21 = tdim2-1

;-- wrong time steps are already set to _FillValue, we need to copy
;-- the data of the non-missing time steps from the original data
   k=0
   m=0
   n=0
   do l=0,tdim21
      if(.not.ismissing(indmis(m)) .and. indmis(m).eq.l) then
         temp2(:,k) = temp2@_FillValue
         s2(:,k)    = s2@_FillValue
         rho2(:,k)  = rho2@_FillValue
         k = k + 1
         m = m + 1
      else
         temp2(:,k) = temp(:,n)
         s2(:,k)    = s(:,n)
         rho2(:,k)  = rho(:,n)
         k = k + 1
         n = n + 1
      end if
   end do

;-- convert time to date strings
   utc_date  =  cd_calendar(time2, 0)
   year      =  utc_date(:,0)
   month     =  utc_date(:,1)

   print("")

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------

;-- open workstation
   wks = gsn_open_wks("ps","hov.ano.${outputIdent}")

;-- color maps
;  gsn_define_colormap(wks,"BlueWhiteOrangeRed")
   gsn_define_colormap(wks,"ViBlGrWhYeOrRe")

;-- set resources
   res                      =  True               ;-- set resource object
   res@gsnDraw              =  False              ;-- don't draw plot, yet
   res@gsnFrame             =  False              ;-- don't advance frame, yet
   res@gsnSpreadColors      =  True               ;-- subset of the whole colormap

;  res@vpWidthF             =  0.98               ;-- set view port width of each plot in panel
;  res@vpHeightF            =  0.40               ;-- set view port height of each plot in panel 
   res@vpWidthF             =  0.6                ;-- set view port width of each plot in panel
   res@vpHeightF            =  0.31               ;-- set view port height of each plot in panel 
                                                  ;--          (optimized value for "ps" output)
;  res@vpXF                 =  1.0
;  res@vpYF                 =  1.0
                                                  ;--          (optimized value for "ps" output)
   res@sfXArray             =  time2              ;-- uses time as plot x-axis
   res@sfYArray             =  depth              ;-- uses depth as plot y-axis

   res@trYReverse           =  True               ;-- reverse the y-axis, if requested
;   res@gsnYAxisIrregular2Linear =  True          ;-- converts irreg depth to linear, if requested

;-- common contour line settings
   res@cnFillOn             =  True               ;-- color for contours
   res@cnLineLabelsOn       =  True               ;-- draw line labels
   res@cnInfoLabelOn        =  False              ;-- don't draw the info label box
   res@cnLineLabelDensityF  =  1.0                ;-- density of labels <1.0 = less, >1.0 = more
   res@cnLineLabelFontHeightF = 0.01              ;-- set line label font size
   res@lbLabelStride        =  2                  ;-- labelling interval at color bar and contour lines

   res@tiMainOn             =  False              ;-- don't draw any title string

;-- common x-axis labeling settings
;  res@tmXBLabelsOn         =  False              ;-- False: x-axis labels removed
;  res@tmXBMode             = "Explicit"          ;-- set x-axis labeling to explicit
;  res@tmXBValues           =  time2              ;-- values for x-axis tickmarks
;  res@tmXBLabels           =  ""+year            ;-- set labels equal to values (type string)
;  res@tmXBLabelStride      =  5                  ;-- draw every 5th label
;  res@tmXBLabelFontHeightF =  0.018              ;-- x-axis font size
;  res@tmXBLabelAngleF      =  40                 ;-- rotate the x-axis labels counter clockwise
;  res@tmXBLabelDeltaF      =  1.0                ;-- move the x-axis labels downward
   res@tmXTOn               =  False              ;-- don't draw tickmarks on top of x-axis
   res@tmXBOn               =  False              ;-- don't draw tickmarks on bottom of x-axis
   res@tmXTLabelsOn         =  False              ;-- don't draw labels on top of x-axis
   res@tmXBLabelsOn         =  False              ;-- don't draw labels on bottom of x-axis
   res@tiXAxisOn            =  False              ;-- don't draw x-axis title        

;-- common y-axis labeling settings
;  res@tmYLMode             = "Explicit"          ;-- set y-axis labeling to explicit
;  res@tmYLValues           =  depth              ;-- values for y-axis tickmarks
;  res@tmYLLabels           =  ""+depth           ;-- set labels equal to values (type string)
;  res@gsnYAxisIrregular2Linear =  True           ;-- converts irreg depth to linear, if requested
;  res@tmYLLabelStride      =  0                  ;-- draw every 5th label
   res@tiYAxisString        = "Depth (m)"         ;-- y-axis title string
   res@tiYAxisFontHeightF   =  0.013              ;-- y-axis font size
   res@tmYLLabelFontHeightF =  0.011              ;-- y-axis mark font size

;-- common labelbar settings
;  res@lbBoxMinorExtentF    =  0.25               ;-- decrease the width of the labelbar
   res@lbLabelFontHeightF   =  0.012              ;-- label bar font
   res@lbOrientation        = "vertical"          ;-- set labelbar orientation
   res@pmLabelBarWidthF     =  0.13               ;-- set labelbar width
   res@pmLabelBarHeightF    =  0.34               ;-- set labelbar height
   ;  this increases the height of the 3 panels, fits DINA4 better!?

;-- assign plot array
   plot = new(3,graphic)

;-------------------------------------------------------------------------------
;-- upper plot - temp (temperature)
;-------------------------------------------------------------------------------
   res0                      =  res               ;-- use res
;  res0@cnFillPalette        = "BlWhRe"           ;-- set colormap
   res0@cnLevelSelectionMode = "ManualLevels"     ;-- set manual contour levels
   res0@cnMinLevelValF       = -3                 ;-- set min contour level
   res0@cnMaxLevelValF       =  3                 ;-- set max contour level
   res0@cnLevelSpacingF      =  0.3               ;-- set contour spacing
;  res0@pmLabelBarOrthogonalPosF =  0.0           ;-- position label bar

   plot(0) = gsn_csm_contour(wks,temp2,res0)

;-------------------------------------------------------------------------------
;-- middle plot - s (salinity)
;-------------------------------------------------------------------------------
   res1                      =  res               ;-- use res
;  res1@cnFillPalette        = "WhiteBlue"        ;-- set colormap
;  res1@cnFillColors         = (/5,10,20,30,40,50,60,70,80,90,100,120/)
   res1@cnLevelSelectionMode = "ManualLevels"     ;-- set manual contour levels
   res1@cnMinLevelValF       = -0.2               ;-- set min contour level
   res1@cnMaxLevelValF       =  0.2               ;-- set max contour level
   res1@cnLevelSpacingF      =  0.02              ;-- set contour spacing
   res1@pmLabelBarOrthogonalPosF =  0.038         ;-- position label bar

   plot(1) = gsn_csm_contour(wks,s2,res1)

;-------------------------------------------------------------------------------
;-- lower plot - rho (potential density)
;-------------------------------------------------------------------------------
   delete(res@tmXBLabelsOn)                       ;-- setting the resource to True won't work
   delete(res@tmXBOn)                             ;-- setting the resource to True won't work

   res2                      =  res               ;-- use res
;  res2@cnFillPalette        = "rainbow"          ;-- set colormap
;  res2@cnFillColors         = (/140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245/)
   res2@cnLevelSelectionMode = "ManualLevels"     ;-- set manual contour levels
   res2@cnMinLevelValF       = -0.6               ;-- set min contour level
   res2@cnMaxLevelValF       =  0.6               ;-- set max contour level
   res2@cnLevelSpacingF      =  0.10              ;-- set contour spacing
   res2@pmLabelBarOrthogonalPosF =  0.0           ;-- position label bar
   res2@tmXBMode             = "Explicit"         ;-- set x-axis labeling to explicit
   res2@tmXBValues           =  time2             ;-- values for x-axis tickmarks
   res2@tmXBLabels           =  ""+year           ;-- set labels equal to values (type string)
   res2@tmXBLabelStride      =  10                ;-- draw every 5th label
   res2@tmXBLabelFontHeightF =  0.012             ;-- x-axis font size
   res2@tmXBLabelAngleF      =  45                ;-- rotate the x-axis labels counter clockwise
   res2@tmXBLabelDeltaF      =  0.5               ;-- move the x-axis labels downward
   res2@tmXBLabelsOn         =  True              ;-- draw the x-axis tickmark labels
   res2@tiXAxisOn            =  True              ;-- draw the x-axis title    
   res2@tiXAxisString        = "Time"             ;-- draw x-axis title string
;  res2@tiXAxisOffsetYF      =  0.00              ;-- move x-axis title string upward
   res2@tiXAxisFontHeightF   =  0.016             ;-- x-axis font size
   res2@tmXBMinorLengthF     =  0.003             ;-- no minor ticks since stride was used above
;  res2@tmXBMajorLengthF     =  0.01
   res2@pmLabelBarOrthogonalPosF = 0.018          ;-- position label bar

   plot(2) = gsn_csm_contour(wks,rho2,res2)

;-------------------------------------------------------------------------------
;-- generate page panel
;-------------------------------------------------------------------------------
  pres                       =  True              ;-- resource object for panel
  pres@gsnPaperOrientation   = "Portrait"         ;-- set paper orientation
  pres@gsnMaximize           = True               ;-- maximize plots (don't set gsnMaximize for the other res's)
; pres@gsnPanelYF            =  (/0.95,.65,.35/)  ;-- adjust middle and lower plot
; pres@gsnPanelYF            =  (/0.92,.61,.30/)  ;-- adjust middle and lower plot
  pres@gsnPanelYF            =  (/0.94,.64,.34/)  ;-- adjust middle and lower plot
; pres@gsnPanelYF            =  (/0.945,.65,.355/)  ;-- adjust middle and lower plot
  pres@gsnPanelTop           =  0.94              ;-- set panel top
  pres@gsnPanelBottom        =  0.05              ;-- set panel bottom
  pres@txFontHeightF         =  0.019             ;-- set text font size
; pres@txString              = "ICON: Global Mean Evolution" ;-- draw title string
                                                  ;-- draw title string including file string
  pres@txString              = "ICON: Global Mean Evolution~C~~Z75~            $outputIdent"
; pres@txPosYF               =  0.98

; tres                       =  True              ;-- resource object for panel
; tres@txFontHeightF         = 0.015              ;
; gsn_text_ndc(wks,"$output",.2,.9,tres)
; gsn_text_ndc(wks,"filename hier",.4,.97,tres)   ;  changes viewport!

  gsn_panel(wks,plot,(/3,1/),pres)

end
EOF


  ncl scr_plot_tsrho.ncl
# rm scr_plot_tsrho.ncl


fi


