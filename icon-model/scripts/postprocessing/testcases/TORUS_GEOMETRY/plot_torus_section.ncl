; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;------------------------------------
; NCL script for plotting a single scalar variable. 
; 
; WARNING: this script should be called through the driver
; "plot_scalar_and_vector", it should not be launched directly from
; the user. See plot_scalar_and_vector -h for details.
;

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

 begin

;---------------------------------------------------------------
; Get environment variables
;---------------------------------------------------------------
   
  varname     = getenv("varname")

  input_file  = getenv("input_file")

  ; time levels
  istep_start = stringtointeger( getenv("time_level_start") )
  istep_end   = stringtointeger( getenv("time_level_end") )
  istep_step  = stringtointeger( getenv("time_level_step") )

  ; fixed level (lat/lon/vertical)
  fix_dim     = getenv("fix_dimension")
  ilev        = stringtointeger( getenv("fix_level") )
 
  colormap    = getenv("colormap")

  ncontours   = stringtointeger( getenv("contour_levels_ncont") )
  if (ncontours.gt.0) then
    contour_start = stringtodouble( getenv("contour_levels_start") )
    contour_end   = stringtodouble( getenv("contour_levels_end") )
  end if

  output_terminal = getenv("output_terminal")
  output_file = getenv("output_file")

;---------------------------------------------------------------
; Open data file and get grid information and fields
;---------------------------------------------------------------

  File = addfile( input_file, "r" )
  
  ; read the coordinates
  x         = File->lon            ; longitude [deg]
  y         = File->lat            ; latitude  [deg]

  if(varname.eq."W")then
    z       = File->lev_2
  else   
    z       = File->lev
  end if

  times       = File->time            ; time      [days]

  ; negative start or end indexes indicate the largest value
  if (istep_start.lt.0) then
    istep_start = dimsizes(times)-1
  end if
  if (istep_end.lt.0) then
    istep_end = dimsizes(times)-1
  end if
  nsteps = (istep_end-istep_start)/istep_step + 1

  ; read the variable
  isteps = new(nsteps,integer)
  do i = 0, nsteps-1
    isteps(i) = istep_start + i*istep_step
  end do

  ;get variable to plot
   var = File->$varname$

  varlongname = var@long_name
  varunit     = var@units

;---------------------------------------------------------------
; Setup resources
;---------------------------------------------------------------

  res = True

  ; coordinates of data points
  if(fix_dim.eq."lev")then
    res@sfXArray = x
    res@sfYArray = y
    res@tiYAxisString    = "latitude"
  else if(fix_dim.eq."lat")then
    res@sfXArray = x
    res@sfYArray = z
    res@tiYAxisString    = "Vertical levels"
    res@trYReverse       = True
  else 
    res@sfXArray = y
    res@sfYArray = z
    res@tiYAxisString    = "Vertical levels"
    res@trYReverse       = True
  end if
  end if

  ;res@cnFillOn  = True  ; use colors
  ;res@cnLinesOn = False ; do not mark contour lines
  ;res@cnLineLabelsOn  = False    ; turns off contour line labels

  ;res@gsnSpreadColors = True ; use the whole color bar

  if (ncontours.gt.0) then
    contour_levs = new(ncontours,double)
    if (ncontours.gt.1) then
      contour_interval = (contour_end-contour_start)/(ncontours-1)
      do i = 0, ncontours-2
        contour_levs(i) = contour_start + i*contour_interval
      end do
    end if
    ; make sure the last element is exact
    contour_levs(ncontours-1) = contour_end
    res@cnLevelSelectionMode = "ExplicitLevels" 
    res@cnLevels = contour_levs
  end if

  ;Labelbar's font size
  res@lbLabelFontHeightF     = 0.016
  res@gsnStringFontHeightF   = 0.018
  res@tmXBLabelFontHeightF   = 0.016
  res@tmYLLabelFontHeightF   = 0.016
 
 ;dashed negative contours
 res@gsnContourPosLineDashPattern = 2
 res@gsnContourZeroLineThicknessF = 0 	    ; eliminates zero contour

;---------------------------------------------------------------
; Make the plot
;---------------------------------------------------------------

  ; create the workspace and 
  wks  = gsn_open_wks(output_terminal,output_file)

  ; make the plot
  gsn_define_colormap(wks,colormap)

  ;adjust axis length
  res@trXMinF = -40
  res@trXMaxF = 40
  res@trYMinF = 30

  if(fix_dim.eq."lev")then
    print("Plotting a horizontal section:")
    do i = 0, nsteps-1
      ;res@gsnCenterString = ;+" step "+isteps(i)
      plot = gsn_csm_contour(wks,var(isteps(i),ilev-1,:,:),res)
    end do
  else if(fix_dim.eq."lat")then
    print("Plotting a vertical section at fix lat:")
    do i = 0, nsteps-1
      ;res@gsnCenterString = ;+" step "+isteps(i)
      plot = gsn_csm_contour(wks,var(isteps(i),:,ilev-1,:),res)
    end do
  else 
    print("Plotting a vertical section at fix lat:")
    do i = 0, nsteps-1
      ;res@gsnCenterString = ;+" step "+isteps(i)
      plot = gsn_csm_contour(wks,var(isteps(i),:,:,ilev-1),res)
    end do
  end if
  end if

end
