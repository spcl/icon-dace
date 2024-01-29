; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  ; Auf Nachfrage: hier muss ich die CDO benutzen, 
  ; um erst mal ein lat-lon Gitter zu haben

  dir ="/home/agassman/ICON/icon-1.0.6_RESTRUCT/experiments/"
  exd ="ags0002gwts"
  gri =dir+"slice_600_grid"
  plt ="/slc_theta_prime_wind00"
  pot = dir+exd+plt

  wks = gsn_open_wks("ps",pot)

  ; get the filenames of interest
  file_nam = systemfunc("ls -1d "+dir+exd+"/"+exd+"_*")
  no_files = dimsizes(file_nam)

  ResC          = True

  ;Contour and shading
  ResC@cnFillOn          = True
  ResC@cnLinesOn         = True
  ResC@cnLevelSelectionMode= "ManualLevels"
  ResC@cnLevelSpacingF   = 0.5e-3
  ResC@cnMinLevelValF    = -4e-3
  ResC@cnMaxLevelValF    = 4e-3 


  ;Color(label) bar
  ResC@lbLabelAutoStride = True
  ResC@lbOrientation     = "Vertical"

  ;Achsenbeschriftung
  ResC@tiYAxisString     = "h [m]"
  ResC@tiXAxisString     = "x [km]"
  ResC@tiMainOn          = True
  ResC@tiMainPosition    = "Right"
  ResC@tiMainSide        = "Bottom"
  ResC@tiMainOffsetYF    = 0.05
  ResC@tiMainFontHeightF = 0.02

  ;Viewport
  ResC@vpHeightF         = 0.3

  ;gsn
  gsn_define_colormap(wks,"posneg_2")
  ResC@gsnSpreadColors     = True
  ResC@gsnSpreadColorStart = 3
  ResC@gsnSpreadColorEnd   = 21
  ResC@gsnContourZeroLineThicknessF = 0.0

  do i = 0,no_files

    if (i.eq.0) then
      system("cdo remapdis,"+gri+" -seltimestep,1 "+file_nam(i)+" tmp.nc")
    end if
    if (i.eq.1) then
      system("cdo remapdis,"+gri+" -seltimestep,2 "+file_nam(i-1)+" tmp.nc")
    end if
    if (i.gt.1) then
      system("cdo remapdis,"+gri+" -seltimestep,1 "+file_nam(i-1)+" tmp.nc")
    end if

    fin =addfile("tmp.nc","r")
 
    ; plot perturbation of potential temperature
    var      = fin->THETA_V(0,:,0,:)
    var(0,:)      = var(0,:)-331.361449957681
    var(1,:)      = var(1,:)-329.676192399636
    var(2,:)      = var(2,:)-327.999505823633
    var(3,:)      = var(3,:)-326.331346638867
    var(4,:)      = var(4,:)-324.671671476230
    var(5,:)      = var(5,:)-323.020437187183
    var(6,:)      = var(6,:)-321.377600842634
    var(7,:)      = var(7,:)-319.743119731824
    var(8,:)      = var(8,:)-318.116951361213
    var(9,:)      = var(9,:)-316.499053453378
    var(10,:)     = var(10,:)-314.889383945912
    var(11,:)     = var(11,:)-313.287900990334
    var(12,:)     = var(12,:)-311.694562950995
    var(13,:)     = var(13,:)-310.109328404001
    var(14,:)     = var(14,:)-308.532156136135
    var(15,:)     = var(15,:)-306.963005143783
    var(16,:)     = var(16,:)-305.401834631870
    var(17,:)     = var(17,:)-303.848604012798
    var(18,:)     = var(18,:)-302.303272905394
    var(19,:)     = var(19,:)-300.765801133857
    var@long_name = var@long_name+" perturbation" 
    var&lon@units = "km"

    x             = fin->lon*5.0/6.0
    y             = fin->mlev*(-500)+10000+250

    ;Coordinates of data points
    ResC@sfXArray          = x
    ResC@sfYArray          = y

    ResC@tiMainString = "U_mean=20 m/s:   "+i*5+" min"

    plot = gsn_csm_contour(wks,var,ResC)

  end do

  ; convert to png
  ;delete(wks)
  ;options=" -geometry 1000x1000 -density 300 -trim "
  ;system("convert"+options+pot+".ps "+pot+".png")
  ;system("display "+pot+".png")

end
