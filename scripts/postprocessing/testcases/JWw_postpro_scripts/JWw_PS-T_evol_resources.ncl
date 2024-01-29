; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;
;-------------------------
; Hui Wan (MPI-M, 2006.06, original version)
; Hui Wan (MPI-M, 2008.07, revised)
; Hui Wan (MPI-M, 2009.04, clean-up)
;-------------------------------------------------------------------
; These are the common plot settings for drawing snapshots of 
; surface pressure and temperature for the Jablonowski-Williamson
; baroclinic wave test.
;-------------------------------------------------------------------

     ResC@cnFillOn       = True
     ResC@cnLinesOn      = True
     ResC@cnLineLabelsOn = False

     FontHeight = 0.019

     ResC@tiXAxisFontHeightF        = FontHeight
     ResC@tiYAxisFontHeightF        = FontHeight
     ResC@tmXBLabelFontHeightF      = FontHeight
     ResC@tmYLLabelFontHeightF      = FontHeight
     ResC@gsnStringFontHeightF      = FontHeight +0.004

     ResC@tmYLLabelDeltaF = -0.6 

     ResC@mpProjection       = "CylindricalEquidistant"

     ResC@mpCenterLonF       = 180.
     ResC@mpMinLatF          =   0.
     ResC@mpMaxLatF          =  90.

     ResC@gsnMajorLatSpacing = 30.
     ResC@gsnMajorLonSpacing = 60.
     ResC@gsnMinorLatSpacing = 10.
     ResC@gsnMinorLonSpacing = 10.

     ResC@mpGeophysicalLineColor = "transparent"
     ResC@mpFillOn               = False

     ResC@pmLabelBarHeightF        = 0.06
     ResC@pmLabelBarWidthF         = 0.75
     ResC@lbLabelAutoStride        = True
     ResC@pmLabelBarOrthogonalPosF = 0.3
     ResC@lbLabelFontHeightF       = FontHeight

     ResC@cnLevelSelectionMode  = "ManualLevels" 

   ;----------------------------------------------------------------
   ; for surface pressure
   ;----------------------------------------------------------------

     psMin1 =  993.
     psMax1 = 1005.
     psInt1 =    1.

     psMin2 =  940.
     psMax2 = 1020.
     psInt2 =   10.

   ;----------------------------------------------------------------
   ; for temperature
   ;----------------------------------------------------------------

     tempMin =  230.
     tempMax =  300.
     tempInt =   10.
