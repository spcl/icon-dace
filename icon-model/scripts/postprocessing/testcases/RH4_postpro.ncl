; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;----------------------------------------------------
; Hui Wan (MPI-Met, 2008-07-05): first version
; Hui Wan (MPI-Met, 2008-07-25): modified for ECHAM
; Hui Wan (MPI-Met, 2008-10-23)
; Hui Wan (MPI-Met, 2009-05-19): rewritten for the ICON test suite
;------------------------------------
; Script type: visualization
;---------------------------------------------------------------
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;---------------------------------------------------------------
 begin

  PlotFormat = getenv("plot_file_format")
  Model      = getenv("Model")
  Data       = getenv("Data")
  PlotPath   = getenv("PlotPath")
  Resolution = getenv("Resolution")
  ConfigStr  = getenv("ConfigStr")
  ExpName    = getenv("ExpName")
  CellType   = getenv("cell_type")
  Day        = getenv("diag_day")

  rad2deg = 45./atan(1.)   ; radians to degrees

;---------------------------------------------------------------
; title of each panel
;---------------------------------------------------------------

  LeftString   = (/"surface pressure (hPa)",                 \
                   "500 hPa geop. height (gpm)",             \
                   "850 hPa zonal wind (m s~S~-1~N~)",       \
                   "850 hPa meridional wind (m s~S~-1~N~)",  \
                   "850 hPa relative vorticity (10~S~-5~N~ s~S~-1~N~)", \
                   "850 hPa divergence (10~S~-6~N~ s~S~-1~N~)",         \
                   "850 hPa temperature (K)",    \
                   "850 hPa omega (Pa s~S~-1~N~)"/)

  RightString    = LeftString
  RightString(:) = ""
  RightString(0) = ConfigStr

;---------------------------------------------------------------
; unit conversion
;---------------------------------------------------------------

     scale = (/ 0.01, 1, 1, 1, 1e5, 1e6, 1, 1/)

;---------------------------------------------------------------
; vertical level
;---------------------------------------------------------------

  pres = (/-999.,50000.,85000.,85000.,85000.,85000.,85000.,85000./)

;---------------------------------------------------------------
; variables
;---------------------------------------------------------------

  if (Model.eq."ECHAM" .or. Model.eq."ECHAM5")  then

      varName = (/"aps","geopoth","u","v","svo","sd","st","omega"/)
  else
      varName = (/"PS","ZF3","U","V","VOR","DIV","T","OMEGA"/)
  end if

;---------------------------------------------------------------
; lat./lon. of cell center (ICOHDC only)
;---------------------------------------------------------------

  if (Model.eq."ICOHDC") then
     if (CellType.eq."3") then
        latName    = "clat"
        lonName    = "clon"
     else
     if (CellType.eq."6") then
        latName    = "vlat"
        lonName    = "vlon"
     else
        print("Wrong cell type. Should be 3 or 6")
     end if
     end if
  end if

;---------------------------------------------------------------
; set the contour levels
;---------------------------------------------------------------

  Min  = (/ 960, 5200,  0, -16, -2.5,  -1, 281.2, -0.020/)
  Max  = (/1030, 5800, 24,  16,  2.5,   1, 282.0,  0.020/)
  Incr = (/   5,   50,  3,   4,  0.5, 0.2,   0.1,  0.005/)

;---------------------------------------------------------------
; colors
;---------------------------------------------------------------

  colorStart = (/   5,   5,   5, 107, 107, 107,   5, 107/)
  colorEnd   = (/ 101, 101, 101, 197, 197, 197, 101, 197/)

;---------------------------------------------------------------
; make plots
;---------------------------------------------------------------

  nvar = dimsizes(varName)
  plot = new( nvar,graphic )

  PlotFileName = PlotPath + ExpName + "_day"+Day

  wks = gsn_open_wks(PlotFormat,PlotFileName)
  gsn_merge_colormaps(wks,"WhBlGrYeRe","ViBlGrWhYeOrRe")

;---------------------------------------
; resources

  ResC                = True

  ResC@gsnDraw        = False
  ResC@gsnFrame       = False

  FontHeight = 0.02

  ResC@tiXAxisFontHeightF        = FontHeight
  ResC@tiYAxisFontHeightF        = FontHeight
  ResC@tmXBLabelFontHeightF      = FontHeight
  ResC@tmYLLabelFontHeightF      = FontHeight
  ResC@gsnLeftStringFontHeightF  = FontHeight + 0.004
  ResC@gsnCenterStringFontHeightF= FontHeight + 0.004
  ResC@gsnRightStringFontHeightF = FontHeight + 0.004

  ResC@cnFillOn             = True
  ResC@cnLinesOn            = True 
  ResC@cnLineLabelsOn       = False

  ResC@lbLabelBarOn             = True
  ResC@pmLabelBarHeightF        = 0.08
  ResC@pmLabelBarWidthF         = 0.75
  ResC@lbLabelAutoStride        = True
  ResC@pmLabelBarOrthogonalPosF = 0.16
  ResC@lbLabelFontHeightF       = FontHeight + 0.002

  ResC@gsnSpreadColors      = True

  ResC@mpProjection       = "CylindricalEquidistant"
  ResC@mpCenterLatF       = 0.
  ResC@mpCenterLonF       = 180.
  ResC@gsnMajorLonSpacing = 60.

  ResC@mpGeophysicalLineColor = "transparent"
  ResC@mpFillOn               = False

  ResC@gsnCenterString = ""

;---------------------------------------------------------------
; read data and make contour plot
;---------------------------------------------------------------

  File = addfile( Data, "r" )

  if (Model.eq."ICOHDC") then
     ResC@sfXArray = File->$lonName$ *rad2deg
     ResC@sfYArray = File->$latName$ *rad2deg
  end if

  do ivar=0,nvar-1

     if (Model.eq."ECHAM".or.Model.eq."ECHAM5") then
      if (ivar.eq.0)
         var = File->$varName(ivar)$(0,:,:)             ;(time,lat,lon)
      else
         var = File->$varName(ivar)$(0,{pres(ivar)},:,:)  ;(time,lev,lat,lon)
      end if

     else
      if (ivar.eq.0)
         var = File->$varName(ivar)$(0,:)               ;(time,cell)
      else
         var = File->$varName(ivar)$(0,{pres(ivar)},:)  ;(time,lev,cell)
      end if
    end if

    var = var*scale(ivar)

    ResC@gsnLeftString  = LeftString(ivar)
    ResC@gsnRightString = RightString(ivar)

    ResC@cnLevelSelectionMode = "ManualLevels"
    ResC@cnMinLevelValF       =  Min(ivar)
    ResC@cnMaxLevelValF       =  Max(ivar)
    ResC@cnLevelSpacingF      =  Incr(ivar)

    ResC@gsnSpreadColorStart  = colorStart(ivar)
    ResC@gsnSpreadColorEnd    = colorEnd(ivar)

    plot(ivar) = gsn_csm_contour_map(wks,var,ResC)

  delete(var)
  end do

;---------------------------------------------------------------
; put all panels on one page
;---------------------------------------------------------------

  ResP                            = True
  ResP@gsnPanelYWhiteSpacePercent = 8
  ResP@gsnPanelXWhiteSpacePercent = 5
  ResP@gsnMaximize = True

  gsn_panel (wks, plot, (/nvar/2,2/), ResP)

end
