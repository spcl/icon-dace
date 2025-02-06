; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;*************************************************
; namelist_postpro.ncl
;
; by: Kristina Froehlich (DWD)
;
; Concepts illustrated:
;   - plots ICON model variables
;   - plots user specified list of namelist parameters
;
;*************************************************
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; load "lm_user.ncl"
;---------------------------------------------------------------

 begin

;get environment variables
  Model      = getenv("MODEL")
  testcase   = getenv("TESTCASE")
  expdir     = getenv("EXP")	
  ConfigStr  = getenv("CONFIG_STRING")
  Diri       = getenv("DIRI")
  FileName   = getenv("FNAM")
  Diro       = getenv("DIRO")
  GraphicFile= getenv("PNAME")
  pfmt       = getenv("PFMT")

  level      = stringtointeger(getenv("LEVEL"))
  timestep   = stringtointeger(getenv("TIMESTEP"))
  
;open files
  File       = addfile(FileName,"r" )                      ;open data file
  attName    = asciiread(Diro+"list.txt",-1,"string")      ;read in ascii file containing namelist parameters

  nattName = dimsizes(attName)

  TitleString = " " 
  TitleOffset = 0.000                                  

  do i = 0,nattName-1

     TitleString = TitleString + attName(i)+" = "+File@$attName(i)$ +", "

     if ((i+1)%3.eq.0) then
       TitleString = TitleString+"~C~"
       TitleOffset = TitleOffset + 0.01
     end if

  end do

;print(TitleString)
;---------------------------------------------------------------
; read in the meteorological field and grid information
;--------------------------------------------------------------

  VarNames = getfilevarnames( File )
 
  rad2deg = 45./atan(1.)              ; radians to degrees 
  
 ; coordinates of global coarse grid
  x = File->clon *rad2deg             ; cell center, lon
  y = File->clat *rad2deg             ; cell center, lat

;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; General Labeling
  RightString   = ConfigStr
  FontHeight   = 0.013
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   wks = gsn_open_wks(pfmt,Diro+GraphicFile)

;---------------------------------------------------------------
; make the plots
;---------------------------------------------------------------

; Set up resources for contour/map plot
  ResC                = True

  ResC@mpProjection         = "CylindricalEquidistant"
  ResC@mpCenterLonF         =  180.         ; choose center lon
  ResC@mpCenterLatF         =  20.          ; choose center lat
  ResC@mpLimitMode          = "LatLon"      ; required 
  ResC@mpFillOn             = False
  ResC@mpMinLatF            = -20.
  ResC@mpMaxLatF            = 90.
;  ResC@mpGreatCircleLinesOn = True
  
  ResC@gsnFrame             = False
  ResC@gsnMaximize          = False
  ResC@gsnSpreadColors      = True
  ResC@gsnStringFontHeightF = FontHeight
  ResC@gsnRightString       = RightString 
  ResC@cnFillOn             = True
  ResC@cnLinesOn            = False
  ResC@cnInfoLabelOn        = False
  ResC@cnFillMode           = "Areafill"
  ResC@cnRasterSmoothingOn  = True

  ResC@sfXArray             = x             ; These are 1D arrays, so a triangular
  ResC@sfYArray             = y             ; mesh will be created internally.

  ResC@tiMainFontHeightF        = FontHeight
  ResC@tiMainOffsetYF           = TitleOffset
  ResC@tiXAxisFontHeightF       = FontHeight
  ResC@tiXAxisOffsetYF          = FontHeight
  ResC@tiMainFont               = 21
   
  ResC@pmLabelBarHeightF        = 0.07
  ResC@pmLabelBarWidthF         = 0.7
  ResC@pmLabelBarOrthogonalPosF = 0.1
  ResC@lbLabelBarOn             = True 
  ResC@lbLabelAutoStride        = True
  ResC@lbTitlePosition          = "Bottom"
  ResC@lbTitleFontHeightF       = FontHeight	
  

;---------------------------------------------------------------
; labels and variable dependent settings
;---------------------------------------------------------------

do i = 0,dimsizes(VarNames) -1

  ; print("Varname="+VarNames(i) )

   if (VarNames(i).eq."T") then

     var = File->$VarNames(i)$(timestep,level,:)

     colormap="BlAqGrYeOrRe"

     gsn_define_colormap(wks,colormap)
    
     ;ResC@cnMinLevelValF  = (/ 220./)
     ;ResC@cnMaxLevelValF  = (/ 330./)
     ;ResC@cnLevelSpacingF = (/  3./)

     ResC@tiMainString         = TitleString
     ResC@lbTitleString        = var@long_name+" in "+var@units
     ResC@tiXAxisString        = "Minimum "+min(var)+" "+ "Maximum "+max(var)
     ResC@gsnLeftString        = var@long_name+" (level "+level+", step "+timestep+")"

     plot = gsn_csm_contour_map(wks,var,ResC)
     frame(wks)

  end if

  if ( VarNames(i) .eq."Q1") then
    var = File->$VarNames(i)$(timestep,level,:)   ; dims: (time,lev,cell)
    scale = 1.e3
    var = var*scale
    var@long_name="specific humidity"
    var@units="g/kg"

    colormap = "BlueGreen14"
    colorStart = 2
    colorEnd   = 15

    gsn_define_colormap(wks,colormap)
    gsn_reverse_colormap(wks)

     ;ResC@cnMinLevelValF  = (/ 2./)
     ;ResC@cnMaxLevelValF  = (/ 22./)
     ;ResC@cnLevelSpacingF = (/  2./)

     ResC@gsnLeftString        = var@long_name+" (level "+level+", step "+timestep+")"
     ResC@tiMainString   = TitleString
     ResC@tiXAxisString  = "Mean "+avg(var)+"   Maximum "+max(var)
     ResC@lbTitleString   = var@long_name+" in "+var@units

     plot = gsn_csm_contour_map(wks,var,ResC)
     frame(wks)
  
    end if

    if ( VarNames(i) .eq."Q2") then
      scale = 1.e3
      var = File->$VarNames(i)$(timestep,level,:)   ; dims: (time,lev,cell)
      var = var*scale
      var@long_name="specific cloud water"
      var@units="g/kg"
      
      colormap = "prcp_2"
      colorStart = 3
      colorEnd   = 15

     gsn_define_colormap(wks,colormap)

     ;ResC@cnMinLevelValF  = (/ 0./)
     ;ResC@cnMaxLevelValF  = (/ 2./)
     ;ResC@cnLevelSpacingF = (/ 0.15/)

     ResC@gsnLeftString        = var@long_name+" (level "+level+", step "+timestep+")"
     ResC@tiMainString   = TitleString
     ResC@tiXAxisString  = "Mean "+avg(var)+"  Maximum "+max(var)
     ResC@lbTitleString  = var@long_name+" in "+var@units
      
     plot = gsn_csm_contour_map(wks,var,ResC)
     frame(wks)
    
  

    end if

      if (VarNames(i) .eq."Q3") then
        var = File->$VarNames(i)$(timestep,icelevel,:)   ; dims: (time,lev,cell)
        scale = 1.e3
        var = var*scale
        var@long_name="specific cloud ice"
        var@units="g/kg"

        colormap = "prcp_2"
        colorStart = 2
        colorEnd   = 15

        gsn_define_colormap(wks,colormap)
        
        ;ResC@cnMinLevelValF  = (/ 0./)
        ;ResC@cnMaxLevelValF  = (/ 0.05/)
        ;ResC@cnLevelSpacingF = (/ 0.0025/)

        ResC@gsnLeftString        = var@long_name+" (level "+level+", step "+timestep+")"
        ResC@tiMainString   = TitleString
        ResC@tiXAxisString  = "Mean "+avg(var)+"  Maximum "+max(var)
        ResC@lbTitleString  = var@long_name+" in "+var@units

        plot = gsn_csm_contour_map(wks,var,ResC)
        frame(wks)

   end if

      if (VarNames(i) .eq."Q7") then ;if (VarNames(i) .eq."Q4") then
        
        var = File->$VarNames(i)$(timestep,level,:)   ; dims: (time,lev,cell)
        scale = 1.e3
        var = var*scale
        var@long_name="specific rain water"
        var@units="g/kg"      
 
        colormap = "prcp_2"
        colorStart = 2
        colorEnd   = 15

        gsn_define_colormap(wks,colormap)

        ResC@cnMinLevelValF  = (/ 0./)
        ResC@cnMaxLevelValF  = (/ 0.05/)
        ResC@cnLevelSpacingF = (/ 0.02/)

        ResC@gsnLeftString        = var@long_name+" (level "+level+", step "+timestep+")"
        ResC@tiMainString   = TitleString
        ResC@tiXAxisString  = "Mean "+avg(var)+"  Maximum "+max(var)
        ResC@lbTitleString  = var@long_name+" in "+var@units

        plot = gsn_csm_contour_map(wks,var,ResC)
        frame(wks)

         print("max QR="+max(var))
         print("min QR="+min(var))
     end if

     if (VarNames(i) .eq."Q5") then
        scale = 1.e3
      var = File->$VarNames(i)$(timestep,icelevel,:)   ; dims: (time,lev,cell)
      var = var*scale
      var@long_name="specific snow"
      var@units="g/kg"

      colormap = "prcp_2"
      colorStart = 2
      colorEnd   = 15

      gsn_define_colormap(wks,colormap)

      ;ResC@cnMinLevelValF  = (/ 0./)
      ;ResC@cnMaxLevelValF  = (/ 0.5/)
      ;ResC@cnLevelSpacingF = (/ 0.05/)

      ResC@gsnLeftString        = var@long_name+" (level "+level+", step "+timestep+")"
      ResC@tiMainString   = TitleString
      ResC@tiXAxisString  = "Mean "+avg(var)+"  Maximum "+max(var)
      ResC@lbTitleString  = var@long_name+" in "+var@units

      plot = gsn_csm_contour_map(wks,var,ResC)
      frame(wks)

    end if

  if (VarNames(i) .eq."PS") then

        var = File->$VarNames(i)$(timestep,:)   ; dims: (time,cell)
        scale = 1.e-2    
        var = var*scale
        var@units="hPa"

        colormap = "testcmap"
        colorStart = 2
        colorEnd   = 200
        
        ;varMin     = (/ 940./)
        ;varMax     = (/ 1040./)
        ;varInt     = (/  5./)
     
        gsn_define_colormap(wks,colormap)

        ;ResC@cnMinLevelValF  = varMin
        ;ResC@cnMaxLevelValF  = varMax 
        ;ResC@cnLevelSpacingF = varInt
        
        ResC@gsnLeftString        = var@long_name+" (level "+level+", step "+timestep+")"
        ResC@tiMainString   = TitleString

        ResC@lbTitleString  = var@long_name+" in "+var@units
        ResC@tiXAxisString  = "Minimum "+min(var)+" "+ "Maximum "+max(var)

        plot = gsn_csm_contour_map(wks,var,ResC)
        frame(wks)

    end if

  if (VarNames(i) .eq."DIV") then
  
        var = File->$VarNames(i)$(timestep,level,:)   ; dims: (time,lev,cell)
        scale = 1.e3    
        var = var*scale*3.6
        var@units="1/h"

        colormap = "testcmap"
        colorStart = 2
        colorEnd   = 200
               
        gsn_define_colormap(wks,colormap)

        ;ResC@cnMinLevelValF  = (/ -0.2/)
        ;ResC@cnMaxLevelValF  = (/ 0.2/)
        ;ResC@cnLevelSpacingF = (/ 0.005/)

        ResC@gsnLeftString        = var@long_name+" (level "+level+", step "+timestep+")"
        ResC@tiMainString   = TitleString

        ResC@tiMainString   = TitleString
        ResC@lbTitleString  = var@long_name+" in "+var@units
        ResC@tiXAxisString  = "Minimum "+min(var)+" "+ "Maximum "+max(var)
        ResS = ResC
        symMinMaxPlt (var,20,False,ResS)
        plot = gsn_csm_contour_map(wks,var,ResS)
        frame(wks)

      end if

      if (VarNames(i) .eq."OMEGA") then

        var = File->$VarNames(i)$(timestep,level,:)   ; dims: (time,lev,cell)
        scale =1.
        var = var*scale
        var@units="Pa/s"

        colormap = "BlWhRe"
        colorStart = 2
        colorEnd   = 102

        gsn_define_colormap(wks,colormap)

        ;ResC@cnMinLevelValF  = (/ -2/)
        ;ResC@cnMaxLevelValF  = (/ 2/)
        ;ResC@cnLevelSpacingF = (/ 0.02/)t

        ResC@gsnLeftString  = var@long_name+" (level "+level+", step "+timestep+")"
        ResC@tiMainString   = TitleString
        ResC@lbTitleString  = var@long_name+" in "+var@units
        ResC@tiXAxisString  = "Minimum "+min(var)+" "+ "Maximum "+max(var)

        ResS = ResC
        symMinMaxPlt (var,20,False,ResS)
        plot = gsn_csm_contour_map(wks,var,ResS)
        frame(wks)

      end if

      if (VarNames(i).eq."V") then

        var = File->$VarNames(i)$(timestep,level,:)   ; dims: (time,lev,cell)
        scale = 1.
        var = var*scale
        var@units="m/s"

        colormap = "testcmap"
        colorStart = 2
        colorEnd   = 200
     
        gsn_define_colormap(wks,colormap)

        ;ResC@cnMinLevelValF  = (/ -40/)
        ;ResC@cnMaxLevelValF  = (/ 40/)
        ;ResC@cnLevelSpacingF = (/ 3/)

        ResC@gsnLeftString        = var@long_name+" (level "+level+", step "+timestep+")"
        ResC@tiMainString   = TitleString
        ResC@lbTitleString  = var@long_name+" in "+var@units
        ResC@tiXAxisString  = "Minimum "+min(var)+" "+ "Maximum "+max(var)
        
        ResS = ResC
        symMinMaxPlt (var,20,False,ResS)
        plot = gsn_csm_contour_map(wks,var,ResS)
        frame(wks)

      end if

      if (VarNames(i).eq."U") then

        var = File->$VarNames(i)$(timestep,level,:)   ; dims: (time,lev,cell)
        scale = 1.
        var = var*scale
        var@units="m/s"

        colormap = "testcmap"
        colorStart = 2
        colorEnd   = 200
     
        gsn_define_colormap(wks,colormap)

        ;ResC@cnMinLevelValF  = (/ -40/)
        ;ResC@cnMaxLevelValF  = (/ 40/)
        ;ResC@cnLevelSpacingF = (/ 3/)

        ResC@gsnLeftString        = var@long_name+" (level "+level+", step "+timestep+")"
        ResC@tiMainString   = TitleString
        ResC@tiMainString   = TitleString
        ResC@lbTitleString  = var@long_name+" in "+var@units
        ResC@tiXAxisString  = "Minimum "+min(var)+" "+ "Maximum "+max(var)
        
        ResS = ResC
        symMinMaxPlt (var,20,False,ResS)
        plot = gsn_csm_contour_map(wks,var,ResS)
        frame(wks)

      end if

      if (VarNames(i).eq."W") then

        var = File->$VarNames(i)$(timestep,level,:)   ; dims: (time,lev,cell)
        scale = 1.
        var = var
        var@units="m/s"

        colormap = "testcmap"
        colorStart = 2
        colorEnd   = 200
        
        ;varMin     = (/ -40/)
        ;varMax     = (/ 40/)
        ;varInt     = (/ 3/)
     
        gsn_define_colormap(wks,colormap)
        ;ResC@cnMinLevelValF  = (/ -40/)
        ;ResC@cnMaxLevelValF  = (/ 40/)
        ;ResC@cnLevelSpacingF = (/ 3/)

        ResC@gsnLeftString        = var@long_name+" (level "+level+", step "+timestep+")"
        ResC@tiMainString   = TitleString
        ResC@lbTitleString  = var@long_name+" in "+var@units
        ResC@tiXAxisString  = "Minimum "+min(var)+" "+ "Maximum "+max(var)
        
        ResS = ResC
        symMinMaxPlt (var,20,False,ResS)
        plot = gsn_csm_contour_map(wks,var,ResS)
        frame(wks)

      end if

  ;if (VarNames(i) .eq."APRL") then
  
  ;      var = File->$VarNames(i)$(timestep,:)   ; dims: (time,cell)
  ;      scale = 8.64e4    
  ;      var = var*scale
  ;      var@units="mm/d"


  ;      colormap   = "WhViBlGrYeOrReWh"
  ;      colorStart = 2
  ;      colorEnd   = 25
               
  ;      gsn_define_colormap(wks,colormap)

  ;      ResC@gsnLeftString        = "APRL (level "+level+", step "+timestep+")"
  ;      ResC@tiMainString   = TitleString

  ;      ResC@tiMainString   = TitleString
  ;      ResC@lbTitleString  = "APRL in mm/d"
  ;      ResC@tiXAxisString  = "Minimum "+min(var)+" "+ "Maximum "+max(var)
  ;      ResS = ResC
  ;      symMinMaxPlt (var,10,False,ResS)
  ;      plot = gsn_csm_contour_map(wks,var,ResS)
  ;      frame(wks)

  ; end if

end do

end


