; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;-----------------------------------------------------------------------------
; Hui Wan (MPI-M, 2011-05, Hamburg)
;-----------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;-----------------------------------------------------------------------------

begin

 PlotFormat = getenv("plot_file_format")
 PlotFile   = getenv("PlotFile")
 wkOrientation = getenv("wkOrientation")

 varName    = getenv("varname")
 LongName   = getenv("LongName")
 Scale      = stringtofloat(getenv("Scale"))

 DataPath    = getenv("DataPath")
 DataFNameP1 = getenv("DataFNameP1")
 DataFNameP2 = getenv("DataFNameP2")
 ifs         = stringtoint(getenv("evol_istart"))
 ife         = stringtoint(getenv("evol_iend"))

;--------------------------------------------------------------- 
; Open data files, inquire dimension sizes and starting time
;--------------------------------------------------------------- 

  cmd      = "ls "+DataFNameP1+"????"+DataFNameP2
  filelist = systemfunc(cmd)  
  nfile    = dimsizes(filelist)

  if (nfile.ne.(ife-ifs+1)) then
    print("!!!!!!!!!!!!!! fatal problem with data files !!!!!!!!!!!!!!!")
    print(" User specified starting and ending file indices: "+ifs+", "+ife)
    print(" But "+nfile+" files were found. Command was")
    print("   >ls "+DataFNameP1+"????"+DataFNameP2)
    print("!!!!!!!!!!!!!!!!!!!!!!!!! aborted !!!!!!!!!!!!!!!!!!!!!!!!!!!")
    exit
  end if

  Files = addfiles(filelist,"r")

  lat   = Files[0]->lat
  nlat  = dimsizes(lat) 

  ntime = dimsizes(Files[:]->time)
  print("Found "+ntime+" time steps in "+nfile+" files")

  tmpunit  = Files[0]->time@units
  refdate  = str_get_field(tmpunit,3," ")+" "+str_get_field(tmpunit,4," ")
  print("Reference date of the time axis is "+refdate)

  timeunit = "days since "+refdate
  delete(tmpunit)

;--------------------------------------------------------------- 
; Read data
;--------------------------------------------------------------- 

  var  = new( (/ntime,nlat/),typeof(Files[0]->$varName$) )
  time = new( ntime,typeof(Files[0]->time) )

  its = 0
  do ifile = 0,nfile-1

     ttmp = Files[ifile]->time
     if (ttmp@calendar.eq."proleptic_gregorian") then
        ttmp@calendar = "gregorian"
     end if

     nts = dimsizes(ttmp)
     ite = its + nts- 1

     ;Read variable
     var(its:ite,:) = (/Files[ifile]->$varName$(:,0,:,0)/) ;(time,lev,lat,lon)

     ;Get relative time 
     time(its:ite) = ut_convert(ttmp,timeunit)

     print("Read file "+sprinti("%0.4i", ifile+ifs) \
          +" (containing "+nts+" steps, "+" global time indices = "+its+":"+ite+")")

     its = ite + 1
     delete(ttmp)
  end do 

  time@units = "days since "+refdate

  var      = var*Scale
  var!0    =  "time"
  var&time =   time
  var!1    =  "lat"
  var&lat  =   lat  

;--------------------------------------------------------------- 
; Make plot
;--------------------------------------------------------------- 

  FontHeight = 0.014

  ResC = True

  ; Contour plot with fixed lat:time aspect ratio 
  ResC@vpHeightF   = 0.3 
  ResC@vpWidthF    = ResC@vpHeightF*1.618*time(ntime-1)/200. 
  if (ResC@vpWidthF.gt.0.8) then
    ResC@vpHeightF = ResC@vpHeightF*0.8/ResC@vpWidthF
    ResC@vpWidthF  = 0.8
  end if
 
  ResC@tiYAxisString        = "Latitude"
  ResC@tiYAxisFontHeightF   = FontHeight
  ResC@tmYLLabelFontHeightF = FontHeight 
  ResC@tmYLMajorLengthF     = 0.01
  ResC@tmYLMinorLengthF     = 0.005
  ResC@tmYRMajorLengthF     = 0.01
  ResC@tmYRMinorLengthF     = 0.005
  ResC@tmYLMajorThicknessF  = 3. 
  ResC@tmYRMajorThicknessF  = 3.
 
  ResC@tmXBLabelFontHeightF = FontHeight 
  ResC@tmXBMajorLengthF     = 0.01
  ResC@tmXBMinorLengthF     = 0.005
  ResC@tmXBMajorThicknessF  = 3. 
  ResC@tmXTMajorThicknessF  = 3. 
 
  ResC@cnFillOn      = True
  ResC@cnLinesOn     = False
  ResC@cnInfoLabelOn = False
 
  ResC@lbLabelBarOn        = True
  ResC@lbOrientation       = "vertical"
  ResC@pmLabelBarWidthF    = 0.07
  ResC@lbLabelFontHeightF  = FontHeight 
  ResC@lbLabelAutoStride   = True
  ResC@gsnSpreadColors     = True
  ResC@gsnSpreadColorStart = 3
  ResC@gsnSpreadColorEnd   = 101

  ResC@tiMainString        = getenv("exp_config") 
  ResC@tiMainFontHeightF   = FontHeight + 0.01
  ResC@gsnLeftString       = LongName +", scale = "+Scale+", mean = "+avg(var)
  ResC@gsnCenterString     = ""
  ResC@gsnRightString      = getenv("resolution") 
  ResC@gsnStringFontHeightF= FontHeight + 0.001
  ResC@tiXAxisJust         = "CenterCenter"
  ResC@tiXAxisString       = "Time (day)"
  ResC@tiXAxisFontHeightF  = FontHeight + 0.001

  ; Set the workstation orientation to "landscape" if the user wants so 
  if ( (.not.ismissing(wkOrientation)) .and. (wkOrientation.eq."landscape") ) then
    PlotFormat@wkOrientation  = "landscape"   ;NCL default is "portrait" 
  end if
  wks = gsn_open_wks(PlotFormat,PlotFile)
  gsn_define_colormap(wks,"WhBlGrYeRe")
 
  contour = gsn_csm_lat_time( wks, add90LatY(var(lat|:,time|:)), ResC )

end
