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
;ncl script to plot a 2D field from a Netcdf file
;
;Pilar Ripodas
;,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

rad2deg = 45./atan(1.)   ; radians to degrees

;nc file

FileD=getenv("FileD")
FileN=getenv("FileN")
PlotFileName=getenv("PlotFileName")
Pvlev=getenv("Plotvertlev")
vlev=stringtoint(Pvlev)
Tstep=getenv("TimS")
tstep=stringtoint(Tstep)
DefBar=getenv("BDefBar")
MaxV=getenv("BMaxV")
MinV=getenv("BMinV")
IntV=getenv("BIntV")
MTitle=getenv("mytitle")
mount=getenv("lmount")
MountD=getenv("MountD")
VorAtCorners = stringtoint(getenv("vorticity_at_corners"))

Units=" "

FileP= FileD + FileN

 File = addfile(FileP,"r")


 f_att=getvaratts(File)
 if (isatt(File, "nroot"))then
  Root=File@nroot
 else 
  if (isatt(File, "grid_root")) then
   Root=File@grid_root
  else
   Root=" "
  end if
 end if
 if (isatt(File, "start_lev"))then
  Bis=File@start_lev
 else 
  if (isatt(File, "grid_level")) then
   Bis=File@grid_level
  else
   Bis=" "
  end if
 end if
 if (isatt(File, "nlev"))then
  Nvl=File@nlev
 else
  if (isatt(File, "run_ctl:nlev"))then
   att="run_ctl:nlev"
   Nvl=File@$att$
  else
   Nvl=" "
  end if
 end if
 if (isatt(File, "i_cell_type"))then
  CellType=File@i_cell_type
 else
  if (isatt(File, "run_ctl:i_cell_type"))then
   att="run_ctl:i_cell_type"
   CellType=File@$att$
  else
   CellType="3"
  end if
 end if
 if (isatt(File, "dtime"))then
  dtime=File@dtime
 else
  if (isatt(File, "run_ctl:dtime"))then
   att="run_ctl:dtime"
   dtime=File@$att$
  else
   dtime=" "
  end if
 end if
 if (MTitle .eq. "nondef") then
  if (isatt(File, "out_expname"))then
   expname=File@out_expname
  else
   if (isatt(File, "io_ctl:out_expname"))then  
    att="io_ctl:out_expname"
    expname=File@$att$
   else
    expname=" "
   end if
  end if
 else
    expname=MTitle
 end if

 print(Bis)
 print(Root)
 print(Nvl)

;Hresol=getenv("Hresol")
;Nvl=getenv("Nvl")
;Opt=getenv("Opt")

 VarN=getenv("VarN")

; Shape of control volume (3 = triangle, 6 = hexagon/pentagon) 

;---------------------------------------------------------------
; 
;---------------------------------------------------------------

 if (CellType.eq."3") then
      if ((VarN .eq."VOR").and.(VorAtCorners.eq.1)) then
         latName    = "vlat" 
         lonName    = "vlon"                                            
      else                                                              
         latName    = "clat"                                            
         lonName    = "clon"                                            
      end if     
      gridName  = "tri"                                                       
 else                                                                 
   if (CellType.eq."6") then                                            
      if ((VarN .eq."VOR").and.(VorAtCorners.eq.1)) then
         latName    = "clat"                                            
         lonName    = "clon"                                            
      else                                                              
         latName    = "vlat"
         lonName    = "vlon"                                            
      end if  
      gridName  = "hex"                                                      
   else                                                                
      print("Wrong cell type. Should be 3 or 6")                        
   end if                                                               
 end if                                                               
                                       

;coordinates
 x = File->$lonName$ *rad2deg
 y = File->$latName$ *rad2deg 
 x = where(x .lt. 0, x+360., x)
;
;Resources (contour)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ResC = True
ResC@gsnDraw         = True
ResC@gsnFrame        = True
ResC@mpProjection   = "CylindricalEquidistant"
ResC@mpLimitMode    = "LatLon"
ResC@mpCenterLonF       = 180.01
ResC@mpMinLonF          = 0.1
ResC@mpMaxLonF          = 359.
;ResC@mpCenterLatF       = 0.05
ResC@mpMinLatF          =  -89.
ResC@mpMaxLatF          =  89.
;;
;ResC@mpCenterLonF       = 165.
;ResC@mpMinLonF          = 120.
;ResC@mpMaxLonF          = 270.
;ResC@mpCenterLatF       = 40.
;ResC@mpMinLatF          =  25.
;ResC@mpMaxLatF          =  75.
;;

ResC@mpGeophysicalLineColor = "transparent"
ResC@mpFillOn               = False
ResC@cnFillOn       = True
ResC@lbLabelBarOn   = True
ResC@lbLabelFontHeightF= 0.012
ResC@lbLabelAutoStride = True
ResC@cnLinesOn      = False
ResC@cnLineLabelsOn  = False
ResC@gsnStringFontHeightF  =0.012
ResC@tmXBLabelFontHeightF  =0.012
ResC@tmYLLabelFontHeightF  =0.012
ResC@tiMainFontHeightF  =   0.015
ResC@gsnSpreadColors  =  True
ResC@gsnMaximize    = True
ResC@gsnPaperOrientation="portrait"


ResC@gsnDraw   =  False
ResC@gsnFrame        = False 
;plot
print(PlotFileName)
wks = gsn_open_wks("eps",PlotFileName)


;gsn_define_colormap(wks,"BlGrYeOrReVi200")
gsn_define_colormap(wks,"gui_default")
;gsn_define_colormap(wks,"JWw_PS_colormap")

 ;printVarSummary(x)

 x@units="degrees_east"
 y@units="degrees_north"
 ResC@sfXArray = x  ;triangular mesh
 ResC@sfYArray = y
; x = where(x .lt.0, x + 360., x)  

 ResC@tiMainString  = expname

 
 if (DefBar .eq. "yes") then
  ResC@cnLevelSelectionMode = "ManualLevels"
  ResC@cnMinLevelValF  = stringtofloat(MinV)
  ResC@cnMaxLevelValF  = stringtofloat(MaxV) 
  ResC@cnLevelSpacingF = stringtofloat(IntV) 
 end if
 ;variable to plot 
 
 day= File ->time(tstep-1)
 iday=doubletoint(day)
 sday=sprinti("%0.8i",iday) 
 
; ResC@gsnCenterString= "day " +sday
 ResC@gsnRightString = "day " +sday + " global  R"+Root+"B0" + Bis + "L"+ Nvl+" "+dtime+"s"

 if (Pvlev .eq. "-1") then
   if (Tstep .eq. "-1")
    var= File ->$VarN$(:)
   else
    var= File ->$VarN$(tstep-1,:)
   end if
   ifsub=isStrSubset(VarN,"pres_sfc")
   if (ifsub) then
    var=var/100.
    Units="hPa"
   end if
   Pvlev="sfc"   
 else

  if (Pvlev .eq. "p") then
   presl=File ->lev
   presl=presl/100.
   var= File ->$VarN$(tstep-1,0,:)
   ResC@gsnLeftString  =  VarN + " " +presl +" hPa"
  else
   print(Tstep)
   print(tstep-1)
   print(Pvlev)
   print(vlev-1)
   var= File ->$VarN$(tstep-1,vlev-1,:)
  end if
 end if
 var@_FillValue=1.e36

 printVarSummary(var)
 varmin=min(var)
 varmax=max(var) 
 svarmin=sprintf("%8.5f",varmin)
 svarmax=sprintf("%8.5f",varmax)
 
 ResC@gsnLeftString  =  VarN + " Lev "+ Pvlev + "    min. "+ svarmin+ " max. "+ svarmax+" "+Units
 
plot=gsn_csm_contour_map(wks,var,ResC)


;Mountain

if (mount .eq. "yes") then
MountF=MountD+"/MRW_mountain_R2B04.nc"
print(MountF)
Filem = addfile(MountF,"r")
varm= Filem->ORO(0,:)
varm@units="m"
varm@long_name ="orography"
;coordinates

 xm = Filem->clon *rad2deg
 ym = Filem->clat *rad2deg 
 xm = where(xm .lt. 0, xm+360., xm)
 xm@units="degrees_east"
 ym@units="degrees_north"
 

 ResCM=True
 ResCM@cnFillOn       = False
 ResCM@cnLinesOn      = True
 ResCM@cnLineLabelsOn = False
 ResCM@cnInfoLabelOn   = False ; do not plot max/min information
 ResCM@gsnRightString = " "
 ResCM@gsnLeftString = " "
 ResCM@gsnDraw   =  False
 ResCM@gsnFrame       = False 
 ResCM@sfXArray = xm  ;triangular mesh
 ResCM@sfYArray = ym
 ResCM@cnLevelSelectionMode = "ManualLevels"
 ResCM@cnMinLevelValF  = 100
 ResCM@cnMaxLevelValF  = 2000 
 ResCM@cnLevelSpacingF = 200
 printVarSummary(xm)
 printVarSummary(ym)
 printVarSummary(varm)

 
plotm=gsn_csm_contour(wks,varm,ResCM)
overlay(plot,plotm)

end if

draw(plot)
frame(wks)  
end
