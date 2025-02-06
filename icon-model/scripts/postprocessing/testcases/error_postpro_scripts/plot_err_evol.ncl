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
; Pilar Ripodas, DWD
; Modified by Constantin Junk, MPI-M, 2010-12-07
;------------------------------------

;Plot time evolution of l1-, l2- and l_ind- error norm

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;;;;;;;;;;;;;;;;;;;

undef("labelres")
procedure labelres(res:logical,legendl)

begin
res@pmLegendDisplayMode="Always"
res@pmLegendSide="Top"
res@pmLegendZone=1
res@pmLegendParallelPosF=0.2
res@pmLegendOrthogonalPosF=0.2
res@pmLegendWidthF=0.1
res@pmLegendHeightF=0.2
res@lgLabelFontHeightF=0.02
res@lgPerimOn=False
res@xyExplicitLegendLabels=legendl
end

undef("log_xy")
procedure log_xy(res)
begin
res@trYLog=True
res@trXLog=True
res@tmYLMinorPerMajor=8
res@tmXBMinorPerMajor=8
end

;;;;;;;; Plotting section   ;;;;;;;;;;;
begin

;get environment variables
 ExpN      = getenv("EXPN")
 legendl   = getenv("hor_res")
 FileN     = getenv("ErrorFile")
 Root      = stringtoint(getenv("Root"))
 bis       = stringtoint(getenv("Bis"))
 Var       = getenv("VarN")
 pname     = getenv("PNAME")
 PltFmt    = getenv("PFMT")
 ConfigStr = getenv("CONFIG_STRING")
 diro      = getenv("DIRO")

;number of columns in Error File - neglect header
 ncols=7
 nheadl=1

;read error file
 table=readAsciiTable(FileN,ncols,"double",nheadl)

;get data for (time) x-axis
 thour=table(:,0)
 tday=thour/24.0
 tdims=dimsizes(thour) 

;initialize arrays and get data for y-axis
 l_n=new((/3,tdims/),double)
; if Var .eq. PS, then the absolute errors are considered
 if (Var .eq. "PS") then
  l_n(0,:)=table(:,1)
  l_n(1,:)=table(:,3)
  l_n(2,:)=table(:,5)
  punit=" (hPa)"
 else
  l_n(0,:)=table(:,2)
  l_n(1,:)=table(:,4)
  l_n(2,:)=table(:,6)
  punit=" norm"
 end if

;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;; plots ;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

 wks = gsn_open_wks(PltFmt,diro+pname)
 plot = new(3,graphic)

do i=1,3                                    ;loop over different norms

 if (i .ne. 3)
   title=Var+"  l"+i+" error"+punit
 else
   title=Var+" l_inf error"+punit
 end if

 ResXY=True
 ResXY@tiXAxisString      = "Time (days)"
 ResXY@tiYAxisString      = title  
;ResXY@trYMinF            = 0.0
;ResXY@trYMaxF            = 0.5 
 ResXY@tiMainString       = " "
 ResXY@gsnRightString     = ConfigStr
 ResXY@gsnDraw            = False
 ResXY@gsnFrame           = False
 ResXY@xyLineColors       =(/"blue"/)
if (Var .eq. "PS") then
 ResXY@trYMinF           = 0.0
 ResXY@trYMaxF           = 0.5
end if

 labelres(ResXY,legendl)

 plt_ind=i-1  
 plot(plt_ind)=gsn_csm_xy(wks,tday,l_n(plt_ind,:),ResXY)
  
end do                                      ;end norm loop
 
 gsn_panel(wks,plot,(/1,3/),False)          ;draw panel plot 

end
