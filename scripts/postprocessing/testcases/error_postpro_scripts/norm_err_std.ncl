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
;
; First version calculating normalized errors.
; First version for the steady state case. The Icon 
; file contain the true solution as initial data
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("f_vweight")
 function f_vweight(hyai[*]:double,hybi[*]:double)
 local eta_i,vweight,vdim,p0
 begin
  p0=100000 ;(Pa)  
  vdim=dimsizes(hyai)-1
  eta_i=new(vdim+1,double)
  eta_i(0:vdim)=hyai(0:vdim)/p0 + hybi(0:vdim)
  vweight=new(vdim,double)
  vweight(0:vdim-1)= eta_i(1:vdim)-eta_i(0:vdim-1)
  return(vweight)
 end
;------------------

undef("cal_err_3D")
 function cal_err_3D(VarT[*][*]:numeric,VarI[*][*]:numeric,tweight[*][*]:numeric,hour)
 local m_err,n_l1,n_l2,n_linf,VarT2,l1_err,l2_err,linf_err, \
       l1_norm_err,l2_norm_err,linf_norm_err, VarDif, VarDifa, VarDif2
 begin
  m_err=new(7,"double")

  ; normalization values
  n_l1=new(1,double)
  n_l2=new(1,double)
  n_linf=new(1,double)
  n_l1=sum(abs(VarT)*tweight)
  VarT2=VarT*VarT
  n_l2=sqrt(sum(VarT2*tweight))
  delete(VarT2)
  n_linf=max(abs(VarT))


  ;calculation of errors
   VarDif=VarI(:,:)-VarT(:,:)
   VarDifa=abs(VarDif)

   l1_err=sum(VarDifa*tweight)
   linf_err=max(VarDifa)
   delete(VarDifa)
   
   VarDif2=VarDif*VarDif
   printVarSummary(VarDif2)
   printVarSummary(tweight)
   l2_err=sum(VarDif2*tweight)
   delete(VarDif2)
   delete(VarDif)
   l2_err=sqrt( l2_err)

   if (n_l1 .gt. 1.e-5) then
    l1_norm_err=l1_err/n_l1
    ;print(l1_norm_err)
    ;print(l1_err)
    ;print(n_l1)
   else
    l1_norm_err=l1_err
    ;print("no")
   end if
   if (n_l2 .gt. 1.e-5) then
    l2_norm_err=l2_err/n_l2
   else
    l2_norm_err=l2_err
   end if
   if (n_linf .gt. 1.e-5) then
    linf_norm_err=linf_err/n_linf
   else
    linf_norm_err=linf_err
   end if

   m_err(0)=hour
   m_err(1)=l1_err
   m_err(2)=l1_norm_err
   m_err(3)=l2_err
   m_err(4)=l2_norm_err
   m_err(5)=linf_err
   m_err(6)=linf_norm_err
  
  return(m_err)
 end
;------------------

undef("cal_err_2D")
 function cal_err_2D(VarT[*]:numeric,VarI[*]:numeric,tweight[*],hour)
 local m_err,n_l1,n_l2,n_linf,VarT2,l1_err,l2_err,linf_err, \
       l1_norm_err,l2_norm_err,linf_norm_err, VarDif, VarDifa, VarDif2
       
 begin
  m_err=new(7,"double")

  ; normalization values
  n_l1=new(1,double)
  n_l2=new(1,double)
  n_linf=new(1,double)
  n_l1=sum(abs(VarT)*tweight)
  VarT2=VarT*VarT
  n_l2=sqrt(sum(VarT2*tweight))
  n_linf=max(abs(VarT))


  ;calculation of errors

   VarDif=VarI(:)-VarT(:)
   VarDifa=abs(VarDif)
   VarDif2=VarDif*VarDif
 
   l1_err=sum(VarDifa*tweight)
   l2_err=sum(VarDif2*tweight)
   l2_err=sqrt( l2_err)
   linf_err=max(VarDifa)

   if (n_l1 .gt. 1.e-35) then
    l1_norm_err=l1_err/n_l1
    ;print(l1_norm_err)
    ;print(l1_err)
    ;print(n_l1)
   else
    l1_norm_err=l1_err
    ;print("no")
   end if
   if (n_l2 .gt. 1.e-35) then
    l2_norm_err=l2_err/n_l2
   else
    l2_norm_err=l2_err
   end if
   if (n_linf .gt. 1.e-35) then
    linf_norm_err=linf_err/n_linf
   else
    linf_norm_err=linf_err
   end if

   m_err(0)=hour
   m_err(1)=l1_err
   m_err(2)=l1_norm_err
   m_err(3)=l2_err
   m_err(4)=l2_norm_err
   m_err(5)=linf_err
   m_err(6)=linf_norm_err
 
  return(m_err)
 end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
begin

rad2deg = 45./atan(1.)   ; radians to degrees

GridFileN=getenv("GridFileN")
IconFileN=getenv("IconFileN")
VarN=getenv("VarN")
VarD=getenv("VarD")
ErrorFile=getenv("ErrorFile")

IFile=addfile(IconFileN,"r")

; get information from the ICON file

f_att=getvaratts(IFile)
 if (isatt(IFile, "nroot"))then
  Root=IFile@nroot
 else 
  if (isatt(IFile, "grid_root")) then
   Root=IFile@grid_root
  else
   Root=" "
  end if
 end if
 if (isatt(IFile, "start_lev"))then
  Bis=IFile@start_lev
 else 
  if (isatt(IFile, "grid_level")) then
   Bis=IFile@grid_level
  else
   Bis=" "
  end if
 end if

if (isatt(IFile, "nlev"))then
  Nvl=IFile@nlev
 else
  Nvl=31
end if

dtime_io=-1.
;printVarSummary(Nvl)
if (isatt(IFile, "dt_data"))then
  dtime_io=IFile@dt_data
 else
  if (isatt(IFile, "io_nml:dt_data"))then
   att="io_nml:dt_data"
   dtime_io=IFile@$att$
  else
   dtime_io=-1.
  end if
end if

print("dtime_io="+dtime_io)

;printVarSummary(dtime)
 if (isatt(IFile, "i_cell_type"))then
  CellType=IFile@i_cell_type
 else
  if (isatt(IFile, "run_ctl:i_cell_type"))then
   att="run_ctl:i_cell_type"
   CellType=IFile@$att$
  else
   CellType="3"
  end if
 end if
 if (isatt(IFile, "dtime"))then
  dtime=IFile@dtime
 else
  dtime=0.
 end if
 if (isatt(IFile, "out_expname"))then
  ExpN=IFile@out_expname
 else
  if (isatt(IFile, "io_ctl:out_expname"))then
   att="io_nml:out_expname"
   ExpN=IFile@$att$
  else
   ExpN="JWs"
  end if
 end if
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Read variable true and model values
;printVarSummary(IFile)
;printVarSummary(dtime)
;print("hallo")
;VarI=IFile->$VarN$
;print("hallo")
;printVarSummary(VarI)

ndims=VarD
if (ndims .eq. 2) then
  VarT=IFile->$VarN$(0,:)   ; (time, cell)
else
 if (ndims .eq. 3) 
  VarT=IFile->$VarN$(0,:,:)   ; (time, mlev, cell)
 ;Read hybrid parameters from ICON File to get the vertical weights

  Ai=IFile->hyai     ;from 1 to Nvl+1 interfaces layers (from 0 to Nvl)
  Bi=IFile->hybi
  vweight=f_vweight(Ai,Bi)
 end if
end if

; set units of pressure to hPa for VarT
if (isatt(VarT,"units")) then
   var_u=VarT@units
   if (var_u .eq. "Pa") then
    VarT=VarT/100.
    VarT@units="hPa"
   end if
end if

;;;;;;;;;
;calculate total weights and normalized them
;;;;;;;;;

;Read horizontal (area) weights from the grid
Grid=addfile(GridFileN,"r")
if (CellType .eq. 3) then
 hweight=Grid->cell_area
 gridtype="tri"
else
 if (CellType .eq. 6) then
  hweight=Grid->dual_area
  gridtype="hex"
 else
  print("wrong CellType")
  exit
 end if
end if
cdim=dimsizes(hweight)
if (ndims .eq. 2) then
 ;2D variable, only horizontal weights
 tweight=hweight/sum(hweight)
else
 if (ndims .eq. 3) then
  ;3D variable, create 2D tweight
  ldim=dimsizes(IFile->mlev)
  tweight=new((/ldim,cdim/),"double")
  tot_w=new(1,double)
  tot_w=0.
  do ilev=0,ldim-1
   do icell=0,cdim-1
    tweight(ilev,icell)=vweight(ilev)*hweight(icell)
    tot_w=tot_w+tweight(ilev,icell)
   end do
  end do
  ;normalized total weights
  tweight=tweight/tot_w
 end if
end if
delete(Grid)
delete(hweight)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; create matrix variable to store the errors
; and calculate errors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

tdim=dimsizes(IFile->time)
mat_err=new((/tdim,7/),"double")

; calculate the errors
if (ndims .eq. 2) then
  do it=0,tdim-1 
   hour=it*dtime_io/3600.
   VarIt=IFile->$VarN$(it,:)
   ; set units of pressure to hPa
   if (isatt(VarIt,"units")) then
    var_u=VarIt@units
    if (var_u .eq. "Pa") then
     VarIt=VarIt/100.
     VarIt@units="hPa"
    end if
   end if
   mat_err(it,:)=cal_err_2D(VarT,VarIt,tweight,hour)
   delete(VarIt)
  end do
else
 if (ndims .eq. 3) then
  print(ndims)
  do it=0,tdim-1 
   hour=it*dtime_io/3600.
   VarIt=IFile->$VarN$(it,:,:)
   ; set units of pressure to hPa
   if (isatt(VarIt,"units")) then
    var_u=VarIt@units
    if (var_u .eq. "Pa") then
     VarIt=VarIt/100.
     VarIt@units="hPa"
    end if
   end if
   mat_err(it,:)=cal_err_3D(VarT,VarIt,tweight,hour)
   delete(VarIt)
  end do
 end if
end if


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;output of errors to a file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

opt=True
opt@fout=ErrorFile
opt@title="     time(hour)         l1_err    l1_norm_err         l2_err    l2_norm_err       linf_err  linf_norm_err  " + VarN +" (" + VarT@units+ ") "   +" R"+Root+"B0"+Bis+" " +gridtype
fmt="7g15.4"
write_matrix(mat_err,fmt,opt)

end