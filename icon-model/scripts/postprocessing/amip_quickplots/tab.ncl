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
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

print (pltdir)
print (quelle)

  datfil = pltdir+"/busy_atm_phy.nc"
  varfil = quelle+"/VARS.txt"

  varnumbers=stringtoint(systemfunc("wc -l " +varfil+ "|cut -f1 -d' '"))
print(varnumbers)
   maxnumber=varnumbers+2

 Vmean  = new((/maxnumber/),double)
 Vindex  = new((/maxnumber/),integer)
 Vname  = new((/maxnumber/),string)
 Vstring  = new((/maxnumber/),string)
 Vcode  = new((/maxnumber/),integer)
 Vindex = -1
 meanV  = 0.
    
 strs = asciiread(varfil,-1,"string")
 delim = " "
 nfields = str_fields_count(strs(0), delim)   ; Count the fields separated
 print(nfields)                               ; by one or more spaces.

 field = 1
 Varname = str_get_field(strs, field, delim)
; print(Varname)              
 field = 2
 Code = stringtointeger(str_get_field(strs, field, delim))
; print(Code)              
 field = 3
 Factor = stringtofloat(str_get_field(strs, field, delim))
; print(Factor)              
 field = 4
 Offset = stringtofloat(str_get_field(strs, field, delim))
; print(Offset)              

  GridInfoFileName = systemfunc("echo $GrdInfoFile")
;  GridInfoFileName ="/pool/data/ICON/grids/private/r2b4_amip/r2b4_amip.nc"
  print(GridInfoFileName)

  GridInfoFile= addfile( GridInfoFileName, "r" )

  rad2deg = 45./atan(1.)            ; radians to degrees

  x = GridInfoFile->clon *rad2deg   ; cell center, lon
  y = GridInfoFile->clat *rad2deg   ; cell center, lat
  cArea = GridInfoFile->cell_area   ; area of grid cell


  File    = addfile (datfil , "r")  

  i=-1
  do k=0,varnumbers-1
    if (isfilevar(File,Varname(k))) then
        i = i+1
        Var = File->$Varname(k)$(time|0,ncells|:)
        Var=(Var*Factor(k))+Offset(k)
        Vmean(i)=dim_avg(Var*cArea)/dim_avg(cArea)
        Vindex(i)=k
        Vname(i)=Varname(k)
        Vcode(i)=Code(k)
;print("RESULT  "+Varname(k)+"   "+Vmean(i)+" vindex  "+Vindex(i))
;print("      min: "+min(Var)+ "      max: "+max(Var))
    end if
  end do

NoVar=i
print(NoVar)

  if (isfilevar(File,"prlr").and.isfilevar(File,"prls")) then
      NoVar=NoVar+1
      Var1 = File->prlr(time|0,ncells|:)
      Var2 = File->prls(time|0,ncells|:)
      Var= (Var1+Var2)*86400.
      Vmean(NoVar)=dim_avg(Var*cArea)/dim_avg(cArea)
      Vcode(NoVar)=142
      Vname(NoVar)="aprl"
   end if

  if (isfilevar(File,"prcr").and.isfilevar(File,"prcs")) then
      NoVar=NoVar+1
      Var1 = File->prcr(time|0,ncells|:)
      Var2 = File->prcs(time|0,ncells|:)
      Var= (Var1+Var2)*86400.
      Vmean(NoVar)=dim_avg(Var*cArea)/dim_avg(cArea)
      Vcode(NoVar)=143
      Vname(NoVar)="aprc"
   end if

print (sprintf("%9.4f",(Vmean(:NoVar)))+"  "+sprinti("%3.0i",(Vcode(:NoVar)))+"  "+(Vname(:NoVar)))
Vstring(:NoVar)= (sprintf("%9.4f",(Vmean(:NoVar)))+"  "+sprinti("%3.0i",(Vcode(:NoVar)))+"  "+(Vname(:NoVar)))
asciiwrite ("global.txt" , Vstring(:NoVar) )  

exit
  
end
