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



undef("define_zon_draw_strings")
 procedure define_zon_draw_strings (nlati,Code,comment,wks)
;************************************************
;  Define min max mean
;  Draw text on plot using NDC coordinates.
;************************************************
local txres

begin
;***compute weights for mean
  nlats=nlati/2
  gs=doubletofloat(gaus(nlats))
  gweight=gs(:,1)
;***compute mean 
  meanV=dim_avg_wgt(Code,gweight,0)

  minV=min(Code)
  maxV=max(Code)
  meanV=dim_avg(meanV)

  statV="min="+minV+"   mean="+meanV(0)+"   max="+maxV

  pltdate = systemfunc("date +%d-%m-%Y")
  print("actual_date: "+pltdate)

  txres               = True                    ; text mods desired
  txres@txFontHeightF = 0.014                   ; font smaller. default big
  gsn_text_ndc(wks,pltdate,0.93,0.08,txres)
  gsn_text_ndc(wks,comment,0.17,0.08,txres)
  txres@txFontHeightF = 0.016 
  gsn_text_ndc(wks,statV,0.6,0.22,txres)
end
;***end procedure 



undef("plot_cont_zon_log")
 procedure plot_cont_zon_log (wks,var,nlati, res_in,pltName,\
                mainTitle,subTitle,comment)
;************************************************

 local plot, res

begin
    GridColorIndex = NhlNewColor(wks,0.4,0.4,0.4) ;MajorGridLineColor grey
    res                      = True             ; plot mods desired
    res@gsnFrame             = False            ; don't advance frame yet
    res@tiMainString         = mainTitle
    res@tiYAxisString        = " "              ; Y-Axis no label
    res@gsnLeftString        = " "  
    res@gsnRightString       = " "  
    res@gsnCenterString      = subTitle

    res@trYLog                = True            ; logarithmic Y-axis
    res@tmXMajorGrid          = True            ; X Axis grid lines
    res@tmXMajorGridLineColor = GridColorIndex  ; color X Axis grid lines
    res@tmYMajorGrid          = True            ; Y Axis grid lines
    res@tmYMajorGridLineColor = GridColorIndex  ; color Y Axis grid lines
    res@tmYLMode = "Explicit"
    res@tmYLValues = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/)
    res@tmYLLabels = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000/) 
    res@tmYRMode = "Automatic"
    res@tmYROn               = False     ; Turn off y-axis (right) major ticks
    res@tmXTOn               = True  
    res=res_in
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.07           ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars
    
    res@vpWidthF  = 0.78		; set width of plot 
    res@vpHeightF = 0.5		; set height of plot

    plot = gsn_csm_pres_hgt(wks,var, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_zon_draw_strings (nlati,var,comment,wks)

    frame(wks)                                    ; now advance frame

end
;***end procedure plot_cont_zon_log

;************************************************
begin


  values    = asciiread("var.txt",7,"string")
  typ = values(0)
  Arun = values(1)
  Brun = values(2)
  meantime = values(3)
  comment = values(4)
  workdir = values(5)        ; plot & input directory
  nlati = stringtointeger(values(6))

  print(meantime)
  print(workdir)
  print(nlati)

;***define plot
   pltType = "png"                        ; x11, eps, ps, pdf, png
   pltPath = workdir+"/"+Arun+"_"+Brun+"_atm_log"



;***define plot
   pltType = "png"                        ; x11, eps, ps, pdf, png
;   pltPath = workdir+"/"+run+"_atm_"
   pltPath = workdir+"/"+Arun+"_"+Brun+"_atm_"


;************************************************
;   Temperature code130 zonal  CMOR: ta
;************************************************

    Cvar="ta"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
    
    mainTitle = "Temperature [C] "+typ
    subTitle  = meantime

if (isfilepresent(fili)) then

;***read code 130
    f       = addfile (fili , "r") 
    d130    = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(d130)

;***open plot
    pltName = pltPath+Cvar+"_"+typ+"_logp_cont"
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    missColor = NhlNewColor(wks,0.5,0.5,0.5)      ; add gray to color map 

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels"   
    res@cnLevels             = (/-8,-4,-2,-1,-.5,.5,1,2,4,8/)
    res@cnFillColors         = (/2,4,6,7,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor 
                      
    plot_cont_zon_log ( wks,d130, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if


;************************************************
;   U code131 zonal 
;************************************************

    Cvar="ua"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
    
    mainTitle = "U [m/s] "+typ
    subTitle  = meantime

if (isfilepresent(fili)) then

;***read code 131
    f       = addfile (fili , "r") 
    d131     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(d131)

;***open plot
    pltName = pltPath+ Cvar+"_"+typ+"_logp_cont"
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    missColor = NhlNewColor(wks,0.5,0.5,0.5)      ; add gray to color map 

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-8,-4,-2,-1,-.5,.5,1,2,4,8/)
    res@cnFillColors         = (/2,4,6,7,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor
                     
    plot_cont_zon_log ( wks,d131, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if



;************************************************
;   V code132 zonal CMOR: va
;************************************************

    Cvar="va"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
    
    mainTitle = "V [m/s] "+typ
    subTitle  = meantime

if (isfilepresent(fili)) then

;***read code 132
    f       = addfile (fili , "r") 
    d132     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(d132)

;***open plot
    pltName = pltPath+ Cvar+"_"+typ+"_logp_cont"
    wks     = gsn_open_wks(pltType, pltName) 

    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    missColor = NhlNewColor(wks,0.5,0.5,0.5)      ; add gray to color map 

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor
                     
    plot_cont_zon_log ( wks,d132, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

end
