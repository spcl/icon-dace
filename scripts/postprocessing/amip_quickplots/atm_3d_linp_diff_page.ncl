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
 procedure define_zon_draw_strings (nlati,Code,comment,wks,plot)
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

;--KMF create text object and add it to the plot
  txres                  =  True                 ; text mods desired
  txres@txFontHeightF    =  0.014                 ;-- smaller font
  txid = gsn_create_text(wks, statV, txres)
  
  amres                  =  True
  amres@amParallelPosF   =  0.0                  ;-- center position
  amres@amOrthogonalPosF =  0.7                  ;-- move downward
  amres@amSide           = "Bottom"
  annoid = gsn_add_annotation(plot, txid, amres)
end
;***end procedure 

undef("plot_cont_zon")
 procedure plot_cont_zon (plot,wks,var,nlati, res_in,pltName,\
                mainTitle,subTitle,comment)
;************************************************

 local plot, res

begin

    res=res_in
    res@gsnFrame             = False            ; don't advance frame yet
    res@gsnDraw              = False             ; don't draw plot yet
    res@tiMainString         = mainTitle
    res@tiYAxisString        = " "              ; Y-Axis no label
    res@gsnCenterString      = subTitle
    res@gsnLeftString        = " " 
    res@gsnRightString       = " " 
    res@trYLog   = False                        ; no logarithmic Y-axis
    res@tmXMajorGrid      = True                ; X Axis grid lines
    res@tmXMajorGridLineColor = "grey"  ; color X Axis grid lines
    res@tmYMajorGrid      = True                ; Y Axis grid lines
    res@tmYMajorGridLineColor = "grey"  ; color Y Axis grid lines
    res@tmYLMode = "Explicit"
    res@tmYLValues = (/100,200,300,400,500,600,700,800,900,1000/)
    res@tmYLLabels = (/100,200,300,400,500,600,700,800,900,1000/) 
    res@tmYRMode = "Automatic"
    res=res_in
    res@cnFillOn             = True               ; turn on color fill
    res@cnLinesOn            = True               ; turn off contour lines
    res@gsnSpreadColors      = False              ; use full colormap
    res@pmLabelBarOrthogonalPosF = 0.07           ; move farther to plot
    res@lbOrientation        = "horizontal"       ; vertical label bars

    res@vpWidthF  = 0.8		; set width of plot 
    res@vpHeightF = 0.5		; set height of plot

    plot = gsn_csm_pres_hgt(wks,var, res)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    define_zon_draw_strings (nlati,var,comment,wks,plot)

end
;***end procedure plot_cont_zon




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

;***define plot and open
   pltType = "pdf"                        ; x11, eps, ps, pdf, png
   pltPath = workdir
   pltName = pltPath+"/"+Arun+"_"+Brun+"_atm_3d_"+typ+"_linp"        ; plot name
   wks     = gsn_open_wks(pltType, pltName) 
;***define color map 
    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    missColor = NhlNewColor(wks,0.5,0.5,0.5)      ; add gray to color map 
;***plot array for panel
   listplots = new(14,graphic)
   noPlot=-1  ;number of listplots



;************************************************
;   T code130 zonal model-ERAin
;************************************************

    Cvar="ta"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
    
    mainTitle = "Temperature [C] "+typ
    subTitle  = meantime

if (isfilepresent(fili)) then
;***read code 131
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = (/-8,-4,-2,-1,-.5,.5,1,2,4,8/)
    res@cnFillColors         = (/2,4,6,7,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor
                       
    plot_cont_zon ( listplots(noPlot),wks,var, nlati, res,pltName,\
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
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/-8,-4,-2,-1,-.5,.5,1,2,4,8/)
    res@cnFillColors         = (/2,4,6,7,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor
                       
    plot_cont_zon ( listplots(noPlot),wks,var, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if



;************************************************
;   V code132 zonal 
;************************************************

    Cvar="va"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
    
    mainTitle = "V [m/s] "+typ
    subTitle  = meantime

if (isfilepresent(fili)) then

;***read code 132
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

 
;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor 
                       
    plot_cont_zon (listplots(noPlot), wks,var, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if


;************************************************
;   hur code157 Relative humidity [%] zonal
;************************************************

    Cvar="hur"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
   
    mainTitle = "Relative humidity [%] "+typ
    subTitle  = meantime

if (isfilepresent(fili)) then
;***read code
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

 
;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/-20,-10.,-5.,-2.,-1,1,2.,5.,10.,20/)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor 
                       
    plot_cont_zon (listplots(noPlot), wks,var, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if


;************************************************
;   Cloud cover code223 zonal
;************************************************

    Cvar="cl"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
   
    mainTitle = "Cloud cover [%] "+typ
    subTitle  = meantime

if (isfilepresent(fili)) then
;***read code
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

 
;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             =  (/-6,-4,-2,-1,-.5,.5,1,2,4,6/)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor 
                       
    plot_cont_zon (listplots(noPlot), wks,var, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

;************************************************
;   Cloud water code153 zonal   CMOR: CLW
;************************************************

    Cvar="clw"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
    
    mainTitle = "Cloud water [mg/kg] "+typ
    subTitle  = meantime

if (isfilepresent(fili)) then
;***read code
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

 
;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             =  (/-10.,-5.,-2.,-1.,-.5,.5,1.,2.,5.,10./)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor 
                       
    plot_cont_zon (listplots(noPlot), wks,var, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if


;************************************************
;   Cloud ice code154 zonal   CMOR: cli
;************************************************

    Cvar="cli"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
   
    mainTitle = "Cloud ice [mg/kg] "+typ
    subTitle  = meantime
if (isfilepresent(fili)) then
;***read code
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

 
;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             =  (/-5.,-2.,-1.,-.5,-.2,.2,.5,1.,2.,5./)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor 
                       
    plot_cont_zon (listplots(noPlot), wks,var, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

;******************************************************************
;   Cloud ice + Cloud water code153 + code154 zonal  CMOR: clw + cli
;******************************************************************

    Cvar="clw_cli"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
    
    mainTitle = "Cloud ice + Cloud water [mg/kg] "+typ
    subTitle  = meantime

if (isfilepresent(fili)) then

;***read code
    f       = addfile (fili , "r") 
    var     = f->clw(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var) 

;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             =  (/-10.,-5.,-2.,-1.,-.5,.5,1.,2.,5.,10./)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor 
                       
    plot_cont_zon (listplots(noPlot), wks,var, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

;************************************************
;   Specific humidity code133 zonal CMOR: hus
;************************************************

    Cvar="hus"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"

    mainTitle = "Specific humidity [g/kg] "+typ
    subTitle  = meantime

if (isfilepresent(fili)) then

;***read code
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var) 

;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             =  (/-2.,-1.,-0.5,-.2,-0.1,0.1,.2,0.5,1.,2./)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor 
                       
    plot_cont_zon (listplots(noPlot), wks,var, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

;***number of listplots
    noPlot=noPlot+1
print("number of plots:  "+noPlot)

;--KMF prepare plots for panel   
   pltres          =  True
   pltres@gsnPanelScalePlotIndex = 1   ; 2rd plot in list, which is the biggest   
                    
   gsn_panel(wks, (/listplots(0),listplots(1),\
                    listplots(2),listplots(3),\
                    listplots(4),listplots(5)/), (/3,2/), pltres)  
   gsn_panel(wks, (/listplots(6),listplots(7),\
                    listplots(8),listplots(9)/), (/3,2/), pltres) 

end
