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

undef("BWdefine_zon_draw_strings")
 procedure BWdefine_zon_draw_strings (nlati,Code,comment,wks,plot)
;************************************************
;  Define min max mean for black and white plot
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

undef("plot_cont_zon_bw")
 procedure plot_cont_zon_bw (plot,wks,var,nlati, BWres_in,pltName,\
                mainTitle,subTitle,comment)
;************************************************

 local plot, BWres

begin
    BWres=BWres_in
    BWres@gsnFrame              = False           ; don't advance frame yet
    BWres@gsnDraw               = False             ; don't draw plot yet 
    BWres@tiMainString          = mainTitle
    BWres@tiYAxisString         = " "             ; Y-Axis no label
    BWres@gsnCenterString       = subTitle
    BWres@gsnLeftString         = " " 
    BWres@gsnRightString        = " " 
    BWres@gsnContourZeroLineThicknessF = 2. 	  ; thickness of zero contour
    BWres@gsnContourNegLineDashPattern = 1        ; negative contours to dash 
    BWres@trYLog                = False           ; no logarithmic Y-axis
    BWres@tmXMajorGrid          = True                ; X Axis grid lines
    BWres@tmXMajorGridLineColor = "grey"          ; color X Axis grid lines
    BWres@tmYMajorGrid          = True            ; Y Axis grid lines
    BWres@tmYMajorGridLineColor = "grey"          ; color Y Axis grid lines
    BWres@tmYLMode = "Explicit"
    BWres@tmYLValues = (/100,200,300,400,500,600,700,800,900,1000/)
    BWres@tmYLLabels = (/100,200,300,400,500,600,700,800,900,1000/) 
    BWres@tmYRMode = "Automatic"

    BWres@cnInfoLabelOn        = False            ; turn off contour label box
    BWres@cnLineLabelFormat    = "#f"             ; no unnecessary 0

    BWres@vpWidthF  = 0.8		; set width of plot 
    BWres@vpHeightF = 0.5		; set height of plot


    plot = gsn_csm_pres_hgt(wks,var, BWres)  ; create plot

;***define statistical Values and draw Values, comment, pltdate
    BWdefine_zon_draw_strings (nlati,var,comment,wks,plot)

end
;***end procedure plot_cont_zon_bw

;************************************************
begin

  values    = asciiread("var.txt",8,"string")
  name = values(0)
  typ = values(1)
  run = values(2)
  meantime = values(3)
  comment = values(4)
  workdir = values(5)        ; plot & input directory
  eratime = values(6)
  nlati = stringtointeger(values(7))

;  print(eratime)
;  print(workdir)


;***define plot and open
   pltType = "pdf"                        ; x11, eps, ps, pdf, png
   pltPath = workdir
   pltName = pltPath+"/"+run+"_atm_3d_"+typ+"_linp"        ; plot name
   wks     = gsn_open_wks(pltType, pltName) 
;***define color map 
    gsn_define_colormap(wks,"amwg_blueyellowred") ; Define color map 
    missColor = NhlNewColor(wks,0.5,0.5,0.5)      ; add gray to color map 
;***plot array for panel
   listplots = new(14,graphic)
   noPlot=-1  ;number of listplots

;************************************************
;   Temperature code130 zonal CMOR: ta
;************************************************

    Cvar="ta"
    fili    = workdir+"/Ubusy_"+Cvar+".nc"
    
    mainTitle = "Temperature [C] "+typ
    subTitle  = run+meantime

if (isfilepresent(fili)) then

;***read code 130
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    BWres=True
    BWres@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
    BWres@cnMinLevelValF       =   -100.         ; set min contour level
    BWres@cnMaxLevelValF       =     60.         ; set max contour level
    BWres@cnLevelSpacingF      =      5.         ; set contour spacing

    plot_cont_zon_bw(listplots(noPlot), wks,var, nlati, BWres,pltName,\
                mainTitle,subTitle,comment )
  end if

;************************************************
;   T code130 zonal model-ERAin
;************************************************

    Cvar="ta"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
    
    mainTitle = "Temperature [C] "+typ
    subTitle  = run+meantime+" - Erain"+eratime

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
    res@cnLevels             = (/-16,-8,-4,-2,-1,1,2,4,8,16/)
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
    fili    = workdir+"/Ubusy_"+Cvar+".nc"
   
    mainTitle = "U [m/s] "+typ
    subTitle  = run+meantime

if (isfilepresent(fili)) then

;***read code 131
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    BWres                      = True             ; plot mods desired
    BWres@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
    BWres@cnMinLevelValF       =   -50.           ; set min contour level
    BWres@cnMaxLevelValF       =    50.           ; set max contour level
    BWres@cnLevelSpacingF      =     5.           ; set contour spacing

    plot_cont_zon_bw(listplots(noPlot), wks,var, nlati, BWres,pltName,\
                mainTitle,subTitle,comment )

  end if

;************************************************
;   U code131 zonal model-ERAin
;************************************************

    Cvar="ua"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
    
    mainTitle = "U [m/s] "+typ
    subTitle  = run+meantime+" - Erain"+eratime

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
    res@cnLevels             = (/-16,-8,-4,-2,-1,1,2,4,8,16/)
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
    fili    = workdir+"/Ubusy_"+Cvar+".nc"
   
    mainTitle = "V [m/s] "+typ
    subTitle  = run+meantime

if (isfilepresent(fili)) then
;***read code 132
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    BWres                      = True             ; plot mods desired
    BWres@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
    BWres@cnMinLevelValF       =   -10.           ; set min contour level
    BWres@cnMaxLevelValF       =    10.           ; set max contour level
    BWres@cnLevelSpacingF      =     0.5          ; set contour spacing

    plot_cont_zon_bw(listplots(noPlot), wks,var, nlati, BWres,pltName,\
                mainTitle,subTitle,comment )
  end if

;************************************************
;   V code132 zonal model-ERAin
;************************************************

    Cvar="va"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
    
    mainTitle = "V [m/s] "+typ
    subTitle  = run+meantime+" - Erain"+eratime

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
    res@cnLevels             = (/-1.,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1./)
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
    fili    = workdir+"/Ubusy_"+Cvar+".nc"
   
    mainTitle = "Relative humidity [%] "+typ
    subTitle  = run+meantime

if (isfilepresent(fili)) then
;***read 
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    BWres                      = True             ; plot mods desired
    BWres@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
    BWres@cnMinLevelValF       =   10.           ; set min contour level
    BWres@cnMaxLevelValF       =    90.           ; set max contour level
    BWres@cnLevelSpacingF      =     10.          ; set contour spacing

    plot_cont_zon_bw(listplots(noPlot), wks,var, nlati, BWres,pltName,\
                mainTitle,subTitle,comment )
  end if

;************************************************
;   hur code157 Relative humidity [%] zonal model-ERAin
;************************************************

    Cvar="hur"
    fili    = workdir+"/Ubusy_d"+Cvar+".nc"
    
    mainTitle = "Relative humidity [%] "+typ
    subTitle  = run+meantime+" - Erain"+eratime

if (isfilepresent(fili)) then

;***read 
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    res                      = True             ; plot mods desired
    res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    res@cnLevels             = (/-30.,-20.,-10.,-5.,-2.,2.,5.,10.,20.,30./)
    res@cnFillColors         = (/2,4,7,8,9,0,11,12,14,15,17/)
    res@cnMissingValFillColor = missColor 
                       
    plot_cont_zon ( listplots(noPlot),wks,var, nlati, res,pltName,\
                mainTitle,subTitle,comment )

    delete(res@cnLevels)
    delete (res@cnFillColors)
  end if

;************************************************
;   Cloud cover code223 zonal
;************************************************

    Cvar="cl"
    fili    = workdir+"/Ubusy_"+Cvar+".nc"
   
    mainTitle = "Cloud cover [%] "+typ
    subTitle  = run+meantime

if (isfilepresent(fili)) then
;***read 
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    BWres                      = True             ; plot mods desired
    BWres@tmYLMode = "Explicit"
    BWres@tmYLValues = (/100,200,300,400,500,600,700,800,900,1000/)
    BWres@tmYLLabels = (/100,200,300,400,500,600,700,800,900,1000/) 

    BWres@cnInfoLabelOn        = False            ; turn off contour label box
    BWres@cnLineLabelFormat    = "#f"             ; no unnecessary 0
    BWres@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
    BWres@cnMinLevelValF       =   5.             ; set min contour level
    BWres@cnMaxLevelValF       =  95.             ; set max contour level
    BWres@cnLevelSpacingF      =   5.             ; set contour spacing


    plot_cont_zon_bw(listplots(noPlot), wks,var, nlati, BWres,pltName,\
                mainTitle,subTitle,comment )
  end if

;************************************************
;   Cloud water code153 zonal   CMOR: CLW
;************************************************

    Cvar="clw"
    fili    = workdir+"/Ubusy_"+Cvar+".nc"
    
    mainTitle = "Cloud water [mg/kg] "+typ
    subTitle  = run+meantime

if (isfilepresent(fili)) then

;***read code 
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    BWres                      = True             ; plot mods desired
    BWres@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    BWres@cnLevels             = (/1,2,5,10,20,50/)

    plot_cont_zon_bw(listplots(noPlot), wks,var, nlati, BWres,pltName,\
                mainTitle,subTitle,comment )
    delete(BWres@cnLevels)
  end if

;************************************************
;   Cloud ice code154 zonal   CMOR: cli
;************************************************

    Cvar="cli"
    fili    = workdir+"/Ubusy_"+Cvar+".nc"
   
    mainTitle = "Cloud ice [mg/kg] "+typ
    subTitle  = run+meantime

if (isfilepresent(fili)) then

;***read code 154
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    delete(var@long_name)                  ; don't draw long_name by default
    delete(var@units)                      ; don't draw unit by default
;    printVarSummary(var) 

;***number of listplots
    noPlot=noPlot+1

;***create plot
    BWres                      = True             ; plot mods desired
    BWres@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    BWres@cnLevels             = (/1,2,5,10,20,50/)

    plot_cont_zon_bw(listplots(noPlot), wks,var, nlati, BWres,pltName,\
                mainTitle,subTitle,comment )

    delete(BWres@cnLevels)
  end if

;******************************************************************
;   Cloud ice + Cloud water code153 + code154 zonal  CMOR: clw + cli
;******************************************************************

    Cvar="clw_cli"
    fili    = workdir+"/Ubusy_"+Cvar+".nc"
    
    mainTitle = "Cloud ice + Cloud water [mg/kg] "+typ
    subTitle  = run+meantime

if (isfilepresent(fili)) then

;***read code
    f       = addfile (fili , "r") 
    var     = f->clw(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var) 

;***number of listplots
    noPlot=noPlot+1

;***create plot
    BWres                      = True             ; plot mods desired
    BWres@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    BWres@cnLevels             = (/1,2,5,10,20,50/)

    plot_cont_zon_bw(listplots(noPlot), wks,var, nlati, BWres,pltName,\
                mainTitle,subTitle,comment )

    delete(BWres@cnLevels)
  end if
;************************************************
;   Specific humidity code133 zonal CMOR: hus
;************************************************

    Cvar="hus"
    fili    = workdir+"/Ubusy_"+Cvar+".nc"

    mainTitle = "Specific humidity [g/kg] "+typ
    subTitle  = run+meantime

if (isfilepresent(fili)) then
;***read code 
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    BWres                      = True             ; plot mods desired
    BWres@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
    BWres@cnLevels             = (/.001,.002,.005,.01,.02,.05,.1,.2,.5,1,2,4,8,16/)

    plot_cont_zon_bw(listplots(noPlot), wks,var, nlati, BWres,pltName,\
                mainTitle,subTitle,comment )

    delete(BWres@cnLevels)
  end if
;************************************************
;   psi (mastrfu) zonal computed from va
;************************************************

    Cvar="mastrfu"
    fili    = workdir+"/Ubusy_"+Cvar+".nc"
   
    mainTitle = "Psi [10~S~9~N~kg/s]  hPa "+typ
    subTitle  = run+meantime

if (isfilepresent(fili)) then
;***read code 132
    f       = addfile (fili , "r") 
    var     = f->$Cvar$(time|0,lev|:,lat|:,lon|0)
    print("ploted: "+mainTitle)
;   printVarSummary(var)

;***number of listplots
    noPlot=noPlot+1

;***create plot
    BWres                      = True             ; plot mods desired
    BWres@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
    BWres@cnMinLevelValF       =   -200.           ; set min contour level
    BWres@cnMaxLevelValF       =    200.           ; set max contour level
    BWres@cnLevelSpacingF      =     10.          ; set contour spacing

    plot_cont_zon_bw(listplots(noPlot), wks,var, nlati, BWres,pltName,\
                mainTitle,subTitle,comment )
  end if

print("number of plots:  "+noPlot)
;--KMF prepare plots for panel   
   pltres          =  True
   pltres@gsnPanelScalePlotIndex = 1   ; 2rd plot in list, which is the biggest   
                    
   gsn_panel(wks, (/listplots(0),listplots(1),\
                    listplots(2),listplots(3),\
                    listplots(4),listplots(5)/), (/3,2/), pltres)  
   gsn_panel(wks, (/listplots(6),listplots(7),\
                    listplots(8),listplots(9),\
                    listplots(10),listplots(11)/), (/3,2/), pltres) 
   gsn_panel(wks, (/listplots(12),listplots(13)/), (/3,1/), pltres)


end