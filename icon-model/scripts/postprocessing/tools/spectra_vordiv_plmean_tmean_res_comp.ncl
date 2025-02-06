; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;================================================;
;  spectra_vordi.ncl
;
;	Example ncl script to produce the
;	kinetic energy plot.
; 	Created by James Kent, Christiane Jablonowski
;	and Paul Ullrich (University of Michigan).  
;
;   modified by levi silvers
;   modified by pilar ripodas: this version uses div and vorticity 
;                  instead of u and v
;                  it uses the internal vertical interpolation to 
;                  pressure levels that is in the model
;                  it calculates the spectra for different days and pressure 
;                   levels and calculates the mean
;                   
;       to use this script simply type
;       >ncl spectra_vordiv.ncl
;

;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
; ================================================;


begin

  re=6371229.    ;Radius of the Earth
;=================================================;
; open file and read in data 
; GIVEN AS (time,lev,lat,long) from 0 to n-1
;=================================================;


; directories containing experiement data
      dirbase = "/e/uhome/mripodas/icon-dev/experiments/"
      exper =(/"real_spectra_R2B04","real_spectra_R2B05","real_spectra_R2B06","real_spectra_R2B07"/)
      rsltn=(/"R2B04L90","R2B05L90","R2B06L90","R2B07L90"/)
      wn_nom=(/125,250,500,1000/)     ; 1000 for R2B07

     ; exper =(/"real_spectra_div_avg_R2B04","real_spectra_div_avg_R2B05","real_spectra_div_avg_R2B06","real_spectra_div_avg_R2B07"/)
     ; rsltn=(/"R2B04L90","R2B05L90","R2B06L90","R2B07L90"/)

     ; wn_nom=(/125,250,500,1000/)     ; 1000 for R2B07

     ; exper =(/"real_spectra_div_avg_R2B04","real_spectra_div_avg_R2B05","real_spectra_div_avg_R2B06","real_spectra_div_avg_R2B07"/)
     ; rsltn=(/"R2B04L90","R2B05L90","R2B06L90","R2B07L90"/)
     ; wn_nom=(/125,250,500,1000/)     ; 1000 for R2B07

     ; exper =(/"real_spectra_R2B06","real_spectra_div_avg_R2B06"/)
     ; rsltn=(/"R2B06L90","R2B06L90_div_avg"/)
     ; wn_nom=(/500,500/)     ; 1000 for R2B07
     ; wn_nom=(/119,239,449/)     ; 1000 for R2B07


      ;exper =(/"real_spectra_R2B07"/)
      ;rsltn=(/"R2B07L90"/)
      ;wn_nom=(/1000/)     ; 1000 for R2B07

      ;exper =(/"real_spectra_R2B04"/)
      ;rsltn=(/"R2B04L90"/)
      ;wn_nom=(/125/)     ; 1000 for R2B07

      n_exp=dimsizes(exper)
      n_res=dimsizes(rsltn)





 legendl=new(n_res,string)

      ; here I should check n_res=n_exp

; set day and level to calculate the spectra
        f_day=10  ;first day
        l_day=15  ;last day
        n_days=l_day-f_day+1
        ;ilev=9    ; 9 is 500 hPa - 11 is 700 hPa
        f_lev=7    ; first level
        l_lev=9      ; last level
        n_lev=l_lev-f_lev+1
 
; The previous values could be passed via environment variables 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  spc=new ( (/n_lev,n_days,n_exp,10000/), "double", 1e-30 )
  spc_plm=new ( (/n_days,n_exp,10000/), "double", 1e-30 )

  spc_plm_tm=new ( (/n_exp,10000/), "double", 1e-30 ) 
  nlat_e=new((/n_exp/),"integer")
  print ("dimsizes of spc are: "+dimsizes(spc))
  spc    = spc@_FillValue
  spc_plm = spc@_FillValue
  spc_plm_tm = spc@_FillValue
  x      = ispan(0,10000-1,1)*1.0
 ; spc!1  = "x"
 ; x!0    = "x"
  max_nlat=0
  do iexp=0,n_exp-1

;        FileN=dirbase+exper(iexp)+"/"+exper(iexp)+"_DOM01_PL_0001.nc"
         FileN=dirbase+exper(iexp)+"/"+"real_spectra"+"_DOM01_PL_0001.nc"
;        FileN=dirbase+exper(iexp)+"/"+"real_spectra_div_avg"+"_DOM01_PL_0001.nc"
  

        print(FileN)
        ;print(exper(iexp))

	f= addfile(FileN,"r")
        print ("dimsizes of f are: "+dimsizes(f))


; parameters

	lat  = f->lat
  	lon  = f->lon
	lev  = f->lev
	nlat   = getfilevardimsizes(f, "lat" )
	nlon   = getfilevardimsizes(f, "lon" )
	nlev   = getfilevardimsizes(f, "lev" )
  
        nlat_e(iexp)=nlat
        if (nlat .gt. max_nlat) then
          max_nlat=nlat
        end if
   ; it is assumed that the output frequency is 1 day


        myf_lev=lev(f_lev)/100.
        myl_lev=lev(l_lev)/100.
        print(myf_lev)
        print(myl_lev)

; For the plot
    output_file="no_divavg_spectra_divvor_res_comp_dmean"+f_day+"-"+l_day+"_plmean"+myf_lev+"-"+myl_lev
    ;wks  = gsn_open_wks("X11","icon-test410KE")   	; output using X11
    wks  = gsn_open_wks("ps",output_file)   	; output using eps
    ;wks  = gsn_open_wks("png","icon-test410KE")   	; output using png

  do iday=0,n_days-1
  do ilev=0,n_lev-1

    
        day=f_day+iday
        mlev=f_lev+ilev    ;model level

        mylevel=lev(mlev)/100.
        print(mylevel)
        print(iday)


     ; Check if Gaussian Grid

	gauGrid = True
	if (abs(lat(1)-lat(0)) .eq. abs(lat(nlat/2+1)-lat(nlat/2)) ) then
       		gauGrid = False
	end if

; getfilevartypes() this function returns the type of each var name listed, U here
;
; new (dimension_sizes,vartype,parameter) new is used to create a new NCL variable 
; and gives the new variable dimension sizes
;
; ispan(start,finish,stride) This function returns an array of values beginning at 
;                            start and ending at finish, stride must be positive 
;
   	;spc    = new ( (/nlat/), getfilevartypes(f,"div"), 1e20 )
 
        DIV   = f->div(day,mlev,:,:)
     	VOR   = f->vor(day,mlev,:,:)

        print(mlev)
        print(iexp)


      
        ;printVarSummary(DIV)
        ;printVarSummary(VOR)
 
     	print ("computing scalar spherical harmonic transform...")
     	if (gauGrid) then
         	abdiv = shagC(DIV(:,:))
         	abvor = shagC(VOR(:,:));  ; Gauss grid
     	else
         	abdiv = shaeC(DIV(:,:))
         	abvor = shaeC(VOR(:,:));   ; Equal spaced grid with both poles
     	end if
        ;printVarSummary(abdiv)
        ;printVarSummary(abvor) 

        ;a triangular truncation can be performed
        ;triangular truncation of coefficients
        ;trunc=159
        ;abdiv=tri_trunC(abdiv, trunc) 
        ;abvor=tri_trunC(abvor, trunc)

     	divr  = abdiv(0,:,:)                		; real coef  (nlat,N)
     	divi  = abdiv(1,:,:)                		; imaginary  (nlat,N)
     	pwr = (divr^2 + divi^2)         		; (nlat,N)  array

        vorr  = abvor(0,:,:)                		; real coef  (nlat,N)
     	vori  = abvor(1,:,:)                		; imaginary  (nlat,N)
     	pwr = pwr + (vorr^2 + vori^2)        	; (nlat,N)  array

     	do n1=1,nlat-1
        	spc(ilev,iday,iexp,n1) = (/ pwr(n1,0) /)
       		do m=1,n1
          		spc(ilev,iday,iexp,n1) = spc(ilev,iday,iexp,n1) + 2.*pwr(n1,m)
       		end do
                spc(ilev,iday,iexp,n1) = spc(ilev,iday,iexp,n1)/4./int2flt(n1)/(n1+1.)
     	end do
        spc(ilev,iday,iexp,1:nlat-1)=spc(ilev,iday,iexp,1:nlat-1)*re*re
        ;print(spc)
   ; spc() is the array that the results need to be saved into.  do we need to create more
; spc() arrays or simply increase the dimensionality?  i also need to check how to overlay
; multiple lines on one ncl plot.   

     delete(divr)
     delete(divi)
     delete(vorr)
     delete(vori)
     delete(pwr)
     delete(DIV)
     delete(VOR)
     delete(abdiv)
     delete(abvor)
    end do  ; ilev
     do n1=1,nlat-1
        spc_plm(iday,iexp,n1)=sum(spc(0:n_lev-1,iday,iexp,n1))/n_lev
     end do
    end do  ;iday
     do n1=1,nlat-1
        spc_plm_tm(iexp,n1)=sum(spc_plm(0:n_days-1,iexp,n1))/n_days
     end do
     spc_plm_tm(iexp,wn_nom(iexp)+1:) = spc@_FillValue
     delete(f)
     delete(lat)
     delete(lon)
     delete(lev)
    end do  ;iexp

	;------------------------------------------------------------
 	; 				plot
 	;------------------------------------------------------------
     	lthick = 2.0
  
        FontHeight = 0.02
 
     	res          = True
     	res@xyXStyle = "Log"
     	res@xyYStyle = "Log"
     	res@tmXBMinorPerMajor = 8   
     	res@tmXTMinorPerMajor = 8   
     	res@tmYRMinorPerMajor = 8   
     	res@tmYLMinorPerMajor = 8 

        res@gsnLeftStringFontHeightF  = FontHeight
        res@gsnRightStringFontHeightF = FontHeight
        res@tiXAxisFontHeightF        = FontHeight
        res@tiYAxisFontHeightF        = FontHeight 

	pltTitle="Real case. DIV/VOR. Model vert. intp. " ; Plot title if required 
	pltTitle=""
     
     	res@tiMainString      = pltTitle
     	res@gsnCenterString   = "KE at "+myf_lev+"-"+myl_lev+" hPa  days "+f_day+"-"+l_day
     	res@tiXAxisString     = "Spherical wavenumber"
     	res@tiYAxisString     = "Kinetic energy (m~S~2~N~ s~S~-2~N~)"
     	res@gsnDraw           = False
     	res@gsnFrame          = False
     	res@xyLineThicknessF  = lthick
   
                           ; create custom color map
     	colors = (/"white","black","OrangeRed4","Red3",     \
     	           "red1","Red3","Red4","Blue",             \
     	"OrangeRed2","OrangeRed3","OrangeRed4","DarkOrange",\
     	"red1","red2","red3","red4", "OrangeRed1",          \
     	"DarkTurquoise","Turquoise","Cyan","blue",          \
     	"gold","greenyellow","GoldenRod","DarkGoldenRod1",  \
     	"lightgoldenrod1","lightyellow","lightgoldenrodyellow","khaki1"/)
     
	
	gsn_define_colormap(wks, colors)
       
     	res@trXMinF = 1       
     	res@trXMaxF = 1500   ; !!This must be larger than max_nlat    
     	res@trYMinF = 1e-6      
     	res@trYMaxF = 1e3   
        ;res@gsnMaximize    = True

        ;legend

        res@pmLegendDisplayMode="Always"
        res@pmLegendSide="Top"
        res@pmLegendZone=1
        res@pmLegendParallelPosF=0.2
        res@pmLegendOrthogonalPosF=0.8
        res@pmLegendWidthF=0.1
        res@pmLegendHeightF=0.2
        res@lgLabelFontHeightF=0.015
        res@lgPerimOn=False
        res@xyExplicitLegendLabels=rsltn
        ; end legend

     	res@xyDashPattern = 0
     	res@xyLineColors  = (/"Red","BLue","orange","brown"/)
     	plot    = gsn_csm_xy (wks,x(1:),spc_plm_tm(:,1:),res) 		; create plot    
    
     	xref    = (/10.,200./)
     	spc_ref = 5*(xref/15)^(-3) 				; make it 6 @ k=15
     	spc_ref = .06*(xref/100.)^(-3) 				; make it 0.8 @ k=40
     
     	polyres                  = True
     	polyres@gsLineThicknessF = lthick
     	polyres@gsLineColor      = "black"
     	gsn_polyline(wks,plot,xref,spc_ref,polyres)             ; add polyline

     	xref    = (/50.,400./)
     	spc_ref = .06*(xref/100.)^(-5./3.) 				; make it 0.8 @ k=40
     
     	polyres                  = True
     	polyres@gsLineThicknessF = lthick
     	polyres@gsLineColor      = "black"
     	gsn_polyline(wks,plot,xref,spc_ref,polyres)             ; add polyline
        draw(plot)
     	frame(wks)                                     		; now advance
     


end