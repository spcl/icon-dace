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
;  spectra_uv.ncl
;
;	Example ncl script to produce the
;	kinetic energy plot.
; 	Created by James Kent, Christiane Jablonowski
;	and Paul Ullrich (University of Michigan).  
;
;   modified by levi silvers
;   modified by pilar ripodas: this version 
;                  it uses the internal vertical interpolation to 
;                  pressure levels that is in the model
;  
;       to use this script simply type
;       >ncl spectra_uv.ncl
;

;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
; ================================================;


begin

;=================================================;
; open file and read in data 
; GIVEN AS (time,lev,lat,long) from 0 to n-1
;=================================================;

; specify the test # and the resolution as well as other help strings
        test="real"
        rsltn="R2B04L90"
; dir is the directory contained experiement data
        dirbase = "/e/uhome/mripodas/icon-dev/experiments/"
        exper ="real_R2B04_spectra_latlon"
        dir = dirbase+exper+"/"
; replace your models file name below
; below is the file name to be read
	;f= addfile("toexps/"+dir+name1+test+"_"+rsltn+"_DOM01_PL_0001.nc","r")
        Name="real_R2B04_spectra_latlon"
        FileN=dir+Name+"_DOM01_PL_0001.nc"

        print(FileN)

; The previous values could be passed via environment variables 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	f= addfile(FileN,"r")
        print ("dimsizes of f are: "+dimsizes(f))


; parameters

	lat  = f->lat
  	lon  = f->lon
	lev  = f->lev
	nlat   = getfilevardimsizes(f, "lat" )
	nlon   = getfilevardimsizes(f, "lon" )
	nlev   = getfilevardimsizes(f, "lev" )

   ; it is assumed that the output frequency is 1 day

    day=(/1,4,5,10,15/)
    ;day=(/4/)    
    ndays=dimsizes(day)

    level=(/7,8,9,11/)
    nlevels=dimsizes(level)

; For the plot
    output_file="spectra_uv_"+rsltn   
    ;wks  = gsn_open_wks("X11","icon-test410KE")   	; output using X11
    wks  = gsn_open_wks("ps",output_file)   	; output using eps
    ;wks  = gsn_open_wks("png","icon-test410KE")   	; output using png



    do jday=0,ndays-1

     do jlev=0,nlevels-1

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
   	spc    = new ( (/nlat/), getfilevartypes(f,"div"), 1e20 )
        print ("dimsizes of spc are: "+dimsizes(spc))
   	spc    = spc@_FillValue
   	x      = ispan(0,nlat-1,1)*1.0
   	spc!0  = "x"
   	x!0    = "x"
 
        U   = f->u(day(jday),level(jlev),:,:)
     	V   = f->v(day(jday),level(jlev),:,:)

        mylevel=lev(level(jlev))/100.
        print(mylevel)
      
        ;printVarSummary(U)
        ;printVarSummary(V)
 
     	print ("computing scalar spherical harmonic transform...")
     	if (gauGrid) then
         	ab = vhagC(U(:,:),V(:,:))  ; Gauss grid
     	else
         	ab = vhaeC(U(:,:),V(:,:));  ; Equal spaced grid with both poles
     	end if
        ;printVarSummary(ab) 

        ;a triangular truncation can be performed
        ;triangular truncation of coefficients
        ;trunc=159
        ;ab(0:1,:,:)=tri_trunC(ab(0:1,:,:), trunc) 
        ;ab(2:3,:,:)=tri_trunC(ab(2:3,:,:), trunc) 
 
     	cr  = ab(0,:,:)                		; real coef  (nlat,nlat)
     	ci  = ab(1,:,:)                		; imaginary  (nlat,nlat)
     	pwr = (cr^2 + ci^2)/2.         		; (nlat,nlat)  array

     	cr  = ab(2,:,:)                		; real coef  (nlat,nlat)
     	ci  = ab(3,:,:)                		; imaginary  (nlat,nlat)
     	pwr = pwr + (cr^2 + ci^2)/2.         	; (nlat,nlat)  array

     	do n1=0,nlat-1
        	spc(n1) = (/ pwr(n1,0) /)
       		do m=1,n1
          		spc(n1) = spc(n1) + 2.* pwr(n1,m)
       		end do
        	spc(n1) = 0.5*spc(n1)
     	end do
   ;  print(spc)


   ; spc() is the array that the results need to be saved into.  do we need to create more
; spc() arrays or simply increase the dimensionality?  i also need to check how to overlay
; multiple lines on one ncl plot.   
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

	pltTitle="Real case. U/V. Model vert. intp. "+rsltn ; Plot title if required 
    
     	res@tiMainString      = pltTitle
     	res@gsnCenterString   = "KE at "+mylevel+" hPa  at day "+day(jday)
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
     	res@trXMaxF = 500       
     	res@trYMinF = 1e-6      
     	res@trYMaxF = 1e3   
        res@gsnMaximize    = True


     	res@xyDashPattern = 0
     	res@xyLineColors  = (/"Red"/)
     	plot    = gsn_csm_xy (wks,x(1:),spc(1:),res) 		; create plot    
    
     	xref    = (/15.,300./)
     	spc_ref = 5*(xref/15)^(-3) 				; make it 6 @ k=15
     
     	polyres                  = True
     	polyres@gsLineThicknessF = lthick
     	polyres@gsLineColor      = "black"
     	gsn_polyline(wks,plot,xref,spc_ref,polyres)             ; add polyline
        draw(plot)
     	frame(wks)                                     		; now advance
     
     end do
    end do

end






