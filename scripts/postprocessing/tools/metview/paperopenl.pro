; ICON
;
; ---------------------------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ---------------------------------------------------------------

pro paperopenl
;-----------------------------------------------------------------------;
; Open postscript landscape page.
;
; Martin Koehler, 2-20-95
;-----------------------------------------------------------------------;

 set_plot,'ps',/copy
 device,filename='plot.ps',/color,/landscape,   $
;   /inches, xsize=10.5,ysize=7.7, xoffset=0.25,yoffset=10.75
    xsize=28.5,ysize=19, xoffset=0.5,yoffset=29.5
;, xoffset=1,yoffset=28.5
 !p.charsize=1.5
 !p.charthick=2
 !p.thick=2
 !x.thick=2
 !y.thick=2

;!p.font=1                          ;True-type fonts
;device, set_font='Times', /tt_font ;see /usr/local/rsi/idl/resource/fonts

 !p.font=0      ;default (hardward postscript font Helvetia, problem with math)
                ;- better for Illustrator (?)
;device, /times ;Times postscript font 
;device, /avantgarde
 device, /helvetica
;device, /palatino

;!p.font=-1     ;vector drawn fonts (not so nice, but compatible with math)
;xyouts,0.5,0.5,/normal,'!6'   ; to set vector drawn font to #6

 end
