; ICON
;
; ---------------------------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ---------------------------------------------------------------

pro scores, direxp, dirref, expnum, expref, inidate, step1, step2, nfor, levtype, stat
;------------------------------------------------------------
; Plot scores from ICON experiments.
;
; run as: .run scores paperopenl paperclose
;         scores,'/scratch/ms/de/deia/icon/dei2/084/201201/metplots',
;                '/scratch/ms/de/deia/icon/dei2/083/201201/metplots',
;                'dei2_084','dei2_083','201201','24','31','ml', 'mean'
;
; Martin Koehler, Sep 2012
;------------------------------------------------------------

print, 'Arguments: ', direxp, ' ', dirref, ' ', expnum, ' ', expref, ' ', $
  inidate, ' ', step1, ' ', step2, ' ', nfor, ' ', levtype, ' ', stat

nlev=90
nlevpl=25
nlevzl=25
nlevlnd=8
nword=19

stepx=24
step1=fix(step1)
step2=fix(step2)


if stat eq 'rms' then begin
  stat_txt   = '_rms'
  diff_txt   = 'rms'
  stat_title = '   -- forecast RMS errors --   '

  amp_ml_rms  = 1000.0    ;* 240/step   ;100.0
  amp_pl_bias = 500.0     ;* 240/step   ;100.0
  amp_pl_rms  = 1000.0    ;* 240/step   ;100.0
  amp_zl_rms  = 1000.0    ;* 240/step   ;40.0
  amp_sfc_bias= 500.0     ;* 240/step   ;12.0
  amp_sfc_rms = 1000.0    ;* 240/step   ;30.0
  refer       = 0.001*0.5 ;* step/240   ;0.01*0.5 = 1%

endif else begin

  stat_txt   = ''
  diff_txt   = 'diff'
  stat_title = '   -- climate mean errors --   '

  amp_ml_rms  = 100.0               ;40.0
  amp_pl_bias = 50.0                ;40.0
  amp_pl_rms  = 100.0               ;40.0
  amp_zl_rms  = 100.0               ;40.0
  amp_sfc_bias= 50.0                ;12.0
  amp_sfc_rms = 100.0               ;30.0
  refer       = 0.01*0.5            ;0.01*0.5 = 1%
endelse

refer_txt = strmid(strcompress(string(refer/0.5*100.0)),1,4) + ' Percent'

title = expnum+' vs. '+expref+stat_title+ $
   nfor+' forecasts from '+inidate+' + '+strcompress(string(step1),/remove_all)+'h'
if step1 ne step2 then begin
  title = title + '-' + strcompress(string(step2),/remove_all)+'h'
endif

if step1 eq step2 then begin
  file1 = direxp+"/scores_"+expnum+stat_txt+".txt"
  file2 = dirref+"/scores_"+expref+stat_txt+".txt"
endif else begin
  file1 = direxp+"/scores_"+expnum+stat_txt+"_all.txt"
  file2 = dirref+"/scores_"+expref+stat_txt+"_all.txt"
endelse
print, 'experiment file: ', file1
print, 'reference file:  ', file2
print, ''
close,/all
openr, 1, file1
openr, 2, file2
spawn,'wc -l '+file1, nlines1
spawn,'wc -l '+file2, nlines2
nlines1=strsplit(nlines1,/extract)
nlines2=strsplit(nlines2,/extract)

var1=make_array(nword,nlines1(0),/string)
var2=make_array(nword,nlines2(0),/string)
tmp=''
for nn=1,nlines1(0) do begin
  readf, 1, tmp
  var1(*,nn-1) = strsplit(tmp,/extract)
end
for nn=1,nlines2(0) do begin
  readf, 2, tmp
  var2(*,nn-1) = strsplit(tmp,/extract)
end

var2d  = ['TOT_PREC' , 'TCC'      , 'HCC'     , 'MCC'     , 'LCC'     , $
          'TQV'      , 'TQC'      , 'TQI'     , $
          'ACCLHFL_S', 'ACCSHFL_S', $
          'ACCSOB_S' , 'ACCTHB_S' , 'ACCSOB_T', 'ACCTHB_T', $
          'PS'       , 'T_2M'     , 'T_G'     , 'U_10M'   , 'V_10M']
nvar2d = n_elements(var2d)

var2obs  = ['ACCSOB_T', 'ACCTHB_T', 'ACCLHFL_S', 'ACCSHFL_S', 'T_2M'     , 'T_G'     , 'SP_10M' ]
obsname  = ['CERES'   , 'CERES'   , 'WHOI'     , 'WHOI'     , 'WHOI'     , 'WHOI'    , 'WHOI'   ]
nvar2obs = n_elements(var2obs)

var3lnd  = ['W_SO', 'T_SO' ]
nvar3lnd = n_elements(var3lnd)

var3d  = ['T', 'U', 'V', 'P', 'QV', 'CC', 'QC', 'QI']   ;, 'Q2', 'CRWC','CSWC'];'QTVAR', 'O3','P''QV', 'QC', 'QI'
nvar3d = n_elements(var3d)

var3pl = ['T', 'U', 'QV', 'FI']   ;, 'U'
nvar3pl= n_elements(var3pl)

var3zl = ['T', 'U', 'V', 'P', 'QC', 'QI']   ;, 'Z'
nvar3zl= n_elements(var3zl)

pres_l = ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100', '150', '200', '250', '300', '400', '500', $
          '600', '700', '800', '850', '900', '925', '950', '1000']

paperopenl
loadct,13


;--------------------------------------------------------------------------------------------

CASE levtype OF

 'ml': BEGIN

  xyouts, 0.15, 0.97, /normal, charsize=1.0, title+'   '+levtype

  for nn=1,nvar3d do begin
   for nk=1,nlev do begin
   ;print,''
   ;print,'searching for variable ', var3d(nn-1)
    mean1_mn=0 & mean2_mn=0 & rms1_mn=0 & rms2_mn=0
    for nstep=step1,step2,stepx do begin
      ind1 = where(var1(0,*) eq var3d(nn-1) and var1(2-1,*) eq levtype and var1(8-1,*) eq nstep and var1(10-1,*) eq nfor and var1(4-1,*) eq nk and var1(19-1,*) eq diff_txt)
      ind2 = where(var2(0,*) eq var3d(nn-1) and var2(2-1,*) eq levtype and var2(8-1,*) eq nstep and var2(10-1,*) eq nfor and var2(4-1,*) eq nk and var2(19-1,*) eq diff_txt)
      ind1=ind1(n_elements(ind1)-1)
      ind2=ind2(n_elements(ind2)-1)
      if ind1 eq -1 or ind2 eq -1 then continue
      mean1_mn= float(var1(16-1,ind1)) + mean1_mn
      mean2_mn= float(var2(16-1,ind2)) + mean2_mn
      rms1_mn = float(var1(18-1,ind1)) + rms1_mn
      rms2_mn = float(var2(18-1,ind2)) + rms2_mn
    end
    mean1 = mean1_mn / ( (step2-step1)/stepx + 1)  ; 10 days
    mean2 = mean2_mn / ( (step2-step1)/stepx + 1)
    rms1  = rms1_mn  / ( (step2-step1)/stepx + 1)
    rms2  = rms2_mn  / ( (step2-step1)/stepx + 1)
    bias = abs(mean1) / ( abs(mean1) + abs(mean2) )    ; [0,1], 0: mean1 is much better, 1: mean2 is much better: 0.5 same
    rms  =     rms1   / (     rms1   +     rms2   )
   ;print, var3d(nn-1), ' L', nk, ' exp',expref, ' bias:', mean2, ' rms:', rms2
    if ( abs(mean1) + abs(mean2) ) lt 1e-8 then begin
      bias = 0.5
    endif
    if (     rms1  +      rms2   ) lt 1e-8 then begin
      rms  = 0.5
    endif
  
    xxx =-0.11+nn*0.125
    yyy =0.93
    vert=0.01
  
    xyouts, xxx+0.035, yyy+0.015, var3d(nn-1), /normal, charsize=0.7                      ;title each column: variable name
    header = 'level    rms      new'
    xyouts, xxx, yyy, header, /normal, charsize=0.7                                       ;header (level, rms, new)
    if stat eq 'rms' then begin
      err1=mean1
      err2=mean2
      err =bias
    endif else begin
      err1=rms1
      err2=rms2
      err =rms
    endelse
    xoff=0.0
    scale = (err-0.5) * (err gt 0.5)
    scale = min([scale, 6.0/amp_ml_rms])
    xr= [0.0,0.01,0.01,0.0 ,0.0] * amp_ml_rms * scale + xxx + xoff
    yr= [0.0,0.0 ,0.01,0.01,0.0] * 0.8            + yyy
    polyfill, xr   +0.04, yr   -nk*vert,           /normal, color=240                     ;positive
    scale = (0.5-err) * (err lt 0.5)
    scale = min([scale, 6.0/amp_ml_rms])
    xr=-[0.0,0.01,0.01,0.0 ,0.0] * amp_ml_rms * scale + xxx + xoff
    polyfill, xr     +0.04, yr     -nk*vert,       /normal, color=160                     ;negative
    plots,    xr(3:4)+0.04, [yr(3)+0.002,yr(4)-0.002]-nk*vert, /normal                    ;zero
    xr= [-0.01,0.01]             * amp_ml_rms * refer + xxx + xoff
    plots,    xr(0:1)+0.04, yr(0:1)-nk*vert+0.005, /normal                                ;reference 1%

    xyouts,   xr(0)       , yr(0)-nk*vert, strcompress(string(nk)), /normal, charsize=0.5 ;level
   ;xyouts,   xr(0)+0.07  , yr(0)-nk*vert, err1,   /normal, charsize=0.4                  ;new error value
    xyouts,   xr(0)+0.09  , yr(0)-nk*vert, string(err1,format='(g8.3)'), /normal, charsize=0.4, alignment=1.0 ;new error value
  end
  end
  
  xr = [-0.01,0.0] * amp_ml_rms * refer
  plots, xr(0:1)+0.50, 0.015+[0.0,0.0],  /normal                                          ;reference 1% (bottom page)
  xyouts,xr(1)  +0.51, 0.010, refer_txt, /normal, charsize=0.6                            ;text refernce value

 END


;--------------------------------------------------------------------------------------------

 'pl': BEGIN

  xyouts, 0.15, 0.97, /normal, charsize=1.0, title+'   '+levtype

  for nn=1,nvar3pl do begin
   for nk=1,nlevpl do begin
   ;print,''
   ;print,'searching for variable ', var3pl(nn-1)
    mean1_mn=0 & mean2_mn=0 & rms1_mn=0 & rms2_mn=0
    for nstep=step1,step2,stepx do begin
      ind1 = where(var1(0,*) eq var3pl(nn-1) and var1(2-1,*) eq levtype and var1(8-1,*) eq nstep and var1(10-1,*) eq nfor and var1(4-1,*) eq nk and var1(19-1,*) eq diff_txt)
      ind2 = where(var2(0,*) eq var3pl(nn-1) and var2(2-1,*) eq levtype and var2(8-1,*) eq nstep and var2(10-1,*) eq nfor and var2(4-1,*) eq nk and var2(19-1,*) eq diff_txt)
      ind1=ind1(n_elements(ind1)-1)
      ind2=ind2(n_elements(ind2)-1)
      if ind1 eq -1 or ind2 eq -1 then continue
      mean1_mn= float(var1(16-1,ind1)) + mean1_mn
      mean2_mn= float(var2(16-1,ind2)) + mean2_mn
      rms1_mn = float(var1(18-1,ind1)) + rms1_mn
      rms2_mn = float(var2(18-1,ind2)) + rms2_mn
    end
    mean1 = mean1_mn / ( (step2-step1)/stepx + 1)  ; 10 days
    mean2 = mean2_mn / ( (step2-step1)/stepx + 1)
    rms1  = rms1_mn  / ( (step2-step1)/stepx + 1)
    rms2  = rms2_mn  / ( (step2-step1)/stepx + 1)
    if stat eq 'rms' then begin
      bias = abs(mean1) / ( abs(mean1) + abs(mean2) ) ; [0,1], 0: mean1 is much better, 1: mean2 is much better: 0.5 same
    endif else begin
      bias = ( abs(mean1) - abs(mean2) ) / (     rms1   +     rms2   ) + 0.5
    endelse
    rms  =     rms1   / (     rms1   +     rms2   )
   ;print, var3pl(nn-1), ' L', nk, ' exp',expref, ' bias:', bias2, 'rms:', rms
    if ( abs(mean1) + abs(mean2) ) lt 1e-8 then begin
      bias = 0.5
    endif
    if (     rms1  +      rms2   ) lt 1e-8 then begin
      rms  = 0.5
    endif
  
    xxx =-0.23+nn*0.25
    yyy =0.90
    vert=0.033
  
    if stat eq 'rms' then begin
      headernew = 'new'
      header    = 'level      rms        ref'
    endif else begin
      headernew = 'new                                    new'
      header    = 'level     bias        ref          level       rms        ref'
    endelse
    xyouts, xxx+0.067, yyy+0.015, headernew, /normal, charsize=0.7,alignment=0.0     ;header (level, rms, new)
    xyouts, xxx-0.007, yyy,       header,    /normal, charsize=0.7                   ;header (level, rms, new)
 
    for nerr=1,2 do begin
      if nerr eq 2 and stat eq 'rms' then continue
      if nerr eq 1 then begin
        err1=mean1
        err2=mean2
        err =bias
        xoff=0.0
        amp = amp_pl_bias
      endif else begin
        err1=rms1
        err2=rms2
        err =rms
        xoff=0.125
        amp = amp_pl_rms
      endelse
      xyouts, xxx+0.035+xoff, yyy+0.025, var3pl(nn-1), /normal, charsize=0.7     ;title each column: variable name
      scale = (err-0.5) * (err gt 0.5)
      scale = min([scale, 6.0/amp])
      xr= [0.0,0.01,0.01,0.0 ,0.0] * amp * scale + xxx + xoff
      yr= [0.0,0.0 ,0.01,0.01,0.0] * 0.8            + yyy
      polyfill, xr   +0.04, yr   -nk*vert,              /normal, color=240       ;positive
      scale = (0.5-err) * (err lt 0.5)
      scale = min([scale, 6.0/amp])
      xr=-[0.0,0.01,0.01,0.0 ,0.0] * amp * scale + xxx + xoff
      polyfill, xr     +0.04, yr     -nk*vert,           /normal, color=160      ;negative
      plots,    xr(3:4)+0.04, [yr(3)+0.004,yr(4)-0.004]-nk*vert,      /normal    ;zero
      xr= [-0.01,0.01]             * amp * refer + xxx + xoff
      plots,    xr(0:1)+0.04, yr(0:1)-nk*vert+0.005, /normal                     ;reference 1%
  
      xyouts,   xr(0)       , yr(0)-nk*vert, pres_l(nk-1),/normal, charsize=0.7  ;level
      xyouts,   xr(0)+0.09  , yr(0)-nk*vert+0.007, string(err1,format='(g8.3)'), /normal, charsize=0.5, alignment=1.0 ;new error value
      xyouts,   xr(0)+0.09  , yr(0)-nk*vert-0.007, string(err2,format='(g8.3)'), /normal, charsize=0.5, alignment=1.0 ;old error value
      end
   end
  end

  xr = [-0.01,0.0] * amp_pl_bias * refer
  plots, xr(0:1)+0.40, 0.015+[0.0,0.0],  /normal                                 ;reference 1% (bottom page)
  xyouts,xr(1)  +0.41, 0.010, refer_txt+' bias', /normal, charsize=0.6           ;text refernce value

  if stat ne 'rms' then begin
  xr = [-0.01,0.0] * amp_pl_rms * refer
  plots, xr(0:1)+0.60, 0.015+[0.0,0.0],  /normal                                 ;reference 1% (bottom page)
  xyouts,xr(1)  +0.61, 0.010, refer_txt+' rms', /normal, charsize=0.6            ;text refernce value
  endif

 END


;--------------------------------------------------------------------------------------------

 'zl': BEGIN

  xyouts, 0.15, 0.97, /normal, charsize=1.0, title+'   '+levtype

  for nn=1,nvar3zl do begin
   for nk=1,nlevzl do begin
   ;print,''
   ;print,'searching for variable ', var3zl(nn-1)
    ind1 = where(var1(0,*) eq var3zl(nn-1) and var1(2-1,*) eq levtype and var1(8-1,*) eq step2 and var1(10-1,*) eq nfor and var1(4-1,*) eq nk and var1(19-1,*) eq diff_txt)
    ind2 = where(var2(0,*) eq var3zl(nn-1) and var2(2-1,*) eq levtype and var2(8-1,*) eq step2 and var2(10-1,*) eq nfor and var2(4-1,*) eq nk and var2(19-1,*) eq diff_txt)
    ind1=ind1(n_elements(ind1)-1)
    ind2=ind2(n_elements(ind2)-1)
    if ind1 eq -1 or ind2 eq -1 then continue
    mean1= float(var1(16-1,ind1))
    mean2= float(var2(16-1,ind2))
    rms1 = float(var1(18-1,ind1))
    rms2 = float(var2(18-1,ind2))
    bias = abs(mean1) / ( abs(mean1) + abs(mean2) )   ; [0,1], 0: mean1 is much better, 1: mean2 is much better: 0.5 same
    rms  =     rms1   / (     rms1   +     rms2   )
   ;print, var3zl(nn-1), ' L', nk, ' exp',expref, ' bias:', bias2, 'rms:', rms2
    if ( abs(mean1) + abs(mean2) ) lt 1e-8 then begin
      bias = 0.5
    endif
    if (     rms1  +      rms2   ) lt 1e-8 then begin
      rms  = 0.5
    endif
  
    xxx =-0.11+nn*0.14
    yyy =0.93
    vert=0.01
  
    xyouts, xxx+0.03, yyy+0.015, var3zl(nn-1), /normal, charsize=0.7                     ;title each column: variable name
    header = 'level    rms      new'
    xyouts, xxx, yyy, header, /normal, charsize=0.7                                      ;header (level, rms, new)
    if stat eq 'rms' then begin
      err1=mean1
      err2=mean2
      err =bias
    endif else begin
      err1=rms1
      err2=rms2
      err =rms
    endelse
    xoff=0.0
    scale = (err-0.5) * (err gt 0.5)
    scale = min([scale, 6.0/amp_zl_rms])
    xr= [0.0,0.01,0.01,0.0 ,0.0] * amp_zl_rms * scale + xxx + xoff
    yr= [0.0,0.0 ,0.01,0.01,0.0] * 0.8            + yyy
    polyfill, xr   +0.04, yr   -nk*vert,           /normal, color=240                     ;positive
    scale=(0.5-err) * (err lt 0.5)
    scale = min([scale, 6.0/amp_zl_rms])
    xr=-[0.0,0.01,0.01,0.0 ,0.0] * amp_zl_rms * scale + xxx + xoff
    polyfill, xr     +0.04, yr     -nk*vert,       /normal, color=160                     ;negative
    plots,    xr(3:4)+0.04, [yr(3)+0.002,yr(4)-0.002]-nk*vert, /normal                    ;zero
    xr= [-0.01,0.01]             * amp_zl_rms * refer + xxx + xoff
    plots,    xr(0:1)+0.04, yr(0:1)-nk*vert+0.005, /normal                                ;reference 1%

    xyouts,   xr(0)     ,   yr(0)-nk*vert, strcompress(string(nk)), /normal, charsize=0.5 ;level
    xyouts,   xr(0)+0.07,   yr(0)-nk*vert, err1,   /normal, charsize=0.4                  ;new error value
   end
  end
  
  xr = [-0.01,0.0] * amp_zl_rms * refer
  plots, xr(0:1)+0.50, 0.015+[0.0,0.0],  /normal                                          ;reference 1% (bottom page)
  xyouts,xr(1)  +0.51, 0.010, refer_txt, /normal, charsize=0.6                            ;text refernce value

 END


;--------------------------------------------------------------------------------------------

 'sfc': BEGIN

  xyouts, 0.15, 0.97, /normal, charsize=1.0, title+'   '+levtype
  
  ; surface parameters (e.g. fluxes) --------------------------------------------------------

  for nn=1,nvar2d do begin
   ;print,''
   ;print,'searching for variable ', var2d(nn-1)
    ind1 = where(var1(0,*) eq var2d(nn-1) and var1(2-1,*) eq levtype and var1(8-1,*) eq step2 and var1(10-1,*) eq nfor and var1(19-1,*) eq diff_txt)
    ind2 = where(var2(0,*) eq var2d(nn-1) and var2(2-1,*) eq levtype and var2(8-1,*) eq step2 and var2(10-1,*) eq nfor and var2(19-1,*) eq diff_txt)
    ind1=ind1(n_elements(ind1)-1)
    ind2=ind2(n_elements(ind2)-1)
    if ind1 eq -1 or ind2 eq -1 then continue
    mean1= float(var1(16-1,ind1))
    mean2= float(var2(16-1,ind2))
    rms1 = float(var1(18-1,ind1))
    rms2 = float(var2(18-1,ind2))
    if stat eq 'rms' then begin
      bias = abs(mean1) / ( abs(mean1) + abs(mean2) ) ; [0,1], 0: mean1 is much better, 1: mean2 is much better: 0.5 same
    endif else begin
      bias = ( abs(mean1) - abs(mean2) ) / (     rms1   +     rms2   ) + 0.5
    endelse
    rms  =     rms1   / (     rms1   +     rms2   )
   ;print, var2d(nn-1), ' exp',expref, ' bias:', bias2, 'rms:', rms2 
  
    xxx = 0.05
    yyy = 0.92
    vert= 0.045
  
    if stat eq 'rms' then begin
      header = 'variable                    rms                new,ref'
    endif else begin
      header = 'variable                    bias               new,ref                     rms               new,ref'
    endelse
    xyouts, xxx, yyy, header, /normal, charsize=0.8

    for nerr=1,2 do begin
      if nerr eq 2 and stat eq 'rms' then continue
      if nerr eq 1 then begin
        err1=mean1
        err2=mean2
        err =bias
        xoff=0.0
        amp =amp_sfc_bias
      endif else begin
        err1=rms1
        err2=rms2
        err =rms
        xoff=0.18
        amp =amp_sfc_rms
      endelse
      scale = (err-0.5) * (err gt 0.5)
      scale = min([scale, 15.0/amp])
      xr= [0.0,0.01,0.01,0.0 ,0.0] * amp * scale + xxx + xoff
      yr= [0.0,0.0 ,0.01,0.01,0.0] * 1.0            + yyy
      polyfill, xr   +0.12 , yr   -nn*vert,               /normal, color=240     ;positive
      scale = (0.5-err) * (err lt 0.5)
      scale = min([scale, 15.0/amp])
      xr=-[0.0,0.01,0.01,0.0 ,0.0] * amp * scale + xxx + xoff
      polyfill, xr     +0.12, yr  -nn*vert,               /normal, color=160     ;negative
      plots,    xr(3:4)+0.12, [yr(3)+0.002,yr(4)-0.002]-nn*vert, /normal         ;zero
      xr= [-0.01,0.01]             * amp * refer + xxx + xoff
      plots,    xr(0:1)+0.12, yr(0:1)-nn*vert+0.005,      /normal                ;reference 1%

      if nerr eq 1 then begin
      xyouts,   xr(0)     ,   yr(0)-nn*vert, var2d(nn-1), /normal, charsize=0.7
      endif
      xyouts,   xr(0)+0.18,   yr(0)-nn*vert+0.008, err1,  /normal, charsize=0.6  ;new error value
      xyouts,   xr(0)+0.18,   yr(0)-nn*vert-0.008, err2,  /normal, charsize=0.6  ;old error value
    end
  end

  ; surface parameters: ICON and IFS versus observations ------------------------------------

  for nn=1,nvar2obs do begin
   ;print,''
   ;print,'searching for variable ', var2obs(nn-1)
    ind1 = where(var1(0,*) eq var2obs(nn-1) and var1(2-1,*) eq levtype and var1(8-1,*) eq step2 and var1(10-1,*) eq nfor and var1(19-1,*) eq 'diff_obs')
    ind2 = where(var2(0,*) eq var2obs(nn-1) and var2(2-1,*) eq levtype and var2(8-1,*) eq step2 and var2(10-1,*) eq nfor and var2(19-1,*) eq 'diff_obs')
    ind1=ind1(n_elements(ind1)-1)
    ind2=ind2(n_elements(ind2)-1)
    if ind1 eq -1 or ind2 eq -1 then continue
    mean1= float(var1(16-1,ind1))
    mean2= float(var2(16-1,ind2))
    rms1 = float(var1(18-1,ind1))
    rms2 = float(var2(18-1,ind2))
    if stat eq 'rms' then begin
      bias = abs(mean1) / ( abs(mean1) + abs(mean2) ) ; [0,1], 0: mean1 is much better, 1: mean2 is much better: 0.5 same
    endif else begin
      bias = ( abs(mean1) - abs(mean2) ) / (     rms1   +     rms2   ) + 0.5
    endelse
    rms  =     rms1   / (     rms1   +     rms2   )
   ;print, var2obs(nn-1), ' exp',expref, ' bias:', bias2, 'rms:', rms2 
  
    xxx = 0.55
    yyy = 0.30
    vert= 0.038
  
    for nerr=1,2 do begin
      if nerr eq 2 and stat eq 'rms' then continue
      if nerr eq 1 then begin
        err1=mean1
        err2=mean2
        err =bias
        xoff=0.0
        amp =amp_sfc_bias
      endif else begin
        err1=rms1
        err2=rms2
        err =rms
        xoff=0.18
        amp =amp_sfc_rms
      endelse
      scale = (err-0.5) * (err gt 0.5)
      scale = min([scale, 15.0/amp])
      xr= [0.0,0.01,0.01,0.0 ,0.0] * amp * scale + xxx + xoff
      yr= [0.0,0.0 ,0.01,0.01,0.0] * 1.0            + yyy
      polyfill, xr   +0.12 , yr   -nn*vert,               /normal, color=240     ;positive
      scale = (0.5-err) * (err lt 0.5)
      scale = min([scale, 15.0/amp])
      xr=-[0.0,0.01,0.01,0.0 ,0.0] * amp * scale + xxx + xoff
      polyfill, xr     +0.12, yr  -nn*vert,               /normal, color=160     ;negative
      plots,    xr(3:4)+0.12, [yr(3)+0.002,yr(4)-0.002]-nn*vert, /normal         ;zero
      xr= [-0.01,0.01]             * amp * refer + xxx + xoff
      plots,    xr(0:1)+0.12, yr(0:1)-nn*vert+0.005,      /normal                ;reference 1%

      if nerr eq 1 then begin
      xyouts,   xr(0)     ,   yr(0)-nn*vert, var2obs(nn-1) + ' - ' + obsname(nn-1), /normal, charsize=0.7
      endif
      xyouts,   xr(0)+0.18,   yr(0)-nn*vert+0.008, err1,  /normal, charsize=0.6
      xyouts,   xr(0)+0.18,   yr(0)-nn*vert-0.008, err2,  /normal, charsize=0.6
    end
  end


  ; land parameters (T and moisture) --------------------------------------------------------

  for nn=1,nvar3lnd do begin
   for nk=1,nlevlnd do begin
   ;print,''
   ;print,'searching for variable ', var3lnd(nn-1)
    ind1 = where(var1(0,*) eq var3lnd(nn-1) and var1(2-1,*) eq levtype and var1(8-1,*) eq step2 and var1(10-1,*) eq nfor and var1(4-1,*) eq nk and var1(19-1,*) eq diff_txt)
    ind2 = where(var2(0,*) eq var3lnd(nn-1) and var2(2-1,*) eq levtype and var2(8-1,*) eq step2 and var2(10-1,*) eq nfor and var2(4-1,*) eq nk and var2(19-1,*) eq diff_txt)
    ind1=ind1(n_elements(ind1)-1)
    ind2=ind2(n_elements(ind2)-1)
    if ind1 eq -1 or ind2 eq -1 then continue
    mean1= float(var1(16-1,ind1))
    mean2= float(var2(16-1,ind2))
    rms1 = float(var1(18-1,ind1))
    rms2 = float(var2(18-1,ind2))
    if stat eq 'rms' then begin
      bias = abs(mean1) / ( abs(mean1) + abs(mean2) ) ; [0,1], 0: mean1 is much better, 1: mean2 is much better: 0.5 same
    endif else begin
      bias = ( abs(mean1) - abs(mean2) ) / (     rms1   +     rms2   ) + 0.5
    endelse
    rms  =     rms1   / (     rms1   +     rms2   )
   ;print, var3lnd(nn-1), ' exp',expref, ' bias:', bias2, 'rms:', rms2 
  
    xxx = 0.55
    yyy = 0.92
    vert= 0.038
  
    if stat eq 'rms' then begin
      header = 'variable                    rms                new,ref'
    endif else begin
      header = 'variable                    bias               new,ref                     rms               new,ref'
    endelse
    xyouts, xxx, yyy, header, /normal, charsize=0.8

    for nerr=1,2 do begin
      if nerr eq 2 and stat eq 'rms' then continue
      if nerr eq 1 then begin
        err1=mean1
        err2=mean2
        err =bias
        xoff=0.0
        amp =amp_sfc_bias
      endif else begin
        err1=rms1
        err2=rms2
        err =rms
        xoff=0.18
        amp =amp_sfc_rms
      endelse
      vertical = ( (nn-1)*nlevlnd + nk ) * vert
      scale = (err-0.5) * (err gt 0.5)
      scale = min([scale, 15.0/amp])
      xr= [0.0,0.01,0.01,0.0 ,0.0] * amp * scale + xxx + xoff
      yr= [0.0,0.0 ,0.01,0.01,0.0] * 1.0            + yyy
      polyfill, xr   +0.12 , yr   -vertical,               /normal, color=240     ;positive
      scale = (0.5-err) * (err lt 0.5)
      scale = min([scale, 15.0/amp])
      xr=-[0.0,0.01,0.01,0.0 ,0.0] * amp * scale + xxx + xoff
      polyfill, xr     +0.12, yr  -vertical,               /normal, color=160     ;negative
      plots,    xr(3:4)+0.12, [yr(3)+0.002,yr(4)-0.002]-vertical, /normal         ;zero
      xr= [-0.01,0.01]             * amp * refer + xxx + xoff
      plots,    xr(0:1)+0.12, yr(0:1)-vertical+0.005,      /normal                ;reference 1%

      if nerr eq 1 then begin
      xyouts,   xr(0)     ,   yr(0)-vertical, var3lnd(nn-1)+'  L'+strcompress(string(nk)), /normal, charsize=0.7
      endif
      xyouts,   xr(0)+0.18,   yr(0)-vertical+0.008, err1,  /normal, charsize=0.6
      xyouts,   xr(0)+0.18,   yr(0)-vertical-0.008, err2,  /normal, charsize=0.6
    end
   end  
  end

  if stat eq 'rms' then begin
    xr = [-0.01,0.0] * amp_sfc_bias * refer
    plots, xr(0:1)+0.40, 0.015+[0.0,0.0],  /normal                                ;reference 1% (bottom page)
    xyouts,xr(1)  +0.41, 0.010, refer_txt+' rms', /normal, charsize=0.6           ;text refernce value
  
  endif else begin

    xr = [-0.01,0.0] * amp_sfc_bias * refer
    plots, xr(0:1)+0.40, 0.015+[0.0,0.0],  /normal                                ;reference 1% (bottom page)
    xyouts,xr(1)  +0.41, 0.010, refer_txt+' bias', /normal, charsize=0.6          ;text refernce value
  
    xr = [-0.01,0.0] * amp_sfc_rms * refer
    plots, xr(0:1)+0.60, 0.015+[0.0,0.0],  /normal                                ;reference 1% (bottom page)
    xyouts,xr(1)  +0.61, 0.010, refer_txt+' rms', /normal, charsize=0.6           ;text refernce value
  endelse

 END

ENDCASE


paperclose


end
