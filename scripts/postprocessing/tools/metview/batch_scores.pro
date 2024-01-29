; ICON
;
; ---------------------------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ---------------------------------------------------------------

; batch execution of ICON scores plotting
;-----------------------------------------------------------------------
; Useage: run as: idl batch_scores -args [direxp] [dirref] [expnum]
; [expref] [inidate] [step] [nfor]
;
; e.g.:
; idl batch_scores -args /scratch/ms/de/deia/icon/dei2/252/201508/metplots \
;                        /scratch/ms/de/deia/icon/dei2/247/201508/metplots dei2_252 dei2_247 20150801 24 31

; Martin Koehler, Sep 2012 & Dec 2012
;-----------------------------------------------------------------------

 args  = COMMAND_LINE_ARGS(count=nargs)
 direxp  = args(0)          ;long(args(0)), fix(args(2))
 dirref  = args(1)
 expnum  = args(2)
 expref  = args(3)
 inidate = args(4)
 step    = args(5)
 nfor    = args(6)

 print,'Input arguments: direxp=    ',direxp
 print,'                 dirref=    ',dirref
 print,'                 expnum=    ',expnum
 print,'                 expref=    ',expref
 print,'                 inidate=   ',inidate
 print,'                 step=      ',step
 print,'                 nfor=      ',nfor

 .run scores

 outdate   = long(inidate) + step/24
 plotfile2 = strcompress(string(outdate),/remove_all)+'00.L1.ps'
 
; mean plots

 plotfile1 = direxp+'/scores.NWP.r2B06L90.'+expnum+'_vs_'+expref

 scores, direxp, dirref, expnum, expref, inidate, step, step, nfor, 'ml', 'mean'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.ml.'  + plotfile2
 
 scores, direxp, dirref, expnum, expref, inidate, step, step, nfor, 'pl', 'mean'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.pl.'  + plotfile2

 scores, direxp, dirref, expnum, expref, inidate, step, step, nfor, 'zl', 'mean'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.zl.'  + plotfile2

 scores, direxp, dirref, expnum, expref, inidate, step, step, nfor, 'sfc', 'mean'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.sfc.' + plotfile2

; RMS plots

 plotfile1 = direxp+'/scores.NWP.r2B06L90.'+expnum+'_rms_vs_'+expref
 
 scores, direxp, dirref, expnum, expref, inidate, step, step, nfor, 'ml', 'rms'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.ml.'  + plotfile2
 
 scores, direxp, dirref, expnum, expref, inidate, step, step, nfor, 'pl', 'rms'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.pl.'  + plotfile2

 scores, direxp, dirref, expnum, expref, inidate, step, step, nfor, 'zl', 'rms'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.zl.'  + plotfile2

 scores, direxp, dirref, expnum, expref, inidate, step, step, nfor, 'sfc', 'rms'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.sfc.' + plotfile2

; RMS time-mean plots

 step1=24
 stepx=24

 plotfile1 = direxp+'/scores.NWP.r2B06L90.'+expnum+'_rms_vs_'+expref
 outdate   = long(inidate) + step1/24
 plotfile2 = strcompress(string(outdate),/remove_all) + "-" + strmid(string(101+step/24,format='(I3)'),1,2) + '_00.L1.ps'
 
 scores, direxp, dirref, expnum, expref, inidate, step1, step, nfor, 'ml', 'rms'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.ml.'  + plotfile2
 
 scores, direxp, dirref, expnum, expref, inidate, step1, step, nfor, 'pl', 'rms'
 spawn,'\cp -f plot.ps ' + plotfile1 + '.pl.'  + plotfile2


 exit

