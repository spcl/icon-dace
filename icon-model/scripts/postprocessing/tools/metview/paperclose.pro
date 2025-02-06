; ICON
;
; ---------------------------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ---------------------------------------------------------------

pro paperclose

 device, /close_file, /helvetica
;spawn,'lpr plot.ps'
;spawn,'rm -f plot.ps'
 set_plot,'x'

 !p.charsize=1 & !p.charthick=1 & !p.font=-1
 !p.thick=1 & !x.thick=1 & !y.thick=1
 !p.position=0
;!x.margin = 0 & !y.margin = 0
 !p.region = 0

 end
