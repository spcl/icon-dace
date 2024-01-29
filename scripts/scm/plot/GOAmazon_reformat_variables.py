# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

class result_container:
    pass

def varmeta(iname, indata1=0, indata2=0, indata3=0, indata4=0, indata5=0, indata6=0, indata7=0):
  ''' define variables for GOAmazon project'''

  print('processing ICON variable name:', iname)

  var = result_container()

  lmetaonly = isinstance(indata1, int)
  print('metadata only?                ', lmetaonly)

  if not lmetaonly: var.data = indata1[:]
 
  '''
   ml_varlist       = 'pres_sfc','accshfl_s','acclhfl_s','accthb_s','accthb_t','accsob_s','accsob_t', 
                      'accthd_s', 'accthu_s', , 'accsod_t', 'accsou_t', 'accsou_t',
                      'accumfl_s','accvmfl_s','cape','rain_gsp','rain_con','snow_gsp','snow_con', 
                      'u_10m','v_10m','t_2m','td_2m','t_g','qv_s',
  
                      'clct','clch','clcm','clcl','tqv_dia','tqc_dia','tqi_dia','tqr','tqs'
   pl_varlist       = 'u','v','w','temp','rho','geopot','rh','clc','tot_qv_dia','tot_qc_dia','tot_qi_dia',
  !                   'ddt_temp_dyn', 'ddt_temp_radsw', 'ddt_temp_radlw',  'ddt_temp_turb', 'ddt_temp_drag', 
  !                   'ddt_temp_pconv','ddt_qv_turb', 'ddt_qc_turb', 'ddt_qi_turb', 'ddt_qv_conv',
  !                   'ddt_temp_gscp','ddt_qv_gscp','ddt_qc_gscp','ddt_qi_gscp'
  '''

# PL variables

#             ICON name                                              GOAmazon variable meta data
  if iname == 'u'         : var.out_name, var.units, var.long_name = 'ua' , 'm/s', 'Zonal wind'
  if iname == 'v'         : var.out_name, var.units, var.long_name = 'va' , 'm/s', 'Meridional wind'
  if iname == 'temp'      : var.out_name, var.units, var.long_name = 'ta' , 'K', 'Temperature'
  if iname == 'tot_qv_dia': var.out_name, var.units, var.long_name = 'hus', 'kg/kg', 'Specific humidity'
  if iname == 'rh'        : var.out_name, var.units, var.long_name = 'hur', 'percent', 'Relative humidity'
  if iname == 'geopot'    : var.out_name, var.units, var.long_name = 'zg' , 'm', 'Geopotential Height (above sea level)'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'wap', 'Pa/s', 'Vertical velocity (pressure)'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'cls', 'fraction', 'Cloud fraction in CLUBB'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'clc', 'fraction', 'Convective cloud cover'
  if iname == 'clc'       : var.out_name, var.units, var.long_name = 'cl' , '%', 'Cloud fraction'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'clws' , 'kg/kg', 'Stratiform cloud liquid water'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'clwc' , 'kg/kg', 'Convective cloud liquid water'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'clis' , 'kg/kg', 'Stratiform cloud ice water'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'clic' , 'kg/kg', 'Convective cloud ice water'
  if iname == 'tot_qc_dia': var.out_name, var.units, var.long_name = 'clw'  , 'kg/kg', 'Grid box averaged cloud liquid amount'
  if iname == 'tot_qi_dia': var.out_name, var.units, var.long_name = 'cli'  , 'kg/kg', 'Grid box averaged cloud ice amount'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'pfull', 'Pa', 'pressure'
  if iname == 'rho'       : var.out_name, var.units, var.long_name = 'rho'  , 'kg/m3', 'Air Density'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'sclwc', 'g/kg', 'Cloud Water Mixing Ratio'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'dclwc', 'kg/kg', 'Mass fraction of deep convective cloud liquid water'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'sclic', 'kg/kg', 'Mass fraction of shallow convective cloud ice water'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'dclic', 'kg/kg', 'Mass fraction of deep convective cloud ice water'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'mcu'  , 'kg/m2/s', 'ZM convection updraft mass flux'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'mcd'  , 'kg/m2/s', 'ZM convection downdraft mass flux'
  if iname == ('ddt_temp_radsw', 'ddt_temp_radlw', 'ddt_temp_turb', 'ddt_temp_drag', 'ddt_temp_pconv', 'ddt_temp_gscp'): 
                            var.out_name, var.units, var.long_name = 'q1', 'K/s', 'T total physics tendency'
                            if not lmetaonly: var.data = var.data + indata2[:] + indata3[:] + indata4[:] + indata5[:] + indata6[:]
  if iname == ('ddt_qv_turb', 'ddt_qv_conv', 'ddt_qv_gscp'):
                            var.out_name, var.units, var.long_name = 'q2', 'kg/kg/s', 'Q total physics tendency'
                            if not lmetaonly: var.data = var.data + indata2[:] + indata3[:]
  if iname == 'ddt_temp_radsw': var.out_name, var.units, var.long_name = 'tntrsw', 'K/s', 'Solar heating rate'
  if iname == 'ddt_temp_radlw': var.out_name, var.units, var.long_name = 'tntrlw', 'K/s', 'Longwave heating rate'
  if iname == ('ddt_temp_radsw', 'ddt_temp_radlw', 'ddt_temp_turb', 'ddt_temp_drag', 'ddt_temp_pconv', 'ddt_temp_gscp','ddt_temp_dyn'): 
                            var.out_name, var.units, var.long_name = 'tnt', 'K/s', 'Total change rate of temperature'
                            if not lmetaonly: var.data = var.data + indata2[:] + indata3[:] + indata4[:] + indata5[:] + indata6[:] + indata7[:]
  if iname == 'ddt_temp_dyn': var.out_name, var.units, var.long_name = 'tnta', 'K/s', 'change rate of temperature due to total advection'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'tntdp', 'K/s' 'T tendency - Zhang-McFarlane moist convection'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'tntsh', 'K/s', 'T tendency - shallow convection'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'tntpbl', 'K/s', 'change rate of temperature due to CLUBB scheme'
  if iname == 'ddt_temp_gscp': var.out_name, var.units, var.long_name = 'tntlscp', 'K/s', 'change rate of temperature due to microphysics'
  if iname ==  ('ddt_temp_radsw', 'ddt_temp_radlw', 'ddt_temp_turb', 'ddt_temp_drag', 'ddt_temp_pconv'): 
                            var.out_name, var.units, var.long_name = 'tntd', 'K/s', 'change rate of temperature due to other diabatic processes'
                            if not lmetaonly: var.data = var.data + indata2[:] + indata3[:] + indata4[:] + indata5[:]
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'tnhus'    , 'kg/kg/s', 'Total change rate of specific humidity'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'tnhusa'   , 'kg/kg/s', 'change rate of specific humidity due to total advection'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'tnhusdp'  , 'kg/kg/s', 'Q tendency - Zhang-McFarlane moist convection'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'tnhussh'  , 'kg/kg/s', 'QV tendency - shallow convection'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'tnhuspbl' , 'kg/kg/s', 'change rate of specific humidity due to CLUBB'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'tnhuslscp', 'kg/kg/s', 'Q tendency - Morrison microphysics'
  if iname == ('ddt_qv_turb', 'ddt_qv_conv'): 
                            var.out_name, var.units, var.long_name = 'tnhusd'   , 'kg/kg/s', 'change rate of specific humidity due to other diabatic processes'      
                            if not lmetaonly: var.data = var.data + indata2[:]

# ML variables

  if iname == 'tot_prec'  : 
                            var.out_name, var.units, var.long_name = 'pr', 'mm/day', 'Total Precipitation'
  if iname == ('rain_con','snow_gsp'): 
                            var.out_name, var.units, var.long_name = 'prc', 'mm/day', 'Convective Precipitation'
                            if not lmetaonly: var.data = var.data + indata2[:]
  if iname == 'tqv_dia'   : var.out_name, var.units, var.long_name = 'prw', 'kg/m2', 'Total (vertically integrated) precipitable water'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'tas', 'K', 'Reference height temperature'
  if iname == 't_g'       : var.out_name, var.units, var.long_name = 'ts', 'K', 'Surface temperature (radiative)'
  if iname == 'accsod_t'  : var.out_name, var.units, var.long_name = 'rsdt',   'W/m2', 'TOA insident shortwave radiation'
  if iname == 'accsou_t'  : var.out_name, var.units, var.long_name = 'rsut',   'W/m2', 'Upwelling solar flux at top of atmosphere'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'rsutcs', 'W/m2', 'Clearsky upwelling solar flux at top of atmosphere'
  if iname == 'accthb_t'  : var.out_name, var.units, var.long_name = 'rlut',   'W/m2', 'Upwelling longwave flux at top of model'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'rlutcs', 'W/m2', 'Clearsky upwelling longwave flux at top of model'
  if iname == ('accsob_s','accsod_s'): 
                            var.out_name, var.units, var.long_name = 'rsus',   'W/m2', 'Surface upward shortwave radiation'
                            if not lmetaonly: var.data = var.data - indata2[:]
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'rsuscs', 'W/m2', 'Surface upward shortwave radiation in clearsky'
  if iname == 'accsod_s'  : var.out_name, var.units, var.long_name = 'rsds',   'W/m2', 'Downwelling solar flux at surface'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'rsdscs', 'W/m2', 'Clearsky downwelling solar flux at surface'
  if iname == 'accthu_s'  : var.out_name, var.units, var.long_name = 'rlus',   'W/m2', 'Surface upward longwave radiation'
  if iname == 'accthd_s'  : var.out_name, var.units, var.long_name = 'rlds',   'W/m2', 'Downwelling longwave flux at surface'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'rldscs', 'W/m2', 'Surface downward longwave radiation in clearsky'
  if iname == 'accshfl_s' : var.out_name, var.units, var.long_name = 'hfss',   'W/m2', 'Surface sensible heat flux'
  if iname == 'acclhfl_s' : var.out_name, var.units, var.long_name = 'hfls',   'W/m2', 'Surface latent heat flux'
  if iname == 'pres_sfc'  : var.out_name, var.units, var.long_name = 'ps', 'Pa', 'Surface pressure'
  if iname == 'u_10m'     : var.out_name, var.units, var.long_name = 'uas', 'm/s', 'Surface u-wind'
  if iname == 'v_10m'     : var.out_name, var.units, var.long_name = 'vas', 'm/s', 'Surface v-wind'
  if iname == 'tqc_dia'   : var.out_name, var.units, var.long_name = 'cllvi', 'kg/m2', 'Total grid-box cloud liquid water path'
  if iname == 'tqi_dia'   : var.out_name, var.units, var.long_name = 'clivi', 'kg/m2', 'Total grid-box cloud ice water path'
  if iname == ('tqc_dia','tqi_dia'): 
                            var.out_name, var.units, var.long_name = 'clwvi', 'kg/m2', 'Total grid-box cloud water path (liquid and ice)'
                            if not lmetaonly: var.data = var.data + indata2[:]
  if iname == ('tqr','tqs'): 
                            var.out_name, var.units, var.long_name = 'cliviall', 'kg/m2', 'Ice + Snow water path'
                            if not lmetaonly: var.data = var.data + indata2[:]
  if iname == 'clct'      : 
                            var.out_name, var.units, var.long_name = 'clt', '%', 'Vertically-integrated total cloud'
                            if not lmetaonly: var.data = var.data
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'qs', 'kg/kg', 'Reference height humidity'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'pblh', 'm', 'PBL height'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'mrsos', 'kg/m2', '0-10cm Soil Moisture'
  if iname == 'tttt'      : var.out_name, var.units, var.long_name = 'orog', 'm', 'Surface altitude'


  if (iname == 'accsod_t' or iname == 'accsou_t' or iname == 'accthb_t' or iname == ('accsob_s','accsod_s') or iname == 'accsod_s' or
      iname == 'accsod_s' or iname == 'accthu_s' or iname == 'accthd_s' or iname == 'accshfl_s'             or iname == 'acclhfl_s') : 
    if not lmetaonly: first            = (var.data[1]               - var.data[0])               / 1800.00 # first
    if not lmetaonly: last             = (var.data[var.data.size-1] - var.data[var.data.size-2]) / 1800.00 # last
    if not lmetaonly: var.data[1:var.data.size-1] = (var.data[2:] - var.data[0:var.data.size-2]) / 3600.00
    if not lmetaonly: var.data[0]                 = first
    if not lmetaonly: var.data[var.data.size-1]   = last

  if (iname == 'tot_prec' or iname == ('rain_con','snow_gsp') ) :
    if not lmetaonly: first            = (var.data[1]               - var.data[0])               * 48 # first
    if not lmetaonly: last             = (var.data[var.data.size-1] - var.data[var.data.size-2]) * 48 # last
    if not lmetaonly: var.data[1:var.data.size-1] = (var.data[2:] - var.data[0:var.data.size-2]) * 24
    if not lmetaonly: var.data[0]                 = first
    if not lmetaonly: var.data[var.data.size-1]   = last

  if (iname == ('accsob_s','accsod_s')) or iname == 'accthb_t' or iname == 'accshfl_s' or iname == 'acclhfl_s' :
    
    if not lmetaonly: var.data = - var.data             

  return var
