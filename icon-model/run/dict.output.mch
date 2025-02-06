# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

# Dictionary for mapping between internal names and GRIB2 shortNames
# needed by GRIB2 write procedures.
#
# GRIB2 shortName   internal name
# --------------------------------
#
  va                 v
  wa                 w
  wap                omega
  ta                 temp
  ps                 pres_sfc
  pmsl               pres_msl
  hus                qv
  pfull              pres
  phalf              pres_ifc
   THETA_V        theta_v          
   DEN            rho              
   DTKE_CON       ddt_tke_pconv    
   FI             geopot           
   HHL            z_ifc            
   VN             vn               
   U              u                
   V              v                
   W              w                
   TKE            tke              
   T              temp             
   P              pres             
   PMSL           pres_msl         
   PS             pres_sfc         
   QV             qv               
   QC             qc               
   QI             qi               
   QR             qr               
   QS             qs               
   QG             qg               
   QH             qh               
   NCCLOUD        qnc              
   NCRAIN         qnr              
   NCICE          qni              
   NCSNOW         qns              
   NCGRAUPEL      qng              
   NCHAIL         qnh              
   T_G            t_g              
   QV_S           qv_s             
   FR_ICE         fr_seaice        
   T_ICE          t_ice            
   H_ICE          h_ice            
   T_SNOW         t_snow           
   FRESHSNW       freshsnow        
   SNOWC          snowfrac_lc      
   W_SNOW         w_snow           
   RHO_SNOW       rho_snow         
   H_SNOW         h_snow           
   HSNOW_MAX      hsnow_max        
   SNOAG          snow_age         
   T_SNOW_M       t_snow_mult      
   RHO_SNOW_M     rho_snow_mult    
   W_SNOW_M       wtot_snow        
   WLIQ_SNOW_M    wliq_snow        
   H_SNOW_M       dzh_snow         
   W_I            w_i              
   W_SO           w_so             
   W_SO_ICE       w_so_ice         
   SMI            smi              
   T_SO           t_so             
   SKT            t_sk             
   T_SEA          t_seasfc         
   Z0             gz0              
   T_MNW_LK       t_mnw_lk         
   T_WML_LK       t_wml_lk         
   H_ML_LK        h_ml_lk          
   T_BOT_LK       t_bot_lk         
   C_T_LK         c_t_lk           
   T_B1_LK        t_b1_lk          
   H_B1_LK        h_b1_lk          
   RELHUM         rh               
   RELHUM_2M      rh_2m            
   RELHUM_2M_L    rh_2m_land       
   TD_2M_L        td_2m_land       
   T_2M_L         t_2m_land        
   T_2M           t_2m             
   T_2M_FILTBIAS  t2m_bias         
   ALB_RAD        albdif           
   ALB_SEAICE     alb_si           
   ASWDIFD_S      asodifd_s        
   ASWDIFU_S      asodifu_s        
   ASWDIR_S       asodird_s        
   HSURF          topography_c     
   VMAX_10M       gust10           
   AER_SS         aer_ss           
   AER_ORG        aer_or           
   AER_BC         aer_bc           
   AER_SO4        aer_su           
   AER_DUST       aer_du           
   ALB_SEAICE     alb_si           
   EVAP_PL        plantevap        
