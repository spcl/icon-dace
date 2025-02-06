! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Contains the setup of variables related to large eddy simulation setup

MODULE mo_les_nml

  USE mo_les_config,          ONLY: les_config
  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_les_namelist
  PUBLIC :: turb_profile_list, turb_tseries_list

  CHARACTER(LEN=7) :: turb_tseries_list(19), turb_profile_list(52) !list of variables
                                !added profiles of LS forcing tendencies (45-51)

CONTAINS
  !-------------------------------------------------------------------------
  !! Read Namelist for LES
  !!
  !! This subroutine 
  !! - reads the Namelist 
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state 
  !!
  SUBROUTINE read_les_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename 
    INTEGER :: istat, funit, jg
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_les_nml: read_les_namelist'

    REAL(wp) :: sst        ! prescribed SST
    REAL(wp) :: psfc       ! prescribed surface pressure
    REAL(wp) :: shflx      ! prescribed sensible heat flux (Km/s)
    REAL(wp) :: lhflx      ! prescribed latent heat flux   (Km/s)
    INTEGER  :: isrfc_type ! 1=fixed sst, 2=fixed flux, 3=fixed buyancy flux
    INTEGER  :: smag_coeff_type  ! 1=Smagorinsky model; 2=set coeff. externally by Km_ext, Kh_ext (for tests)

    REAL(wp) :: ufric      ! friction velocity

    LOGICAL  :: is_dry_cbl  !special case for CBL testcase

    !For isrf_type==3
    REAL(wp) :: bflux      !Buoyancy flux
    REAL(wp) :: tran_coeff !Surface transfer coefficient in units of velocity (m/s)

    !Some parameters
    REAL(wp) :: smag_constant
    REAL(wp) :: turb_prandtl
    REAL(wp) :: km_min         !min mass weighted turbulent viscosity
    REAL(wp) :: Km_ext         !externally set constant kinematic viscosity (m2/s)
    REAL(wp) :: Kh_ext         !externally set constant diffusion coeff.    (m2/s)
    REAL(wp) :: max_turb_scale !max turbulence length scale
    REAL(wp) :: min_sfc_wind  !min sfc wind in free convection limit

    !Scheme for vertical discretization
    INTEGER :: vert_scheme_type !1=explicit, 2=implicit

    !Parameters for additional diagnostic output
    LOGICAL  :: ldiag_les_out                    !.TRUE. to turn it on
    REAL(wp) :: avg_interval_sec, sampl_freq_sec !averaging and sampling time
    CHARACTER(MAX_CHAR_LENGTH) :: expname        !name of experiment for naming the file
    LOGICAL  :: les_metric

    NAMELIST/les_nml/ sst, psfc, shflx, lhflx, isrfc_type, ufric, is_dry_cbl, &
         smag_constant, turb_prandtl, bflux, tran_coeff,   &
         vert_scheme_type, avg_interval_sec, sampl_freq_sec,  &
         expname, ldiag_les_out, km_min, min_sfc_wind, les_metric, &
         max_turb_scale, Km_ext, Kh_ext, smag_coeff_type

    !-----------------------
    ! 1. default settings
    !-----------------------
    sst          = 300._wp
    psfc         = -999._wp
    shflx        = 0.1_wp 
    lhflx        = 0._wp 
    isrfc_type   = 1 
    ufric        = -999._wp 

    is_dry_cbl   = .FALSE.

    !parameters
    smag_constant    = 0.23_wp
    turb_prandtl     = 0.33333333333_wp
    km_min           = 0.001_wp
    Km_ext           = 75.0_wp
    Kh_ext           = 75.0_wp
    smag_coeff_type  = 1
    max_turb_scale   = 300._wp
    min_sfc_wind     = 1._wp !Default from Holstag and Boville 1991

    bflux       = 0.0007_wp
    tran_coeff  = 0.02_wp

    vert_scheme_type = 2 !implicit

    !output parameters
    ldiag_les_out = .FALSE. 
    expname  = 'ICOLES'
    avg_interval_sec = 900._wp
    sampl_freq_sec   = 60._wp

    turb_profile_list = (/                                                     &
      'u      ','v      ','w      ','th     ','exner  ','rho    ','qv     ',   & !1-7
      'qc     ','wu     ','wv     ','wth    ','wqv    ','wqc    ','ww     ',   & !8-14
      'thth   ','qvqv   ','qcqc   ','uu     ','vv     ','kh     ','km     ',   & !15-21
      'thv    ','wthv   ','wqvd   ','wthd   ','wqcd   ','bruvais','mechprd',   & !22-28
      'wud    ','wvd    ','wthsfs ','rh     ','clc    ','qi     ','qs     ',   & !29-35
      'qr     ','qg     ','qh     ','lwf    ','swf    ','dt_t_sw','dt_t_lw',   & !36-42
      'dt_t_tb','dt_t_mc','dthls_w','dqls_w ','dthls_h','dqls_h ','nt_thl ',   & !43-49 
      'nt_qt  ','wfls   ','tke    ' /)                                           !50-52

    turb_tseries_list = (/                                          &
      'ccover ','shflx  ','lhflx  ','ustress','vstress','tsfc   ',  & !1-6
      'qsfc   ','hbl    ','psfc   ','swf_tom','lwf_tom','swf_sfc',  & !7-12
      'lwf_sfc','precp_t','precp_r','precp_s','precp_g','precp_h',  & !13-18
      'precp_i' /)                                                    !19

    ! grid metric terms in the les diffusion
    les_metric = .FALSE.


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('les_nml')
      READ(funit,NML=les_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('les_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, les_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, les_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, les_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    DO jg = 1 , max_dom
      les_config(jg)% sst          =  sst
      les_config(jg)% psfc         =  psfc
      les_config(jg)% shflx        =  shflx
      les_config(jg)% lhflx        =  lhflx
      les_config(jg)% isrfc_type   =  isrfc_type
      les_config(jg)% ufric        =  ufric
      les_config(jg)% is_dry_cbl   =  is_dry_cbl
      les_config(jg)% smag_constant     =  smag_constant
      les_config(jg)% turb_prandtl      =  turb_prandtl
      les_config(jg)% rturb_prandtl     =  1._wp/turb_prandtl
      les_config(jg)% bflux             =  bflux
      les_config(jg)% tran_coeff        =  tran_coeff
      les_config(jg)% vert_scheme_type  =  vert_scheme_type
      les_config(jg)% ldiag_les_out     =  ldiag_les_out
      les_config(jg)% expname           =  expname
      les_config(jg)% avg_interval_sec  =  avg_interval_sec
      les_config(jg)% sampl_freq_sec    =  sampl_freq_sec
      les_config(jg)% km_min            =  km_min
      les_config(jg)% Km_ext            =  Km_ext
      les_config(jg)% Kh_ext            =  Kh_ext
      les_config(jg)% smag_coeff_type   =  smag_coeff_type
      les_config(jg)% max_turb_scale    =  max_turb_scale
      les_config(jg)% min_sfc_wind      =  min_sfc_wind
      les_config(jg)% les_metric        =  les_metric
    END DO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=les_nml)                    
      CALL store_and_close_namelist(funit,'les_nml') 
    ENDIF
    ! 7. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=les_nml)

  END SUBROUTINE read_les_namelist

END MODULE mo_les_nml
