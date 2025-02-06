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

MODULE mo_io_config
  USE mo_cdi,                     ONLY: FILETYPE_NC2
  USE mo_exception,               ONLY: finish, message
  USE mo_impl_constants,          ONLY: max_dom, max_echotop, max_wshear, max_srh
  USE mo_io_units,                ONLY: filename_max
  USE mo_kind,                    ONLY: wp
  USE mo_parallel_config,         ONLY: num_restart_procs
  USE mo_run_config,              ONLY: dtime
  USE mo_util_string,             ONLY: int2string
  USE mtime,                      ONLY: max_timedelta_str_len
  USE mo_name_list_output_config, ONLY: is_variable_in_output, &
    &                                   is_variable_in_output_dom
  USE mo_lnd_nwp_config,          ONLY: groups_smi

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Derived type
  !--------------------------------------------------------------------------

  ! from namelist

  LOGICAL :: lkeep_in_sync              ! if .true., sync stream after each timestep
  REAL(wp):: dt_diag                    ! diagnostic output timestep [seconds]
  REAL(wp):: gust_interval(max_dom)     ! time interval [seconds] over which maximum wind gusts are taken
  REAL(wp):: ff10m_interval(max_dom)    ! time interval [seconds] over which ff10m is averaged
  REAL(wp):: celltracks_interval(max_dom)  ! time interval [seconds] over which extrema of cell track vars are taken
                                           !  (LPI_MAX, UH_MAX, VORW_CTMAX, W_CTMAX, DBZ_CTMAX)
  CHARACTER(len=max_timedelta_str_len) :: precip_interval(max_dom)   ! time interval over which precipitation variables are accumulated
  CHARACTER(len=max_timedelta_str_len) :: totprec_d_interval(max_dom)! time interval over which the special tot_prec_d is accumulated
  CHARACTER(len=max_timedelta_str_len) :: runoff_interval(max_dom)   ! time interval over which runoff variables are accumulated
  CHARACTER(len=max_timedelta_str_len) :: sunshine_interval(max_dom) ! time interval over which sunshine duration is accumulated
  CHARACTER(len=max_timedelta_str_len) :: melt_interval(max_dom)     ! time interval over which snow melt rate is accumulated
  CHARACTER(len=max_timedelta_str_len) :: maxt_interval(max_dom)     ! time interval for tmax_2m, tmin_2m
  REAL(wp):: dt_lpi                     ! calling frequency [seconds] of lpi diagnosis for hourly maximum calculation
  REAL(wp):: dt_celltracks              ! calling frequency [seconds] of celltrack diagnosis for hourly maximum calculation
                                        ! this pertains to the following variables: tcond_max/tcond10_max, uh_max, vorw_ctmax, w_ctmax
  REAL(wp):: dt_radar_dbz               ! calling frequency [seconds] of radar reflectivity diagnosis for hourly maximum calculation
  REAL(wp):: dt_hailcast                ! calling frequency [seconds] of hailcast for dhail_* diagnostics
  REAL(wp):: wdur_min_hailcast          ! minimal updraft persistence [seconds] for hailcast to be activated

  REAL(wp):: dt_checkpoint              ! timestep [seconds] for triggering new restart file

  INTEGER :: inextra_2d                 ! number of extra output fields for debugging
  INTEGER :: inextra_3d                 ! number of extra output fields for debugging

  LOGICAL :: lflux_avg                  ! if .FALSE. the output fluxes are accumulated
                                        ! from the beginning of the run
                                        ! if .TRUE. the output fluxex are average values
                                        ! from the beginning of the run, except of
                                        ! TOT_PREC that would be accumulated

  INTEGER :: itype_dursun               ! if 0 the sunshine duration is counted if >120W/m^2
                                        ! if 1 the sunshine duration is counted only
                                        ! if direct radiation > 200 W/m^2 and relative sunshine duration in % is computed

  INTEGER :: itype_convindices          ! if 1 CAPE_MU/CIN_MU are approximated via the CAPE/CIN of the parcel with maximum equivalent temperature
                                        ! if 2 the full computation is done

  INTEGER :: itype_hzerocl              ! Specifies height of freezing level if T < 0 Celsius in the whole atmospheric column
                                        ! 1: set hzerocl to orography height (default)
                                        ! 2: set hzerocl to -999.0_wp (undef)
                                        ! 3: set hzerocl to extrapolated value below ground (assuming -6.5 K/km)

  INTEGER :: itype_pres_msl             ! Specifies method for computation of mean sea level pressure
  INTEGER :: itype_rh                   ! Specifies method for computation of relative humidity

  CHARACTER(LEN=filename_max) :: &
    &        output_nml_dict,    &      !< maps variable names onto the internal ICON names.
    &        netcdf_dict                !< maps internal variable names onto names in output file (NetCDF only).

  LOGICAL :: linvert_dict               !< inverts columns in output_nml_dict (allows using the same dictionary file as for input)

  LOGICAL :: lnetcdf_flt64_output       !< if .TRUE. floating point valued NetCDF output
                                        !  is written in 64-bit instead of 32-bit accuracy

  INTEGER, PARAMETER :: uh_max_nlayer = 3
  REAL(wp):: uh_max_zmin(uh_max_nlayer), &  !< minimum and maximum height in m MSL for the vertical
             uh_max_zmax(uh_max_nlayer)     !  integration of uh_max (output vars uh_max<_low,_med>)
  LOGICAL :: luh_max_out(max_dom, uh_max_nlayer) = .false.

  TYPE t_echotop_meta
    REAL(wp) :: time_interval           !< time interval [seconds] over which echotops are maximized/minimized
                                        !   (ECHOTOP, ECHOTOPinM) for a domain
    INTEGER  :: nechotop                !< number of actual echotop levels for domain
    REAL(wp) :: dbzthresh(max_echotop)  !< actual echotop levels from namelist
  END TYPE t_echotop_meta

  TYPE(t_echotop_meta) :: echotop_meta(max_dom)

  REAL(wp):: wshear_uv_heights(max_wshear) !< heights AGL for which wind shear components WSHEAR_U, WSHEAR_V are desired
  INTEGER :: n_wshear                   !< actual number of wshear heights given in namelist
  
  REAL(wp):: srh_heights(max_srh)       !< heights AGL for which storm relative helicity calculations are desired
  INTEGER :: n_srh                      !< actual number of srh heights given in namelist
  
  ! Derived type to collect logical variables indicating if optional diagnostics are requested for output
  TYPE t_var_in_output
    LOGICAL :: pres_msl     = .FALSE. !< Flag. TRUE if computation of mean sea level pressure desired
    LOGICAL :: omega        = .FALSE. !< Flag. TRUE if computation of vertical velocity desired
    LOGICAL :: rh           = .FALSE. !< Flag. TRUE if computation of relative humidity desired
    LOGICAL :: pv           = .FALSE. !< Flag. TRUE if computation of potential vorticity desired
    LOGICAL :: sdi2         = .FALSE. !< Flag. TRUE if computation of supercell detection index desired
    LOGICAL :: lpi          = .FALSE. !< Flag. TRUE if computation of lightning potential index desired
    LOGICAL :: lpi_max      = .FALSE. !< Flag. TRUE if computation of max. of lightning potential index desired
    LOGICAL :: lpi_con      = .FALSE. !< Flag. TRUE if computation of convective lightning potential index desired
    LOGICAL :: mlpi_con     = .FALSE. !< Flag. TRUE if computation of modified convective lightning potential index desired
    LOGICAL :: lpi_con_max  = .FALSE. !< Flag. TRUE if computation of maximum convective lightning potential index desired
    LOGICAL :: mlpi_con_max = .FALSE. !< Flag. TRUE if computation of maximum modified convective lightning potential index desired
    LOGICAL :: lfd_con      = .FALSE. !< Flag. TRUE if computation of lighting flash density desired
    LOGICAL :: lfd_con_max  = .FALSE. !< Flag. TRUE if computation of maximum lighting flash density  desired
    LOGICAL :: koi          = .FALSE. !< Flag. TRUE if computation of convection index
    LOGICAL :: ceiling      = .FALSE. !< Flag. TRUE if computation of ceiling height desired
    LOGICAL :: vis          = .FALSE. !< Flag. TRUE if computation of visibility desired
    LOGICAL :: hbas_sc      = .FALSE. !< Flag. TRUE if computation of height of base from shallow convection desired
    LOGICAL :: htop_sc      = .FALSE. !< Flag. TRUE if computation of height of top  from shallow convection desired
    LOGICAL :: inversion_height = .FALSE. !< Flag. TRUE if computation of height of top  from shallow convection desired
    LOGICAL :: twater       = .FALSE. !< Flag. TRUE if computation of total column integrated water desired
    LOGICAL :: q_sedim      = .FALSE. !< Flag. TRUE if computation of specific content of precipitation particles desired
    LOGICAL :: dbz          = .FALSE. !< Flag. TRUE if computation of radar reflectivity is desired
    LOGICAL :: dbz850       = .FALSE. !< Flag. TRUE if computation of radar reflectivity in approx. 850 hPa is desired
    LOGICAL :: dbzcmax      = .FALSE. !< Flag. TRUE if computation of radar reflectivity column maximum is desired
    LOGICAL :: dbzctmax     = .FALSE. !< Flag. TRUE if computation of radar reflectivity column and time maximum is desired
    LOGICAL :: dbzlmx_low   = .FALSE. !< Flag. TRUE if computation of radar reflectivity layer maximum [500,2500] m AGL is desired
    LOGICAL :: echotop      = .FALSE. !< Flag. TRUE if computation of echo tops in hPa of radar reflectivity is desired
    LOGICAL :: echotopinm   = .FALSE. !< Flag. TRUE if computation of echo tops in m MSL of radar reflectivity is desired
    LOGICAL :: smi          = .FALSE. !< Flag. TRUE if computation of soil moisture index desired
    LOGICAL :: tcond_max    = .FALSE. !< Flag. TRUE if computation of total column-integrated condensate desired
    LOGICAL :: tcond10_max  = .FALSE. !< Flag. TRUE if computation of total column-integrated condensate above z(T=-10 degC) desired
    LOGICAL :: uh_max_low   = .FALSE. !< Flag. TRUE if computation of updraft helicity (0 - 3000 m) desired
    LOGICAL :: uh_max_med   = .FALSE. !< Flag. TRUE if computation of updraft helicity (2000 - 5000 m) desired
    LOGICAL :: uh_max       = .FALSE. !< Flag. TRUE if computation of updraft helicity (2000 - 8000 m) desired
    LOGICAL :: vorw_ctmax   = .FALSE. !< Flag. TRUE if computation of maximum rotation amplitude desired
    LOGICAL :: w_ctmax      = .FALSE. !< Flag. TRUE if computation of maximum updraft track desired
    LOGICAL :: dursun       = .FALSE. !< Flag. TRUE if computation of sunshine duration is required
    LOGICAL :: dursun_m     = .FALSE. !< Flag. TRUE if computation of maximum sunshine duration is required
    LOGICAL :: dursun_r     = .FALSE. !< Flag. TRUE if computation of relative sunshine duration is required
    LOGICAL :: res_soilwatb = .FALSE. !< Flag. TRUE if computation of residuum of soil water is desired
    LOGICAL :: snow_melt    = .FALSE. !< Flag. TRUE if computation of snow melt desired
    LOGICAL :: dhail_mx     = .FALSE. !< Flag. TRUE if computation of maximum expected hail diameter desired
    LOGICAL :: dhail_av     = .FALSE. !< Flag. TRUE if computation of average expected hail diameter desired
    LOGICAL :: dhail_sd     = .FALSE. !< Flag. TRUE if computation of standard deviation of hail diameter desired
    LOGICAL :: wshear_u     = .FALSE. !< Flag. TRUE if computation of vertical U wind shear components is desired
    LOGICAL :: wshear_v     = .FALSE. !< Flag. TRUE if computation of vertical V wind shear components is desired
    LOGICAL :: lapserate    = .FALSE. !< Flag. TRUE if computation of T(500hPa) - T(850hPa) is desired
    LOGICAL :: mconv        = .FALSE. !< Flag. TRUE if computation of low level moisture convergence is desired
    LOGICAL :: srh          = .FALSE. !< Flag. TRUE if computation of storm relative helicity (SRH) is desired
    LOGICAL :: tot_pr_max   = .FALSE. !< Flag. TRUE if computation of time max precipitation rate is desired
    LOGICAL :: cloudtop     = .FALSE. !< Flag. TRUE if computation of CLOUDTOP is desired
    LOGICAL :: si           = .FALSE. !< Flag. TRUE if computation of SI is desired
    LOGICAL :: sli          = .FALSE. !< Flag. TRUE if computation of SLI is desired
    LOGICAL :: swiss12      = .FALSE. !< Flag. TRUE if computation of SWISS12 is desired
    LOGICAL :: swiss00      = .FALSE. !< Flag. TRUE if computation of SWISS00 is desired
    LOGICAL :: cape_mu      = .FALSE. !< Flag. TRUE if computation of most unstable CAPE is desired
    LOGICAL :: cin_mu       = .FALSE. !< Flag. TRUE if computation of most unstable convective inhibition MU is desired
    LOGICAL :: hpbl         = .FALSE. !< Flag. TRUE if computation of boundary layer height is desired
    LOGICAL :: aod_550nm    = .FALSE. !< Flag. TRUE if computation of aerosol optical depth at 550 nm is desired
    LOGICAL :: cape_3km     = .FALSE. !< Flag. TRUE if computation of CAPE 3KM is desired
    LOGICAL :: lfc_ml       = .FALSE. !< Flag. TRUE if computation of the Level of Free Convection is desired
    LOGICAL :: lcl_ml       = .FALSE. !< Flag. TRUE if computation of the Lifted Condensation Level is desired
    LOGICAL :: ddt_temp_drag = .FALSE. !< Flag. TRUE if temp-tend. from sso+gravity-wave-drag+Rayleigh-frict. is required
    ! add vars for global mean claclulations
    LOGICAL :: tas_gmean    = .FALSE. !< Flag. TRUE if computation of global mean T2m 
    LOGICAL :: rsdt_gmean   = .FALSE. !< Flag. TRUE if computation of global mean toa downward short wave rad
    LOGICAL :: rsut_gmean   = .FALSE. !< Flag. TRUE if computation of global mean toa upward short wave rad
    LOGICAL :: rlut_gmean   = .FALSE. !< Flag. TRUE if computation of global mean toa upward long wave rad
    LOGICAL :: prec_gmean   = .FALSE. !< Flag. TRUE if computation of global mean precipitation
    LOGICAL :: evap_gmean   = .FALSE. !< Flag. TRUE if computation of global mean evapotranspiration
    LOGICAL :: pme_gmean    = .FALSE. !< Flag. TRUE if computation of global mean P-E
    LOGICAL :: radtop_gmean = .FALSE. !< Flag. TRUE if computation of global mean net radiation
    !
    ! diagnostics for the horizontal wind tendencies in the dynamical core
    LOGICAL :: ddt_vn_dyn  = .FALSE. !< Flag. TRUE if the storage of ddt_vn_dyn is required
    LOGICAL :: ddt_ua_dyn  = .FALSE. !<                              ddt_ua_dyn
    LOGICAL :: ddt_va_dyn  = .FALSE. !<                              ddt_va_dyn
    LOGICAL :: ddt_vn_dmp  = .FALSE. !< Flag. TRUE if the storage of ddt_vn_dmp is required
    LOGICAL :: ddt_ua_dmp  = .FALSE. !<                              ddt_ua_dmp
    LOGICAL :: ddt_va_dmp  = .FALSE. !<                              ddt_va_dmp
    LOGICAL :: ddt_vn_hdf  = .FALSE. !< Flag. TRUE if the storage of ddt_vn_hdf is required
    LOGICAL :: ddt_ua_hdf  = .FALSE. !<                              ddt_ua_hdf
    LOGICAL :: ddt_va_hdf  = .FALSE. !<                              ddt_va_hdf
    LOGICAL :: ddt_vn_adv  = .FALSE. !< Flag. TRUE if the storage of ddt_vn_adv is required
    LOGICAL :: ddt_ua_adv  = .FALSE. !<                              ddt_ua_adv
    LOGICAL :: ddt_va_adv  = .FALSE. !<                              ddt_va_adv
    LOGICAL :: ddt_vn_cor  = .FALSE. !< Flag. TRUE if the storage of ddt_vn_cor is required
    LOGICAL :: ddt_ua_cor  = .FALSE. !<                              ddt_ua_cor
    LOGICAL :: ddt_va_cor  = .FALSE. !<                              ddt_va_cor
    LOGICAL :: ddt_vn_pgr  = .FALSE. !< Flag. TRUE if the storage of ddt_vn_pgr is required
    LOGICAL :: ddt_ua_pgr  = .FALSE. !<                              ddt_ua_pgr
    LOGICAL :: ddt_va_pgr  = .FALSE. !<                              ddt_va_pgr
    LOGICAL :: ddt_vn_phd  = .FALSE. !< Flag. TRUE if the storage of ddt_vn_phd is required
    LOGICAL :: ddt_ua_phd  = .FALSE. !<                              ddt_ua_phd
    LOGICAL :: ddt_va_phd  = .FALSE. !<                              ddt_va_phd
    LOGICAL :: ddt_vn_iau  = .FALSE. !< Flag. TRUE if the storage of ddt_vn_iau is required
    LOGICAL :: ddt_ua_iau  = .FALSE. !<                              ddt_ua_iau
    LOGICAL :: ddt_va_iau  = .FALSE. !<                              ddt_va_iau
    LOGICAL :: ddt_vn_ray  = .FALSE. !< Flag. TRUE if the storage of ddt_vn_ray is required
    LOGICAL :: ddt_ua_ray  = .FALSE. !<                              ddt_ua_ray
    LOGICAL :: ddt_va_ray  = .FALSE. !<                              ddt_va_ray
    LOGICAL :: ddt_vn_grf  = .FALSE. !< Flag. TRUE if the storage of ddt_vn_grf is required
    LOGICAL :: ddt_ua_grf  = .FALSE. !<                              ddt_ua_grf
    LOGICAL :: ddt_va_grf  = .FALSE. !<                              ddt_va_grf
  END TYPE t_var_in_output

  TYPE(t_var_in_output), ALLOCATABLE :: var_in_output(:)
  
  ! derived variables
  !
  INTEGER, PARAMETER :: read_netcdf_broadcast_method  = 1
  INTEGER, PARAMETER :: read_netcdf_distribute_method = 2
  INTEGER :: default_read_method = 2

  INTEGER :: restart_file_type = FILETYPE_NC2

  LOGICAL :: write_initial_state = .true.
  LOGICAL :: write_last_restart  = .false.
  INTEGER :: timeSteps_per_outputStep    = 0

  INTEGER :: n_chkpt           ! number of timesteps between successive checkpoint events
  INTEGER :: n_diag            ! number of timesteps between successive tot_int diag events


  ! currently used by hydrostatic model only
  LOGICAL :: l_outputtime      ! if .true., output is written at the end of the time step.

  LOGICAL :: lmask_boundary(max_dom)    ! flag: true, if interpolation zone should be masked *in output*

  CHARACTER(LEN = 256) :: restart_write_mode

  ! When using the restart write mode "dedicated proc mode", it is
  ! possible to split the restart output into several files, as if
  ! "nrestart_streams" * "num_io_procs" restart processes were
  ! involved. This speeds up the read-in process, since all the files
  ! may then be read in parallel.
  INTEGER :: nrestart_streams

  ! Allows checkpointing (followed by stopping) during runtime triggered by a file named 'stop_icon' in the workdir
  LOGICAL :: checkpoint_on_demand

  ! constants to communicate which restart writing MODULE to USE
  ENUM, BIND(C)
    ENUMERATOR :: kSyncRestartModule = 1, kAsyncRestartModule, kMultifileRestartModule
  END ENUM

  CHARACTER(*), PARAMETER :: modname = "mo_io_config"

  INTEGER, PARAMETER :: ALL_WORKERS_INVOLVED = -1


CONTAINS

  !! Precomputation of derived type collecting logical variables indicating whether
  !! optional diagnostics are requested in the output namelists
  !!
  !! Replaces repeated calculations of the same that used to be scattered around various places in the model code
  !!
  SUBROUTINE init_var_in_output(n_dom, lnwp)

    INTEGER, INTENT(in)  :: n_dom  ! number of model domains
    LOGICAL, INTENT(in)  :: lnwp   ! true if ICON runs in NWP mode, implying that the full set of variables 
                                   ! needs to be computed

    INTEGER :: jg, jgr

    ALLOCATE(var_in_output(n_dom))
    !$ACC ENTER DATA CREATE(var_in_output)

    DO jg=1,n_dom
      var_in_output(jg)%pres_msl = is_variable_in_output(var_name="pres_msl") .OR. &
        &                          is_variable_in_output(var_name="psl_m")
      var_in_output(jg)%omega    = is_variable_in_output(var_name="omega")    .OR. &
        &                          is_variable_in_output(var_name="wap_m")
      var_in_output(jg)%res_soilwatb = is_variable_in_output_dom(var_name="resid_wso", jg=jg)
      var_in_output(jg)%ddt_temp_drag = is_variable_in_output_dom(var_name="ddt_temp_drag", jg=jg)
    END DO


    IF (lnwp) THEN
      DO jg=1,n_dom
        var_in_output(jg)%rh          = is_variable_in_output_dom(var_name="rh", jg=jg)
        var_in_output(jg)%pv          = is_variable_in_output_dom(var_name="pv", jg=jg)
        var_in_output(jg)%sdi2        = is_variable_in_output_dom(var_name="sdi2", jg=jg)
        var_in_output(jg)%lpi         = is_variable_in_output_dom(var_name="lpi", jg=jg)
        var_in_output(jg)%lpi_max     = is_variable_in_output_dom(var_name="lpi_max", jg=jg)
        var_in_output(jg)%lpi_con     = is_variable_in_output_dom(var_name="lpi_con", jg=jg)
        var_in_output(jg)%mlpi_con    = is_variable_in_output_dom(var_name="mlpi_con", jg=jg)
        var_in_output(jg)%lpi_con_max = is_variable_in_output_dom(var_name="lpi_con_max", jg=jg)
        var_in_output(jg)%mlpi_con_max= is_variable_in_output_dom(var_name="mlpi_con_max", jg=jg)
        var_in_output(jg)%lfd_con     = is_variable_in_output_dom(var_name="lfd_con", jg=jg)
        var_in_output(jg)%lfd_con_max = is_variable_in_output_dom(var_name="lfd_con_max", jg=jg)
        var_in_output(jg)%koi         = is_variable_in_output_dom(var_name="koi", jg=jg)
        var_in_output(jg)%ceiling     = is_variable_in_output_dom(var_name="ceiling", jg=jg)
        var_in_output(jg)%vis         = is_variable_in_output_dom(var_name="vis", jg=jg)
        var_in_output(jg)%hbas_sc     = is_variable_in_output_dom(var_name="hbas_sc", jg=jg)
        var_in_output(jg)%htop_sc     = is_variable_in_output_dom(var_name="htop_sc", jg=jg)
        var_in_output(jg)%inversion_height = is_variable_in_output_dom(var_name="inversion_height", jg=jg)
        var_in_output(jg)%twater      = is_variable_in_output_dom(var_name="twater", jg=jg)
        var_in_output(jg)%q_sedim     = is_variable_in_output_dom(var_name="q_sedim", jg=jg)
        var_in_output(jg)%tcond_max   = is_variable_in_output_dom(var_name="tcond_max", jg=jg)
        var_in_output(jg)%tcond10_max = is_variable_in_output_dom(var_name="tcond10_max", jg=jg)
        var_in_output(jg)%uh_max_low  = is_variable_in_output_dom(var_name="uh_max_low", jg=jg)
        var_in_output(jg)%uh_max_med  = is_variable_in_output_dom(var_name="uh_max_med", jg=jg)
        var_in_output(jg)%uh_max      = is_variable_in_output_dom(var_name="uh_max", jg=jg)
        var_in_output(jg)%vorw_ctmax  = is_variable_in_output_dom(var_name="vorw_ctmax", jg=jg)
        var_in_output(jg)%w_ctmax     = is_variable_in_output_dom(var_name="w_ctmax", jg=jg)
        var_in_output(jg)%dbz         = is_variable_in_output_dom(var_name="dbz", jg=jg)
        var_in_output(jg)%dbz850      = is_variable_in_output_dom(var_name="dbz_850", jg=jg)
        var_in_output(jg)%dbzcmax     = is_variable_in_output_dom(var_name="dbz_cmax", jg=jg)
        var_in_output(jg)%dbzctmax    = is_variable_in_output_dom(var_name="dbz_ctmax", jg=jg)
        var_in_output(jg)%dbzlmx_low  = is_variable_in_output_dom(var_name="dbzlmx_low", jg=jg)
        var_in_output(jg)%echotop     = is_variable_in_output_dom(var_name="echotop", jg=jg)
        var_in_output(jg)%echotopinm  = is_variable_in_output_dom(var_name="echotopinm", jg=jg)
        var_in_output(jg)%smi         = is_variable_in_output_dom(var_name="smi", jg=jg)
        var_in_output(jg)%dursun      = is_variable_in_output_dom(var_name="dursun", jg=jg)
        var_in_output(jg)%dursun_m    = is_variable_in_output_dom(var_name="dursun_m", jg=jg)
        var_in_output(jg)%dursun_r    = is_variable_in_output_dom(var_name="dursun_r", jg=jg)
        var_in_output(jg)%snow_melt   = is_variable_in_output_dom(var_name="snow_melt", jg=jg)
        var_in_output(jg)%dhail_mx    = is_variable_in_output_dom(var_name="dhail_mx", jg=jg)
        var_in_output(jg)%dhail_av    = is_variable_in_output_dom(var_name="dhail_av", jg=jg)
        var_in_output(jg)%dhail_sd    = is_variable_in_output_dom(var_name="dhail_sd", jg=jg)
        var_in_output(jg)%wshear_u    = is_variable_in_output_dom(var_name="wshear_u", jg=jg)
        var_in_output(jg)%wshear_v    = is_variable_in_output_dom(var_name="wshear_v", jg=jg)
        var_in_output(jg)%lapserate   = is_variable_in_output_dom(var_name="lapse_rate", jg=jg)
        var_in_output(jg)%mconv       = is_variable_in_output_dom(var_name="mconv", jg=jg)
        var_in_output(jg)%srh         = is_variable_in_output_dom(var_name="srh", jg=jg)
        var_in_output(jg)%tot_pr_max  = is_variable_in_output_dom(var_name="tot_pr_max", jg=jg)
        var_in_output(jg)%cape_mu     = is_variable_in_output_dom(var_name="cape_mu", jg=jg)
        var_in_output(jg)%cin_mu      = is_variable_in_output_dom(var_name="cin_mu", jg=jg)
        var_in_output(jg)%cape_3km    = is_variable_in_output_dom(var_name="cape_3km", jg=jg)
        var_in_output(jg)%lfc_ml      = is_variable_in_output_dom(var_name="lfc_ml", jg=jg)
        var_in_output(jg)%lcl_ml      = is_variable_in_output_dom(var_name="lcl_ml", jg=jg)
        var_in_output(jg)%si          = is_variable_in_output_dom(var_name="si", jg=jg)
        var_in_output(jg)%sli         = is_variable_in_output_dom(var_name="sli", jg=jg)
        var_in_output(jg)%swiss12     = is_variable_in_output_dom(var_name="swiss12", jg=jg)
        var_in_output(jg)%swiss00     = is_variable_in_output_dom(var_name="swiss00", jg=jg)
        var_in_output(jg)%cloudtop    = is_variable_in_output_dom(var_name="cloudtop", jg=jg)
        var_in_output(jg)%hpbl        = is_variable_in_output_dom(var_name="hpbl", jg=jg)
        var_in_output(jg)%aod_550nm   = is_variable_in_output_dom(var_name="aod_550nm", jg=jg)

        ! add vars for global mean calculations
        var_in_output(jg)%tas_gmean   = is_variable_in_output_dom(var_name="tas_gmean", jg=jg)
        var_in_output(jg)%rsdt_gmean  = is_variable_in_output_dom(var_name="rsdt_gmean", jg=jg)
        var_in_output(jg)%rsut_gmean  = is_variable_in_output_dom(var_name="rsut_gmean", jg=jg)
        var_in_output(jg)%rlut_gmean  = is_variable_in_output_dom(var_name="rlut_gmean", jg=jg)
        var_in_output(jg)%prec_gmean  = is_variable_in_output_dom(var_name="prec_gmean", jg=jg)
        var_in_output(jg)%evap_gmean  = is_variable_in_output_dom(var_name="evap_gmean", jg=jg)
        var_in_output(jg)%pme_gmean   = is_variable_in_output_dom(var_name="pme_gmean", jg=jg)
        var_in_output(jg)%radtop_gmean= is_variable_in_output_dom(var_name="radtop_gmean", jg=jg)

        ! Check for special case: SMI is not in one of the output lists but it is part of a output group.
        ! In this case, the group can not be checked, as the connection between SMI and the group will be
        ! established during the add_var call. However, add_var for SMI will only be called if l_smi =.true.
        ! As a crutch, a character array containing the output groups of SMI from mo_lnd_nwp_config is used
        ! here and also at the add_var call.
        ! The check loops through the output groups. It has to be checked if l_smi is already .true., to not
        ! overwrite an existing .true. with a false. 

        IF (.NOT. var_in_output(jg)%smi) THEN 
          ! Check for output groups containing SMI
          DO jgr = 1,SIZE(groups_smi)
            IF (.NOT. var_in_output(jg)%smi) THEN
              var_in_output(jg)%smi = is_variable_in_output_dom(&
                                      var_name='group:'//TRIM(groups_smi(jgr)), jg=jg)
            END IF
          END DO
        END IF
      END DO
    END IF


    ! diagnostics for the horizontal wind tendencies in the dynamical core
    DO jg=1,n_dom
      var_in_output(jg)%ddt_vn_dyn    = is_variable_in_output_dom(var_name="ddt_vn_dyn", jg=jg)
      var_in_output(jg)%ddt_ua_dyn    = is_variable_in_output_dom(var_name="ddt_ua_dyn", jg=jg)
      var_in_output(jg)%ddt_va_dyn    = is_variable_in_output_dom(var_name="ddt_va_dyn", jg=jg)
      var_in_output(jg)%ddt_vn_dmp    = is_variable_in_output_dom(var_name="ddt_vn_dmp", jg=jg)
      var_in_output(jg)%ddt_ua_dmp    = is_variable_in_output_dom(var_name="ddt_ua_dmp", jg=jg)
      var_in_output(jg)%ddt_va_dmp    = is_variable_in_output_dom(var_name="ddt_va_dmp", jg=jg)
      var_in_output(jg)%ddt_vn_hdf    = is_variable_in_output_dom(var_name="ddt_vn_hdf", jg=jg)
      var_in_output(jg)%ddt_ua_hdf    = is_variable_in_output_dom(var_name="ddt_ua_hdf", jg=jg)
      var_in_output(jg)%ddt_va_hdf    = is_variable_in_output_dom(var_name="ddt_va_hdf", jg=jg)
      var_in_output(jg)%ddt_vn_adv    = is_variable_in_output_dom(var_name="ddt_vn_adv", jg=jg)
      var_in_output(jg)%ddt_ua_adv    = is_variable_in_output_dom(var_name="ddt_ua_adv", jg=jg)
      var_in_output(jg)%ddt_va_adv    = is_variable_in_output_dom(var_name="ddt_va_adv", jg=jg)
      var_in_output(jg)%ddt_vn_cor    = is_variable_in_output_dom(var_name="ddt_vn_cor", jg=jg)
      var_in_output(jg)%ddt_ua_cor    = is_variable_in_output_dom(var_name="ddt_ua_cor", jg=jg)
      var_in_output(jg)%ddt_va_cor    = is_variable_in_output_dom(var_name="ddt_va_cor", jg=jg)
      var_in_output(jg)%ddt_vn_pgr    = is_variable_in_output_dom(var_name="ddt_vn_pgr", jg=jg)
      var_in_output(jg)%ddt_ua_pgr    = is_variable_in_output_dom(var_name="ddt_ua_pgr", jg=jg)
      var_in_output(jg)%ddt_va_pgr    = is_variable_in_output_dom(var_name="ddt_va_pgr", jg=jg)
      var_in_output(jg)%ddt_vn_phd    = is_variable_in_output_dom(var_name="ddt_vn_phd", jg=jg)
      var_in_output(jg)%ddt_ua_phd    = is_variable_in_output_dom(var_name="ddt_ua_phd", jg=jg)
      var_in_output(jg)%ddt_va_phd    = is_variable_in_output_dom(var_name="ddt_va_phd", jg=jg)
      var_in_output(jg)%ddt_vn_iau    = is_variable_in_output_dom(var_name="ddt_vn_iau", jg=jg)
      var_in_output(jg)%ddt_ua_iau    = is_variable_in_output_dom(var_name="ddt_ua_iau", jg=jg)
      var_in_output(jg)%ddt_va_iau    = is_variable_in_output_dom(var_name="ddt_va_iau", jg=jg)
      var_in_output(jg)%ddt_vn_ray    = is_variable_in_output_dom(var_name="ddt_vn_ray", jg=jg)
      var_in_output(jg)%ddt_ua_ray    = is_variable_in_output_dom(var_name="ddt_ua_ray", jg=jg)
      var_in_output(jg)%ddt_va_ray    = is_variable_in_output_dom(var_name="ddt_va_ray", jg=jg)
      var_in_output(jg)%ddt_vn_grf    = is_variable_in_output_dom(var_name="ddt_vn_grf", jg=jg)
      var_in_output(jg)%ddt_ua_grf    = is_variable_in_output_dom(var_name="ddt_ua_grf", jg=jg)
      var_in_output(jg)%ddt_va_grf    = is_variable_in_output_dom(var_name="ddt_va_grf", jg=jg)
    END DO

    !$ACC UPDATE DEVICE(var_in_output) ASYNC(1)

  END SUBROUTINE init_var_in_output


  !! Set up derived components of the I/O config state
  !!
  !! Set up derived components of the I/O config state. This routine is
  !! called, after all namelists have been read and a synoptic consistency
  !! check has been done.
  !!
  SUBROUTINE configure_io()

    !-----------------------------------------------------------------------

    ! number of timesteps between successive checkpoint events
    n_chkpt = n_checkpoints()

    ! number of timesteps between successive tot_int diag events
    n_diag  = n_diags()

  END SUBROUTINE configure_io



  !----------------------------------------------------------------------------------
   FUNCTION n_checkpoints()

     INTEGER :: n_checkpoints

     n_checkpoints = NINT(dt_checkpoint/dtime)  ! write restart files
     IF (n_checkpoints == 0) n_checkpoints = HUGE(1)
   END FUNCTION n_checkpoints
  !----------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------
   FUNCTION n_diags()

     INTEGER :: n_diags

     n_diags  = MAX(1,NINT(dt_diag/dtime)) ! number of: diagnose of total integrals
   END FUNCTION n_diags
  !----------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------

   FUNCTION is_checkpoint_time(current_step, n_checkpoints, n_steps) RESULT(l_checkpoint)
     INTEGER, INTENT(IN)            :: current_step, n_checkpoints
     INTEGER, INTENT(IN), OPTIONAL  :: n_steps

     LOGICAL              :: l_checkpoint

     IF (n_checkpoints == 0) THEN
       l_checkpoint = .FALSE.
     ELSE
       IF (PRESENT(n_steps)) THEN
         IF ( MOD(current_step,n_checkpoints)==0 .AND. current_step/=n_steps ) THEN
           l_checkpoint = .TRUE.
         ELSE
           l_checkpoint = .FALSE.
         END IF
       ELSE
         IF ( MOD(current_step,n_checkpoints)==0 ) THEN
           l_checkpoint = .TRUE.
         ELSE
           l_checkpoint = .FALSE.
         END IF
       END IF
     END IF
   END FUNCTION is_checkpoint_time
  !----------------------------------------------------------------------------------

  !! Decides about diagnostic computation of total integrals, which
  !! is performed in "supervise_total_integrals_nh"
  !! Total integrals are computed
  !! - at the first time step (or the first time step after restart)
  !! - if (MOD(current_step,n_diag) == 0)
  !! - at the very last time step
  !!
  FUNCTION is_totint_time(current_step, restart_step, n_diag, n_steps)

    INTEGER, INTENT(IN) :: current_step !< current time step number
    INTEGER, INTENT(IN) :: restart_step !< time step for which the restart file was produced
                                        !< rfile_step+1: first step after restart
    INTEGER, INTENT(IN) :: n_steps      !< total number of time steps
    INTEGER, INTENT(IN) :: n_diag       !< number of timesteps between successive calls

    LOGICAL :: is_totint_time           ! Result

    ! local variables
    INTEGER :: kstep                    ! time step number relative to restart step

    kstep = current_step - restart_step

    is_totint_time = ((kstep == 1)                        .OR. &
      &              (MOD(current_step,n_diag) == 0)      .OR. &
      &              (kstep==n_steps))                    .AND.&
      &              (current_step > 0 )
  END FUNCTION is_totint_time

  ! Inquire the parameters for restart writing.
  !
  ! opt_dedicatedProcCount: The number of processes that are split off
  !                         the MPI_Communicator to serve as dedicated
  !                         restart processes.
  !
  ! opt_restartProcCount:   The number of processes that actually
  !                         perform restart writing. This IS never zero,
  !                         especially NOT when opt_dedicatedProcCount
  !                         IS zero.
  !
  ! opt_restartModule:      The id of the MODULE that handles the restart
  !                         writing. Possible values are
  !                         kSyncRestartModule, kAsyncRestartModule, AND
  !                         kMultifileRestartModule.
  !
  ! opt_lDedicatedProcMode: Whether the restart processes are split
  !                         off the MPI_Communicator OR are a subset
  !                         of mpi_comm_work.
  !
  SUBROUTINE restartWritingParameters(opt_dedicatedProcCount, &
    &                                 opt_restartProcCount,   &
    &                                 opt_restartModule,      &
    &                                 opt_lDedicatedProcMode, &
    &                                 opt_nrestart_streams)

    INTEGER, INTENT(OUT), OPTIONAL :: opt_dedicatedProcCount, &
      &                               opt_restartProcCount,   &
      &                               opt_restartModule,      &
      &                               opt_nrestart_streams
    LOGICAL, INTENT(OUT), OPTIONAL :: opt_lDedicatedProcMode

    LOGICAL, SAVE :: cacheValid = .FALSE., lDedicatedProcMode = .FALSE.
    INTEGER, SAVE :: dedicatedProcCount = -1, &
      &              restartProcCount   = -1, &
      &              restartModule      = -1, &
      &              nrestartStreams    = -1
    CHARACTER(:), ALLOCATABLE :: errorMessage
    CHARACTER(*), PARAMETER :: routine = modname//":restartWritingParameters"

    ! If this IS the first CALL of this FUNCTION, analyze the namelist
    ! parameters AND cache the RESULT, sanity checking the settings.
    IF(.NOT.cacheValid) THEN
      IF(num_restart_procs < 0) THEN
        ! No matter what, negative process counts are illegal.
        errorMessage = "illegal value of namelist parameter num_restart_procs: value must not be negative"
        
      ELSE IF(restart_write_mode == "") THEN
        ! No restart_write_mode given, so we fall back to the old
        ! behavior of switching between sync/async restart mode based
        ! on num_restart_procs for backwards compatibility.
        dedicatedProcCount = num_restart_procs
        restartProcCount   = MAX(1, num_restart_procs)
        lDedicatedProcMode = num_restart_procs > 0
        nrestartStreams    = 1
        IF(lDedicatedProcMode) THEN
          restartModule = kAsyncRestartModule
        ELSE
          restartModule = kSyncRestartModule
        END IF
        
      ELSE IF(restart_write_mode == "sync") THEN
        IF(num_restart_procs /= 0) THEN
          errorMessage = "inconsistent namelist parameters: num_restart_procs must be zero OR unset &
            &for the given restart_write_mode"
        END IF
        dedicatedProcCount = 0
        restartProcCount   = 1
        restartModule      = kSyncRestartModule
        lDedicatedProcMode = .FALSE.
        nrestartStreams    = 1
        
      ELSE IF(restart_write_mode == "async") THEN
        IF(num_restart_procs == 0) THEN
          errorMessage = "inconsistent namelist parameters: num_restart_procs must be non-zero for &
            &the given restart_write_mode"
        END IF
        dedicatedProcCount = num_restart_procs
        restartProcCount   = num_restart_procs
        restartModule      = kAsyncRestartModule
        lDedicatedProcMode = .TRUE.
        nrestartStreams    = 1

      ELSE IF(restart_write_mode == "joint procs multifile") THEN
        dedicatedProcCount = 0
        
        ! if not set otherwise (num_restart_procs=0), set the
        ! number of restart PEs to the number of worker
        ! PEs... this is done later, since the number of workers
        ! is not yet available.
        IF (num_restart_procs == 0) THEN
          restartProcCount = ALL_WORKERS_INVOLVED
        ELSE
          restartProcCount = num_restart_procs
        END IF
        restartModule      = kMultifileRestartModule
        lDedicatedProcMode = .FALSE.
        nrestartStreams    = 1
        
      ELSE IF(restart_write_mode == "dedicated procs multifile") THEN
        IF(num_restart_procs == 0) THEN
          errorMessage = "inconsistent namelist parameters: num_restart_procs must be non-zero for &
            &the given restart_write_mode"
        END IF
        dedicatedProcCount = num_restart_procs
        restartProcCount   = num_restart_procs
        restartModule      = kMultifileRestartModule
        lDedicatedProcMode = .TRUE.
        nrestartStreams    = nrestart_streams
        
      ELSE
        errorMessage = "illegal value of namelist parameter restart_write_mode: expected one of 'sync', 'async', &
          &'joint procs multifile', or 'dedicated procs multifile'"
        
      END IF
      IF(ALLOCATED(errorMessage)) THEN
        CALL finish(routine, errorMessage//" (got restart_write_mode = '"//TRIM(restart_write_mode)// &
          &"', num_restart_procs = "//TRIM(int2string(num_restart_procs))//")")
      END IF
      cacheValid = .TRUE.
    END IF
    
    ! Ok, the cache IS up to date. What info did the caller want again?
    IF(PRESENT(opt_dedicatedProcCount))  opt_dedicatedProcCount = dedicatedProcCount;
    IF(PRESENT(opt_restartProcCount))    opt_restartProcCount   = restartProcCount;
    IF(PRESENT(opt_restartModule))       opt_restartModule      = restartModule;
    IF(PRESENT(opt_lDedicatedProcMode))  opt_lDedicatedProcMode = lDedicatedProcMode;
    IF(PRESENT(opt_nrestart_streams))    opt_nrestart_streams   = nrestartStreams

  END SUBROUTINE restartWritingParameters

END MODULE mo_io_config
