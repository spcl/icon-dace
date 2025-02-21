!
! Contains subroutines for initializing the AES physics package.
!
!
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

!----------------------------
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

!----------------------------

MODULE mo_aes_phy_init

  ! infrastructure
  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish, message, message_text, print_value
  USE mtime,                   ONLY: datetime, OPERATOR(>), OPERATOR(==)
  USE mo_read_interface,       ONLY: openInputFile, closeFile, read_2D, &
    &                                t_stream_id, on_cells
  USE mo_timer,                ONLY: timers_level, timer_start, timer_stop, &
    &                                timer_prep_aes_phy

  ! model configuration
  USE mo_impl_constants,       ONLY: min_rlcell_int, grf_bdywidth_c, vname_len
  USE mo_parallel_config,      ONLY: nproma
  USE mo_master_config,        ONLY: isrestart
  USE mo_run_config,           ONLY: ltestcase, msg_level,                        &
    &                                iqv, iqc, iqi, iqs, iqr, iqg, iqm_max,       &
    &                                iqh, iqni,iqnr,iqns,iqng,iqnh, iqnc,ininact, &
    &                                iqt, io3, ico2, ich4, in2o, ntracer
  USE mo_advection_config,     ONLY: advection_config

  ! horizontal grid and indices
  USE mo_model_domain,         ONLY: t_patch
  USE mo_loopindices,          ONLY: get_indices_c

  ! test cases
  USE mo_nh_testcases_nml,     ONLY: nh_test_name, ape_sst_case, th_cbl, tpe_temp
  USE mo_ape_params,           ONLY: ape_sst
  USE mo_physical_constants,   ONLY: tmelt, Tf, albedoW, amd, amo3, zemiss_def

  USE mo_sea_ice_nml,          ONLY: albi

  ! aes physics
  USE mo_aes_phy_config,       ONLY: eval_aes_phy_config, eval_aes_phy_tc, print_aes_phy_config, &
    &                                aes_phy_config, aes_phy_tc, dt_zero
  USE mo_aes_phy_memory,       ONLY: prm_field, t_aes_phy_field, &
    &                                prm_tend,  t_aes_phy_tend

  ! radiation
  USE mo_aes_rad_config,       ONLY: eval_aes_rad_config, print_aes_rad_config, aes_rad_config

  ! vertical diffusion
  USE mo_aes_vdf_config,       ONLY: eval_aes_vdf_config, print_aes_vdf_config, aes_vdf_config
  USE mo_turb_vdiff,           ONLY: vdiff_init
  ! USE mo_vdf_diag_smag,        ONLY: sfc_exchange_coefficients

  ! land surface
  USE mo_jsb_model_init,       ONLY: jsbach_init
  USE mo_jsb_interface,        ONLY: jsbach_get_var
  ! USE mo_phy_schemes,          ONLY: register_exchange_coefficients_procedure
  USE mo_ext_data_state,       ONLY: ext_data

  ! carbon cycle
  USE mo_ccycle_config,        ONLY: print_ccycle_config, ccycle_config

  ! cumulus convection
  USE mo_convect_tables,       ONLY: init_convect_tables
  USE mo_aes_convect_tables,   ONLY: init_aes_convect_tables => init_convect_tables

  ! cloud optical properties
  USE mo_aes_cop_config,       ONLY: print_aes_cop_config, aes_cop_config

  ! two-moment bulk microphysics
  USE mo_2mom_mcrph_driver,    ONLY: two_moment_mcrph_init

  ! cloud cover
  USE mo_aes_cov_config,       ONLY: eval_aes_cov_config, print_aes_cov_config

  ! WMO tropopause
  USE mo_aes_wmo_config,       ONLY: eval_aes_wmo_config, print_aes_wmo_config

  ! air-sea-land interface
  USE mo_aes_sfc_indices,      ONLY: nsfc_type, iwtr, iice, ilnd, init_sfc_indices

  ! for AMIP boundary conditions
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, calculate_time_interpolation_weights
  USE mo_bc_sst_sic,           ONLY: read_bc_sst_sic, bc_sst_sic_time_interpolation
  USE mo_bc_greenhouse_gases,  ONLY: read_bc_greenhouse_gases, bc_greenhouse_gases_time_interpolation, &
    &                                bc_greenhouse_gases_file_read
  USE mo_bc_aeropt_splumes,    ONLY: setup_bc_aeropt_splumes

  ! for 6hourly sst and ice data
  USE mo_reader_sst_sic,       ONLY: t_sst_sic_reader
  USE mo_interpolate_time,     ONLY: t_time_intp

  USE mo_rte_rrtmgp_setup,     ONLY: rte_rrtmgp_basic_setup
  USE mo_rte_rrtmgp_interface, ONLY: pressure_scale, droplet_scale
  USE mo_lcariolle,            ONLY: lcariolle_init_o3, lcariolle_init, &
    &l_cariolle_initialized_o3, t_avi, t_time_interpolation

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: init_aes_phy_params, init_aes_phy_external, init_aes_phy_field
  PUBLIC  :: init_o3_lcariolle
  PUBLIC  :: sst_intp, sic_intp, sst_sic_reader

  CHARACTER(len=*), PARAMETER :: modname = 'mo_aes_phy_init'
  TYPE(t_sst_sic_reader), TARGET :: sst_sic_reader
  TYPE(t_time_intp)      :: sst_intp
  TYPE(t_time_intp)      :: sic_intp
  REAL(wp), ALLOCATABLE  :: sst_dat(:,:,:,:)
  REAL(wp), ALLOCATABLE  :: sic_dat(:,:,:,:)

CONTAINS
  !>
  !! Top-level routine for initialization of AES physics.
  !! It calls a series of subroutines to initialize tunable parameters,
  !! lookup tables, and the physics state vectors "prm_field" and "prm_tend".
  !!
  SUBROUTINE init_aes_phy_params( p_patch )

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch(:)

    INTEGER :: nhydromet, ntrac
    INTEGER :: ng, jg
    INTEGER :: nlev

    LOGICAL :: lany

    ! Shortcuts to components of aes_rad_config
    !
    IF (timers_level > 1) CALL timer_start(timer_prep_aes_phy)

    !-------------------------------------------------------------------
    ! Initialize parameters and lookup tables
    !-------------------------------------------------------------------

    ng   = SIZE(p_patch)
    nlev = p_patch(1)%nlev

    ! Diagnostics (all time steps)
    ! ----------------------------
    
    ! For surface processes
    ! nsfc_type, iwtr, etc. are set in this subroutine.
    ! See mo_sfc_indices.f90 for further details.
    !
    CALL init_sfc_indices( nh_test_name )
    
    ! Lookup tables for saturation vapour pressure
    !
    CALL init_convect_tables
    CALL init_aes_convect_tables 


    ! AES physics time control
    ! --------------------------

    ! Evaluate the AES physics configuration variables aes_phy_config(:)
    ! and the derived time control variables aes_phy_tc(:) on all grids
    ! and for all controled processes.
    !
    CALL  eval_aes_phy_config(ng)
    CALL  eval_aes_phy_tc    (ng)
    CALL print_aes_phy_config(ng)


    ! Set tracer indices for physics
    ! ------------------------------

    CALL init_aes_phy_tracer(ng)

    ! Parameterizations (with time control)
    ! -------------------------------------

    ! radiation
    !
    lany=.FALSE.
    DO jg = 1,ng
       lany = lany .OR. (aes_phy_tc(jg)%dt_rad > dt_zero)
    END DO
    IF (lany) THEN
      !
      ! Radiation configuration
      CALL  eval_aes_rad_config(ng)
      CALL print_aes_rad_config(ng)
      !
      ! Cloud optical properties
      CALL print_aes_cop_config(ng)
      !
      ! Radiation constants for gas and cloud optics
      CALL rte_rrtmgp_basic_setup(nproma, nlev, pressure_scale, droplet_scale,       &
        &                    aes_cop_config(1)%cinhoml1 ,aes_cop_config(1)%cinhoml2, &
        &                    aes_cop_config(1)%cinhoml3 ,aes_cop_config(1)%cinhomi)
    END IF

    ! vertical turbulent mixing and surface
    !
    lany=.FALSE.
    DO jg = 1,ng
       lany = lany .OR. (aes_phy_tc(jg)%dt_vdf > dt_zero)
    END DO
    IF (lany) THEN
      !
      CALL  eval_aes_vdf_config
      CALL print_aes_vdf_config(ng)
      !
      ! Allocate memory for the tri-diagonal solver needed by the implicit
      ! time stepping scheme; Compute time-independent parameters.
      !
      ! vdiff diffuses only water vapor (index iqv), two hydro meteors
      ! cloud water (index iqc) and cloud ice (index iqi), and further
      ! tracers, which are supposed to be gases or suspended particles.
      ! These additional ntrac tracers are supposed to be stored with
      ! indices in the range [iqt,ntracer].
      !
      ! Precipitating hydrometeors (rain, snow, graupel) are not diffused
      ! 
      nhydromet = 2              ! diffuse two hydro meteor specied: cloud water and ice
      ntrac = ntracer - iqt + 1  ! and ntrac further species
      !
      CALL vdiff_init( nhydromet, ntrac )
      !
      ! JSBACH land processes
      !

      IF (ilnd <= nsfc_type .AND. ANY(aes_phy_config(:)%ljsb)) THEN
        DO jg=1,ng
          IF (aes_phy_config(jg)%ljsb) THEN 
            CALL jsbach_init(jg)
          END IF
        END DO ! jg
        ! IF (aes_vdf_config(1)%use_tmx) THEN
        !   CALL register_exchange_coefficients_procedure(sfc_exchange_coefficients)
        ! END IF
      END IF ! 

    ENDIF

    ! carbon cycle
    !
    CALL print_ccycle_config

    ! cloud microphysics (mig - graupel)
    !
    ! No setting of parameters necessary - nothing to do
    !
    ! cloud microphysics (two moment bulk microphysics)
    !
    lany=.FALSE.
    DO jg = 1,ng
       lany = lany .OR. (aes_phy_tc(jg)%dt_two > dt_zero)
    END DO
    IF (lany) THEN
      jg=1
      ! original atm_phy_nwp_config(jg)%inwp_gscp equals 4
      CALL two_moment_mcrph_init(igscp=4, msg_level=msg_level )
    END IF

    ! cloud cover diagnostics
    !
    CALL  eval_aes_cov_config(ng)
    CALL print_aes_cov_config(ng)

    ! WMO tropopause diagnostics
    !
    CALL  eval_aes_wmo_config(ng)
    CALL print_aes_wmo_config(ng)

    ! Cariolle linearized o3 chemistry
    !
    lany=.FALSE.
    DO jg = 1,ng
       lany = lany .OR. (aes_phy_tc(jg)%dt_car > dt_zero)
    END DO
    IF (lany) THEN
      IF (io3 < iqt .OR. ntracer < io3) THEN
        CALL finish('init_aes_phy: mo_aes_phy_init.f90', &
                   &'cannot find an ozone tracer - abort')
      END IF
      IF (ng > 1) THEN
        CALL finish('init_aes_phy: mo_aes_phy_init.f90', &
                   &'Cariolle initialization not ready for ng>1')
      END IF
      CALL lcariolle_init()
    END IF

    IF (timers_level > 1) CALL timer_stop(timer_prep_aes_phy)

  END SUBROUTINE init_aes_phy_params


  SUBROUTINE init_aes_phy_tracer(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg, jt
    LOGICAL             :: lany
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_aes_phy_tracer'

    ! Set the indices for specific tracers, if they occur among the named tracers.
    !
    iqm_max = 0 ! number of water species tracers
    !
    DO jt=1,advection_config(1)%nname
       !
       SELECT CASE (TRIM(advection_config(1)%tracer_names(jt)))
       CASE('qv','hus')
          iqv=jt
          iqm_max=iqm_max+1
          advection_config(:)%cfstd_names(jt)='specific_humidity'
          advection_config(:)% long_names(jt)='specific humidity'
       CASE('qc','clw')
          iqc=jt
          iqm_max=iqm_max+1
          advection_config(:)%cfstd_names(jt)='cloud_liquid_water'
          advection_config(:)% long_names(jt)='cloud liquid water'
       CASE('qi','cli')
          iqi=jt
          iqm_max=iqm_max+1
          advection_config(:)%cfstd_names(jt)='cloud_ice'
          advection_config(:)% long_names(jt)='cloud ice'
       CASE('qr')
          iqr=jt
          iqm_max=iqm_max+1
          advection_config(:)%cfstd_names(jt)='rain'
          advection_config(:)% long_names(jt)='rain'
       CASE('qs')
          iqs=jt
          iqm_max=iqm_max+1
          advection_config(:)%cfstd_names(jt)='snow'
          advection_config(:)% long_names(jt)='snow'
       CASE('qg')
          iqg=jt
          iqm_max=iqm_max+1
          advection_config(:)%cfstd_names(jt)='graupel'
          advection_config(:)% long_names(jt)='graupel'
       CASE('qh')
          iqh=jt
          iqm_max=iqm_max+1
          advection_config(:)%cfstd_names(jt)='hail'
          advection_config(:)% long_names(jt)='hail'
       CASE('qnc')
          iqnc=jt
       CASE('qni')
          iqni=jt
       CASE('qnr')
          iqnr=jt
       CASE('qns')
          iqns=jt
       CASE('qng')
          iqng=jt
       CASE('qnh')
          iqnh=jt
       CASE('ninact')
          ininact=jt
       CASE('o3')
          io3=jt
          advection_config(:)%cfstd_names(jt)='ozone'
          advection_config(:)% long_names(jt)='ozone'
       CASE('co2')
          ico2=jt
          advection_config(:)%cfstd_names(jt)='carbon_dioxide'
          advection_config(:)% long_names(jt)='carbon dioxide'
       CASE('ch4')
          ich4=jt
          advection_config(:)%cfstd_names(jt)='methane'
          advection_config(:)% long_names(jt)='methane'
       CASE('n2o')
          in2o=jt
          advection_config(:)%cfstd_names(jt)='nitrous_oxide'
          advection_config(:)% long_names(jt)='nitrous oxide'
       CASE DEFAULT
          advection_config(:)%cfstd_names(jt)='tracer_'//TRIM(advection_config(1)%tracer_names(jt)(:vname_len-7))
          advection_config(:)% long_names(jt)='tracer '//TRIM(advection_config(1)%tracer_names(jt)(:vname_len-7))
       END SELECT
    END DO

    iqt=iqm_max+1
    ! Is "graupel" cloud microphysics active?
    ! Then iqv, iqc, iqi, iqr, iqs, and iqg must be non-zero and in {1,2,3,4,5,6}
    lany=.FALSE.
    DO jg = 1,ng
       lany = lany .OR. (aes_phy_tc(jg)%dt_mig > dt_zero)
    END DO
    IF (lany) THEN
       IF (MIN(iqv,iqc,iqi,iqr,iqs,iqg) == 0) THEN
          CALL finish(routine,           &
               &      'For "Graupel" cloud microphysics, the 6 tracers '// &
               &      'qv/hus, qc/clw, qi/cli, qr, qs, and qg must be ' // &
               &      'included in transport_nml/tracer_names')
       END IF
       IF (MAX(iqv,iqc,iqi,iqr,iqs,iqg) > 6) THEN
          CALL finish(routine,            &
               &      'For "Graupel" cloud microphysics, the 6 tracers ' // &
               &      'qv/hus, qc/clw, qi/cli, qr, qs, and qg must be '  // &
               &      'among the first 6 included in transport_nml/tracer_names')
       END IF
       IF (iqm_max > 6) THEN
          CALL print_value('ATTENTION! '   // &
               &           '"Graupel" cloud microphysics is used with more than 6 ' //       &
               &           'tracers for mass fractions of water substances in air: iqm_max', &
               &           iqm_max, routine=routine)
       END IF
    END IF

    ! Is Two moment bulk cloud microphysics active?
    ! Then iqv, iqc, iqi, iqr, iqs, iqg, iqh, iqnc, iqni, iqnr, iqns, iqng, iqnh, and ininact
    ! must be non-zero and in {1,...,14}
    lany=.FALSE.
    DO jg = 1,ng
       lany = lany .OR. (aes_phy_tc(jg)%dt_two > dt_zero)
    END DO
    IF (lany) THEN
       IF (MIN(iqv,iqc,iqi,iqr,iqs,iqg,iqh,iqnc,iqni,iqnr,iqns,iqng,iqnh,ininact) == 0) THEN
          CALL finish('mo_aes_phy_init:init_aes_phy_tracer',                    &
               &      'For two moment bulk microphysics, 14 tracers must be '// &
               &      'included in transport_nml/tracer_names')
       END IF
       IF (MAX(iqv,iqc,iqi,iqr,iqs,iqg,iqh,iqnc,iqni,iqnr,iqns,iqng,iqnh,ininact) > 14) THEN
          CALL finish('mo_aes_phy_init:init_aes_phy_tracer',                    &
               &      'For two moment bulk microphysics, 14 tracers must be '// &
               &      'among the first 14 included in transport_nml/tracer_names')
       END IF
       IF (iqm_max > 7) THEN
          CALL print_value('mo_aes_phy_init:init_aes_phy_tracer: ATTENTION! '       // &
               &           'two moment bulk microphysics is used with more than 7 ' // &
               &           'tracers for mass fractions of water substances in air: iqm_max',iqm_max)
       END IF
    END IF

    ! Is Cariolle's linearized ozone chemistry active?
    ! Then io3 must be non-zero.
    lany=.FALSE.
    DO jg = 1,ng
       lany = lany .OR. (aes_phy_tc(jg)%dt_car > dt_zero)
    END DO
    IF (lany) THEN
       IF (io3 == 0) THEN
          CALL finish(routine,           &
               &      'For the linearized ozone chemistry of Cariolle, '// &
               &      'the tracer o3 must be included in transport_nml'// &
               &      '/tracer_names')
       END IF
    END IF

    CALL message('','')
    CALL message('','Tracer configuration')
    CALL message('','====================')
    CALL message('','')
    CALL message('','total number of tracers')
    CALL print_value('ntracer',ntracer)
    CALL message('','index variables defined for active tracers')
    IF (iqv  > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqv))//'"  : iqv    ',iqv )
    IF (iqc  > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqc))//'"  : iqc    ',iqc )
    IF (iqi  > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqi))//'"  : iqi    ',iqi )
    IF (iqr  > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqr))//'"   : iqr    ',iqr )
    IF (iqs  > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqs))//'"   : iqs    ',iqs )
    IF (iqg  > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqg))//'"   : iqg    ',iqg )
    IF (iqh  > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqh))//'"   : iqh    ',iqh )
    IF (iqnc > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqnc))//'"  : iqnc   ',iqc )
    IF (iqni > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqni))//'"  : iqni   ',iqi )
    IF (iqnr > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqnr))//'"  : iqnr   ',iqr )
    IF (iqns > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqns))//'"  : iqns   ',iqs )
    IF (iqng > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqng))//'"  : iqng   ',iqg )
    IF (iqnh > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(iqnh))//'"  : iqnh   ',iqh )
    IF (ininact > 0) CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(ininact))//'"  : ininact',iqh )
    IF (io3  > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(io3))//'"  : io3    ',io3 )
    IF (ico2 > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(ico2))//'" : ico2   ',ico2)
    IF (ich4 > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(ich4))//'" : ich4   ',ich4)
    IF (in2o > 0)    CALL print_value('tracer "'//TRIM(advection_config(1)%tracer_names(in2o))//'" : in2o   ',in2o)
    CALL message('','last  index for tracers of mass fraction of water species in air')
    CALL print_value('iqm_max',iqm_max)
    CALL message('','first index for other tracers')
    CALL print_value('iqt    '    ,iqt    )
    CALL message('','number of other tracers')
    CALL print_value('ntrac  ',ntracer-iqt+1)

  END SUBROUTINE init_aes_phy_tracer


  SUBROUTINE init_aes_phy_external( p_patch, mtime_current)

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch(:)
    TYPE(datetime),  INTENT(in), POINTER    :: mtime_current !< Date and time information

    INTEGER :: ng, jg
    LOGICAL :: lany
    TYPE(t_stream_id) :: stream_id

    CHARACTER(len=26+2+3) :: land_frac_fn
    CHARACTER(len=26+2+3) :: land_phys_fn

    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_aes_phy_external'

    IF (timers_level > 1) CALL timer_start(timer_prep_aes_phy)

    ! external data on ICON grids:

    ng = SIZE(p_patch)

    DO jg= 1,ng

      IF (ng > 1) THEN
        WRITE(land_frac_fn, '(a,i2.2,a)') 'bc_land_frac_DOM', jg, '.nc'
        WRITE(land_phys_fn, '(a,i2.2,a)') 'bc_land_phys_DOM', jg, '.nc'
      ELSE
        land_frac_fn = 'bc_land_frac.nc'
        land_phys_fn = 'bc_land_phys.nc'
      ENDIF

      IF (ilnd <= nsfc_type) THEN

        ! land, glacier and lake masks
        !
        IF (aes_vdf_config(jg)%use_tmx) THEN
          WRITE(message_text,'(2a)') 'Read notsea and glac from file ', TRIM(land_frac_fn)
        ELSE
          WRITE(message_text,'(2a)') 'Read notsea, glac and lake from file ', TRIM(land_frac_fn)
        END IF
        CALL message(routine, message_text)
        !
        ! notsea land-sea mask is already read in mo_ext_data_init:read_ext_data_atm because it is needed
        ! early on for topography smoothing
        prm_field(jg)%lsmask(:,:) = ext_data(jg)%atm%fr_land(:,:)
        !
        CALL openInputFile(stream_id, land_frac_fn, p_patch(jg))
        CALL read_2D(stream_id=stream_id, location=on_cells, &
             &          variable_name='glac',               &
             &          fill_array=prm_field(jg)% glac(:,:))
        IF (aes_phy_config(jg)%llake) THEN
          CALL read_2D(stream_id=stream_id, location=on_cells, &
               &          variable_name='lake',               &
               &          fill_array=prm_field(jg)% alake(:,:))
        ELSE
          ! If running without lakes, set lake fraction to zero.
          prm_field(jg)%alake(:,:) = 0._wp
        END IF
        CALL closeFile(stream_id)

        ! For security:
        !  - this statement was necessary for coupled ruby-0 model using land-data created 2020-06
        !  - with land-date created 2021-07 it is not used anymore:
        !prm_field(jg)%lsmask(:,:) = MERGE(1._wp, prm_field(jg)%lsmask(:,:), &
        !                                 prm_field(jg)%lsmask(:,:) > 1._wp - 10._wp*EPSILON(1._wp))
        !
        ! At this point, %lsmask is the fraction of land (incl. glacier and lakes) in the grid box.
        !
        IF (aes_phy_config(jg)%llake) THEN
          !
          ! Substract lake fraction from %lsmask at inner land points (no ocean fraction, %lsmask==1).
          ! Elsewhere (fractional coastal points and ocean points), set alake to zero.
          WHERE (prm_field(jg)%lsmask(:,:) >= 1._wp) ! Inner land point
            prm_field(jg)%lsmask(:,:) = prm_field(jg)%lsmask(:,:) - prm_field(jg)%alake(:,:)
          ELSE WHERE
            prm_field(jg)%alake(:,:) = 0._wp
          END WHERE
        END IF

        ! roughness length and background albedo
        !
        IF (aes_phy_tc(jg)%dt_vdf > dt_zero .OR. aes_phy_tc(jg)%dt_rad > dt_zero) THEN
          !
          CALL openInputFile(stream_id, land_phys_fn, p_patch(jg))
          !
          IF (aes_phy_tc(jg)%dt_vdf > dt_zero .AND. .NOT. aes_vdf_config(jg)%use_tmx) THEN
            !
            WRITE(message_text,'(2a)') 'Read roughness_length from file: ', TRIM(land_phys_fn)
            CALL message(routine, message_text)
            !
            CALL read_2D(stream_id=stream_id, location=on_cells, &
                  &       variable_name='roughness_length',      &
                  &       fill_array=prm_field(jg)% z0m(:,:))
            !
          END IF
          !
          IF (aes_phy_tc(jg)%dt_rad > dt_zero) THEN
            !
            WRITE(message_text,'(2a)') 'Read albedo           from file: ', TRIM(land_phys_fn)
            CALL message(routine, message_text)
            !
            CALL read_2D(stream_id=stream_id, location=on_cells, &
                 &       variable_name='albedo',                &
                 &       fill_array=prm_field(jg)% alb(:,:))
            !
            ! Here surface emissivity should be read from an external file.
            ! But currently this is not available. Instead a default constant
            ! is used as source.
            WRITE(message_text,'(2a)') 'Use default surface emissivity zemiss_def from mo_physical_constants'
            CALL message(routine, message_text)
            !
            prm_field(jg)% emissivity(:,:) = zemiss_def
            !
          END IF
          !
          CALL closeFile(stream_id)
          !
        END IF

      END IF ! (ilnd <= nsfc_type)

    END DO ! jg

    ! external data:


    ! for radiation
    !
    ! Read file for simple plumes aerosol distributions
    !
    lany=.FALSE.
    DO jg = 1,ng
       lany = lany .OR. (aes_phy_tc(jg)%dt_rad > dt_zero)
    END DO
    IF (lany) THEN
      !
      ! parameterized simple plumes of tropospheric aerosols
      !
      IF (ANY(aes_rad_config(:)%irad_aero == 18) .OR. &
        & ANY(aes_rad_config(:)%irad_aero == 19)) THEN
        CALL setup_bc_aeropt_splumes
      END IF

      !
    END IF


    ! for radiation and carbon cycle
    !
    ! Read scenario file for concentrations of CO2, CH4, N2O, CFC11 and CFC12
    ! if radiation is used with any of the gases from the greenhouse gases file
    ! or if the carbon cycle is used with prescribed co2 from this file.
    lany=.FALSE.
    DO jg = 1,ng
       lany = lany .OR. ( aes_phy_tc(jg)%dt_rad > dt_zero .AND.         &
            &             ( aes_rad_config(jg)%irad_co2   == 3 .OR.     &
            &               aes_rad_config(jg)%irad_ch4   == 3 .OR.     &
            &               aes_rad_config(jg)%irad_ch4   ==13 .OR.     &
            &               aes_rad_config(jg)%irad_n2o   == 3 .OR.     &
            &               aes_rad_config(jg)%irad_n2o   ==13 .OR.     &
            &               aes_rad_config(jg)%irad_cfc11 == 3 .OR.     &
            &               aes_rad_config(jg)%irad_cfc12 == 3      ) ) &
            &      .OR. ( ccycle_config(jg)%iccycle  == 2   .AND.       &
            &             ccycle_config(jg)%ico2conc == 4               )
    END DO
    IF (lany) THEN
      !
      ! scenario of well mixed greenhouse gases, horizontally constant
      !
      ! read annual means
      IF (.NOT. bc_greenhouse_gases_file_read) THEN
        CALL read_bc_greenhouse_gases('bc_greenhouse_gases.nc')
      END IF
      ! interpolate to the current date and time, placing the annual means at
      ! the mid points of the current and preceding or following year, if the
      ! current date is in the 1st or 2nd half of the year, respectively.
      CALL bc_greenhouse_gases_time_interpolation(mtime_current)
      !
    END IF

    
    ! for radiation and vertical diffusion
    !
    ! Read sea surface temperature, sea ice concentration and depth
    ! Note: For coupled runs, this is only used for initialization of surface temperatures
    !
    IF (iwtr <= nsfc_type .AND. iice <= nsfc_type) THEN
      !
      IF (.NOT. isrestart()) THEN
        !
        IF (.NOT.  ANY(aes_phy_config(:)%lsstice)) THEN
          !
          ! interpolation weights for linear interpolation of monthly means to the current time
          current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_current)
          !
          DO jg= 1,ng
            !
            IF (aes_phy_tc(jg)%dt_rad > dt_zero .OR. aes_phy_tc(jg)%dt_vdf > dt_zero) THEN
              !
              CALL read_bc_sst_sic(mtime_current%date%year, p_patch(jg))
              !
              CALL bc_sst_sic_time_interpolation(current_time_interpolation_weights                   ,&
                   &                             prm_field(jg)%ts_tile(:,:,iwtr)                      ,&
                   &                             prm_field(jg)%seaice (:,:)                           ,&
                   &                             prm_field(jg)%siced  (:,:)                           ,&
                   &                             p_patch(jg)                                          ,&
                   &                             prm_field(jg)%lsmask(:,:) + prm_field(jg)%alake(:,:) < 1._wp ,&
                   &                             .TRUE. )
              !
            END IF
            !
          END DO
          !
          !
        ELSE
          !
          ! READ 6-hourly sst values (dyamond+- setup, preliminary)
          CALL sst_sic_reader%init(p_patch(1), 'sst-sic-runmean_G.nc')
          CALL sst_intp%init(sst_sic_reader, mtime_current, "SST")
          CALL sst_intp%intp(mtime_current, sst_dat)
          WHERE (sst_dat(:,1,:,1) > 0.0_wp)
            prm_field(1)%ts_tile(:,:,iwtr) = sst_dat(:,1,:,1)
          END WHERE
          !
          CALL sic_intp%init(sst_sic_reader, mtime_current, "SIC")
          CALL sic_intp%intp(mtime_current, sic_dat)
          prm_field(1)%seaice(:,:) = sic_dat(:,1,:,1)
          prm_field(1)%seaice(:,:) = MERGE(0.99_wp, prm_field(1)%seaice(:,:), prm_field(1)%seaice(:,:) > 0.99_wp)
          prm_field(1)%seaice(:,:) = MERGE(0.0_wp, prm_field(1)%seaice(:,:), prm_field(1)%seaice(:,:) <= 0.01_wp)

          ! set ice thickness
          WHERE (prm_field(1)%seaice(:,:) > 0.0_wp)
            prm_field(1)%siced(:,:) = MERGE(2.0_wp, 1.0_wp, p_patch(1)%cells%center(:,:)%lat > 0.0_wp)
          ELSEWHERE
            prm_field(1)%siced(:,:) = 0.0_wp
          ENDWHERE
        !
        END IF
        !
      END IF
      !
    END IF

    IF (timers_level > 1) CALL timer_stop(timer_prep_aes_phy)

  END SUBROUTINE init_aes_phy_external


  !-------------
  !>
  !! Loop over all grid levels and give proper values to some components
  !! of the state vectors "prm_field" and "prm_tend".
  !!
  SUBROUTINE init_aes_phy_field( p_patch        ,&
    &                              temp           )

!FIXME: PGI + OpenMP produce error in this routine... check correctness of parallel code

    TYPE(t_patch)    ,INTENT(in) :: p_patch
    REAL(wp)         ,INTENT(in) :: temp          (:,:,:)

    ! local variables and pointers

    INTEGER  :: jg, ncd, rls, rle, jb, jbs, jbe, jc, jcs, jce
    REAL(wp) :: zlat

    TYPE(t_aes_phy_field),POINTER :: field => NULL()
    TYPE(t_aes_phy_tend) ,POINTER :: tend  => NULL()

      jg = p_patch%id
    
      field => prm_field(jg)
      tend  => prm_tend (jg)

      ! Inquire current grid level and the total number of grid cells
      ncd = MAX(1,p_patch%n_childdom)
      rls = grf_bdywidth_c+1
      rle = min_rlcell_int
      
      jbs     = p_patch%cells%start_blk(rls,  1)
      jbe     = p_patch%cells%  end_blk(rle,ncd)

      ! Assign initial values for some components of the "field" and
      ! "tend" state vectors.
!FIXME: PGI + OpenMP produce error in this routine... check correctness of parallel code
!$OMP PARALLEL WORKSHARE
      !
      ! constant-in-time fields
      ! initial and re-start
      !
      field%      clon(:,  :) = p_patch% cells% center(:,:)% lon
      field%      clat(:,  :) = p_patch% cells% center(:,:)% lat
      field% areacella(:,  :) = p_patch% cells%   area(:,:)
      field%    coriol(:,  :) = p_patch% cells%    f_c(:,:)
      !

!$OMP END PARALLEL WORKSHARE

      ! in case of restart, reset output fields of unused parameterizations,
      ! to their intial value
      !
      IF (isrestart()) THEN
         !
         IF ( aes_phy_tc(jg)%dt_rad == dt_zero ) THEN
            field% rld_rt      (:,:,:) = 0.0_wp
            field% rlu_rt      (:,:,:) = 0.0_wp
         END IF
         !
         IF ( aes_phy_tc(jg)%dt_two == dt_zero .AND. aes_phy_tc(jg)%dt_mig == dt_zero) THEN
            field% rsfl (:,:) = 0.0_wp
            field% ssfl (:,:) = 0.0_wp
         END IF
         !
      END IF

      ! set initial co2 flux to 0 everywhere if no restart
      IF (.NOT. isrestart()) THEN
         field% co2_flux_tile(:,:,:) = 0.0_wp
      END IF

      ! vertical diffusion
      IF (.NOT. isrestart() .AND. .NOT. aes_vdf_config(jg)%use_tmx) THEN
        IF (aes_phy_tc(jg)%dt_vdf > dt_zero) THEN
          IF (iwtr<=nsfc_type) field% z0m_tile(:,:,iwtr) = aes_vdf_config(jg)%z0m_oce
          IF (iice<=nsfc_type) field% z0m_tile(:,:,iice) = aes_vdf_config(jg)%z0m_ice
          ! IF (aes_vdf_config(jg)%use_tmx) THEN
          !   IF (iwtr<=nsfc_type) field% z0h_tile(:,:,iwtr) = aes_vdf_config(jg)%z0m_oce
          !   IF (iice<=nsfc_type) field% z0h_tile(:,:,iice) = aes_vdf_config(jg)%z0m_ice
          ! END IF
!!!          IF (iwtr<=nsfc_type) field% z0m_tile(:,:,iwtr) = aes_vdf_config(jg)%z0m_oce
!!!          IF (iice<=nsfc_type) field% z0m_tile(:,:,iice) = aes_vdf_config(jg)%z0m_ice
          IF (ilnd<=nsfc_type) THEN
            ! IF (aes_vdf_config(jg)%use_tmx) THEN
            !   field% z0m_tile(:,:,ilnd) = field%z0m(:,:) ! or maybe a larger value?
            !   field% z0h_tile(:,:,ilnd) = field%z0m(:,:) ! or maybe a larger value?
            ! ELSE
            ! IF (.NOT. aes_vdf_config(jg)%use_tmx) THEN
              field% z0m_tile(:,:,ilnd) = field%z0m(:,:) ! or maybe a larger value?
              field% z0h_lnd(:,:)       = field%z0m(:,:) ! or maybe a larger value?
            ! END IF
          END IF
        END IF
      END IF

      ! surface properties
      !
      ! Initialize some variables for water, ice and land tiles
      ! This can be overridden by the testcases below
      ! initial and re-start
      IF (iwtr <= nsfc_type) THEN
        !
        field% albvisdir_tile(:,:,iwtr) = albedoW ! albedo in the visible range for direct radiation
        field% albnirdir_tile(:,:,iwtr) = albedoW ! albedo in the NIR range for direct radiation
        field% albvisdif_tile(:,:,iwtr) = albedoW ! albedo in the visible range for diffuse radiation
        field% albnirdif_tile(:,:,iwtr) = albedoW ! albedo in the NIR range for diffuse radiation
        field% albedo_tile   (:,:,iwtr) = albedoW
        !
      END IF

      IF (ilnd <= nsfc_type) THEN
        !
        IF (.NOT. isrestart()) THEN
          IF (.NOT. ltestcase) THEN
            field%ts_tile(:,:,ilnd) = field%ts_tile(:,:,iwtr)
          END IF
        END IF
        !
        ! initial and re-start
        field% albvisdir_tile(:,:,ilnd) = field%alb(:,:)    ! albedo in the visible range for direct radiation
        field% albnirdir_tile(:,:,ilnd) = field%alb(:,:)    ! albedo in the NIR range for direct radiation
        field% albvisdif_tile(:,:,ilnd) = field%alb(:,:)    ! albedo in the visible range for diffuse radiation
        field% albnirdif_tile(:,:,ilnd) = field%alb(:,:)    ! albedo in the NIR range for diffuse radiation
        field% albedo_tile   (:,:,ilnd) = field%alb(:,:)
        !
      END IF

      IF (iice <= nsfc_type) THEN
        !
        IF (.NOT. isrestart()) THEN
          field%ts_tile(:,:,iice) = field%ts_tile(:,:,iwtr)
        END IF
        !
        ! initial and re-start
        IF(.NOT. isrestart() .OR. .NOT. aes_vdf_config(jg)%use_tmx) THEN
          field% albvisdir_tile(:,:,iice) = albi    ! albedo in the visible range for direct radiation
          field% albnirdir_tile(:,:,iice) = albi    ! albedo in the NIR range for direct radiation
          field% albvisdif_tile(:,:,iice) = albi    ! albedo in the visible range for diffuse radiation
          field% albnirdif_tile(:,:,iice) = albi    ! albedo in the NIR range for diffuse radiation
        END IF
        field% albedo_tile   (:,:,iice) = albi
        !
        IF (.NOT. isrestart()) THEN
          ! The ice model should be able to handle different thickness classes,
          ! but for AMIP we ONLY USE one ice class.
          IF (.NOT. aes_vdf_config(jg)%use_tmx) THEN
            field% albvisdir_ice(:,:,:) = albi ! albedo in the visible range for direct radiation
            field% albnirdir_ice(:,:,:) = albi ! albedo in the NIR range for direct radiation
            field% albvisdif_ice(:,:,:) = albi ! albedo in the visible range for diffuse radiation
            field% albnirdif_ice(:,:,:) = albi ! albedo in the NIR range for diffuse radiation
          END IF
          field% Tsurf(:,:,:) = Tf
          field% T1   (:,:,:) = Tf
          field% T2   (:,:,:) = Tf
          WHERE (field%seaice(:,:) > 0.0_wp)
             field% hs   (:,1,:) = 0.1_wp       ! set initial snow depth on sea ice
          ELSEWHERE
             field% hs   (:,1,:) = 0.0_wp
          ENDWHERE
          field% hi   (:,1,:) = field%siced(:,:)
          field% conc (:,1,:) = field%seaice(:,:)
        END IF
        !
      END IF

      ! For idealized test cases

      SELECT CASE (nh_test_name)

      CASE('APE','APE_aes','RCEhydro','RCE_glb','RCE_Tconst','RCE_Tprescr','aes_bubble','CBL_flxconst','RCEMIP_analytical')
        ! Note that there is only one surface type in this case !!!
        !
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          IF (.NOT. isrestart()) THEN
            DO jc = jcs,jce
              zlat = p_patch%cells%center(jc,jb)%lat
              field% ts_tile(jc,jb,iwtr) = ape_sst(ape_sst_case,zlat)
            END DO
          END IF
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% seaice(jcs:jce,jb) = 0._wp   ! zero sea ice fraction
          field% emissivity(jcs:jce,jb) = zemiss_def ! use default emissivity
          !
        END DO
!$OMP END PARALLEL DO

      CASE('aes_bubble_land')
        ! Note that there is only one surface type for land in this case !!!


        !
        IF (.NOT. aes_phy_config(jg)%ljsb) THEN 
          CALL finish('mo_aes_phy_init:init_aes_phy_field', 'aes_bubble_land testcase but JSBACH not activated (ljsb)')
        END IF

        IF (.NOT. isrestart()) THEN
          CALL jsbach_get_var('seb_t',            jg, tile='land', arr2d=field%ts_tile       (:,:,ilnd))
          IF (.NOT. aes_vdf_config(jg)%use_tmx) THEN
            CALL jsbach_get_var('turb_rough_m',     jg, tile='veg',  arr2d=field%z0m_tile    (:,:,ilnd))
          END IF
          CALL jsbach_get_var('rad_alb_vis_soil', jg, tile='veg',  arr2d=field%albvisdif_tile(:,:,ilnd))
          CALL jsbach_get_var('rad_alb_nir_soil', jg, tile='veg',  arr2d=field%albnirdif_tile(:,:,ilnd))
        END IF

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 

          IF (.NOT. isrestart()) THEN
            DO jc = jcs,jce
              field% ts            (jc,jb     ) = field%ts_tile(jc,jb,ilnd)
              IF (.NOT. aes_vdf_config(jg)%use_tmx) THEN
                field% z0m         (jc,jb     ) = field%z0m_tile(jc,jb,ilnd)
                field% z0h_lnd     (jc,jb     ) = field%z0m(jc,jb) / 3._wp
              END IF
              field% albvisdir_tile(jc,jb,ilnd) = field%albvisdif_tile(jc,jb,ilnd)
              field% albnirdir_tile(jc,jb,ilnd) = field%albnirdif_tile(jc,jb,ilnd)
              field% albvisdif     (jc,jb     ) = field%albvisdif_tile(jc,jb,ilnd)
              field% albnirdif     (jc,jb     ) = field%albnirdif_tile(jc,jb,ilnd)
              field% albvisdir     (jc,jb     ) = field%albvisdif_tile(jc,jb,ilnd)
              field% albnirdir     (jc,jb     ) = field%albnirdif_tile(jc,jb,ilnd)
              field% alb(jc,jb) = 0.5_wp * (field%albvisdif(jc,jb) + field%albnirdif(jc,jb))
            END DO
          END IF
          !
          ! initial and re-start
          field% frac_tile (jcs:jce,jb,ilnd) = 1._wp
          field% lsmask    (jcs:jce,jb     ) = 1._wp   ! land everywhere
          field% alake     (jcs:jce,jb     ) = 0._wp   ! zero lake fraction
          field% glac      (jcs:jce,jb     ) = 0._wp   ! zero glacier fraction
          field% seaice    (jcs:jce,jb     ) = 0._wp   ! zero sea ice fraction
          field% emissivity(jcs:jce,jb     ) = zemiss_def ! use default emissivity
          !
        END DO
!$OMP END PARALLEL DO





      CASE('RCE') !Note that there is only one surface type in this case
        !
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          IF (.NOT. isrestart()) THEN
            DO jc = jcs,jce
              zlat = p_patch%cells%center(jc,jb)%lat
              field% ts_tile(jc,jb,iwtr) = th_cbl(1)
            END DO
          END IF
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% seaice(jcs:jce,jb) = 0._wp   ! zeor sea ice fraction
          field% emissivity(jcs:jce,jb) = zemiss_def ! use default emissivity
          !
        END DO
!$OMP END PARALLEL DO

      CASE('APEi')
        ! The same as APE, except that whenever SST reaches tmelt, we put
        ! 1m-thick ice with a concentration of 0.9 on top
        !
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          IF (.NOT. isrestart()) THEN
            DO jc = jcs,jce
              zlat = p_patch%cells%center(jc,jb)%lat
              ! SST must reach Tf where there's ice. It may be better to modify ape_sst it self.
              field% ts_tile (jc,jb,iwtr) = ape_sst(ape_sst_case,zlat) + Tf
              field% ts_tile (jc,jb,iice) = Tf + tmelt ! K
              field% Tsurf   (jc,1, jb  ) = Tf         ! degC
              field% T1      (jc,1, jb  ) = Tf         ! degC
              field% T2      (jc,1, jb  ) = Tf         ! degC
              field% hs      (jc,1, jb  ) = 0._wp
              IF ( field% ts_tile(jc,jb,iwtr) <= Tf + tmelt ) THEN
                field% Tsurf (jc,1,jb) = field% ts_tile(jc,jb,iice) - tmelt
                field% conc  (jc,1,jb) = 0.9_wp
                field% hi    (jc,1,jb) = 1.0_wp
                field% seaice(jc,  jb) = field%conc(jc,1,jb)
              ELSE
                field% conc  (jc,1,jb) = 0._wp
                field% hi    (jc,1,jb) = 0._wp
                field% seaice(jc,  jb) = field%conc(jc,1,jb)
              ENDIF
            END DO
          END IF
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% emissivity(jcs:jce,jb) = zemiss_def ! use default emissivity
          !
        END DO
!$OMP END PARALLEL DO
        IF (.NOT. isrestart() .AND. .NOT. aes_vdf_config(jg)%use_tmx) THEN
          field% albvisdir_ice(:,:,:) = albi    ! albedo in the visible range for direct radiation
          field% albnirdir_ice(:,:,:) = albi    ! albedo in the NIR range for direct radiation
          field% albvisdif_ice(:,:,:) = albi    ! albedo in the visible range for diffuse radiation
          field% albnirdif_ice(:,:,:) = albi    ! albedo in the NIR range for diffuse radiation
        END IF
       
      CASE('APEc','APEc_nh')
        ! The same as APEi, except we initialize with no ice and don't modify the surface
        ! temperature. This is meant for a coupled run.
        !
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          IF (.NOT. isrestart()) THEN
            DO jc = jcs,jce
              zlat = p_patch%cells%center(jc,jb)%lat
              field% ts_tile (jc,jb,iwtr) = ape_sst(ape_sst_case,zlat)
              field% ts_tile (jc,jb,iice) = Tf + tmelt ! K
              field% Tsurf   (jc,1,jb)    = Tf         ! degC
              field% T1      (jc,1,jb)    = Tf         ! degC
              field% T2      (jc,1,jb)    = Tf         ! degC
              field% hs      (jc,1,jb)    = 0._wp
              field% conc    (jc,1,jb)    = 0._wp
              field% hi      (jc,1,jb)    = 0._wp
              field% seaice  (jc,  jb)    = 0._wp
            END DO
          END IF
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% emissivity(jcs:jce,jb) = zemiss_def ! use default emissivity
          !
        END DO
!$OMP END PARALLEL DO
        IF (.NOT. isrestart() .AND. .NOT. aes_vdf_config(jg)%use_tmx) THEN
          field% albvisdir_ice(:,:,:) = albi  ! albedo in the visible range for direct radiation
          field% albnirdir_ice(:,:,:) = albi  ! albedo in the NIR range for direct radiation
          field% albvisdif_ice(:,:,:) = albi  ! albedo in the visible range for diffuse radiation
          field% albnirdif_ice(:,:,:) = albi  ! albedo in the NIR range for diffuse radiation
        END IF
     
      CASE('TPEc', 'TPEo') !Note that there is only one surface type (ilnd) in this case
        !
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 1._wp   ! land fraction = 1
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% emissivity(jcs:jce,jb) = zemiss_def ! use default emissivity
          !
          IF (.NOT. isrestart()) THEN
            field% seaice(jcs:jce,jb) = 0._wp   ! zeor sea ice fraction
            field% ts_tile(jcs:jce,jb,ilnd) = tpe_temp
          END IF
          !
        END DO
!$OMP END PARALLEL DO

      CASE('JWw-Moist','LDF-Moist','jabw_m')
        !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          !
          CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
          IF (jcs>jce) CYCLE 
          !
          ! Set the surface temperature to the same value as the lowest model
          ! level above surface. For this test case, currently we assume
          ! there is no land or sea ice.
          !
          IF (.NOT. isrestart()) THEN
            field% ts_tile(jcs:jce,jb,iwtr) = temp(jcs:jce, p_patch%nlev,jb)
          END IF
          !
          ! initial and re-start
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% seaice(jcs:jce,jb) = 0._wp   ! zero sea ice fraction
          field% emissivity(jcs:jce,jb) = zemiss_def ! use default emissivity
        END DO
!$OMP END DO  NOWAIT
!$OMP END PARALLEL

      END SELECT

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
      ! initial and re-start
      DO jb = jbs,jbe
        !
        CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
        IF (jcs>jce) CYCLE 
        !
        ! Set surface tiling fractions, wrt. the cell area
        DO jc = jcs,jce
          field% sftlf (jc,jb) = field% lsmask(jc,jb) + field% alake(jc,jb) ! land incl. lakes
          field% sftgif(jc,jb) = field% lsmask(jc,jb) * field% glac (jc,jb) ! land ice
          field% sftof (jc,jb) = 1._wp - field% sftlf(jc,jb)                ! ocean
        END DO
        !
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! Settings for total surface
      ! (after tile masks and variables potentially have been overwritten by testcases above)

      ! initial and re-start
      IF (iwtr <= nsfc_type) THEN
        field%albedo   (:,:) = albedoW
      ELSE
        field%albedo   (:,:) = field%alb(:,:)
      END IF

      IF (.NOT. isrestart()) THEN
        !
        IF (iwtr <= nsfc_type) THEN
          field%ts       (:,:) = field%ts_tile(:,:,iwtr)
          field%albvisdir(:,:) = albedoW
          field%albvisdif(:,:) = albedoW
          field%albnirdir(:,:) = albedoW
          field%albnirdif(:,:) = albedoW
        ELSE IF (.NOT. (ltestcase .AND. TRIM(nh_test_name) == 'aes_bubble_land')) THEN
          field%ts       (:,:) = field%ts_tile(:,:,ilnd)
          field%albvisdir(:,:) = field%alb(:,:)
          field%albvisdif(:,:) = field%alb(:,:)
          field%albnirdir(:,:) = field%alb(:,:)
          field%albnirdif(:,:) = field%alb(:,:)
        END IF
        !
        field%ts_rad     (:,:) = field%ts(:,:)
        field%ts_rad_rt  (:,:) = field%ts(:,:)
        !
      END IF

      NULLIFY( field,tend )

  END SUBROUTINE init_aes_phy_field


  !-------------
  !>
  !! Initialize the O3 tracer from the Cariolle initial ozone field.
  !! Initialize ozone mass mixing ratios for Cariolle scheme. 
  !! An approximative initialization that considers the atmosphere as being dry is enough.
  !!
  SUBROUTINE init_o3_lcariolle( mtime_current  ,&
    &                           p_patch        ,&
    &                           pres           ,&
    &                           o3              )

    TYPE(datetime)   ,POINTER    :: mtime_current
    TYPE(t_patch)    ,INTENT(in) :: p_patch
    REAL(wp)         ,INTENT(in) :: pres          (:,:,:)
    REAL(wp)         ,INTENT(out):: o3            (:,:,:)

    ! local variables
    INTEGER  :: ncd, rls, rle, jb, jbs, jbe, jcs, jce

    ! Variables for Cariolle ozone scheme
    TYPE(t_avi)                        :: avi
    REAL(wp), TARGET                   :: latc  (nproma)
    REAL(wp), TARGET                   :: pfull (nproma, p_patch%nlev)
    REAL(wp)                           :: vmr_o3(nproma, p_patch%nlev)
    TYPE(t_time_interpolation)         :: time_interpolation
    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights

    avi%ldown=.TRUE.
    current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_current)
    time_interpolation%imonth1=current_time_interpolation_weights%month1_index
    time_interpolation%imonth2=current_time_interpolation_weights%month2_index
    time_interpolation%weight1=current_time_interpolation_weights%weight1
    time_interpolation%weight2=current_time_interpolation_weights%weight2

    ! Inquire current grid level and the total number of grid cells
    ncd = MAX(1,p_patch%n_childdom)
    rls = grf_bdywidth_c+1
    rle = min_rlcell_int
    jbs     = p_patch%cells%start_blk(rls,  1)
    jbe     = p_patch%cells%  end_blk(rle,ncd)
    !
    !$OMP PARALLEL DO PRIVATE(jb,jcs,jce,pfull,avi,latc,vmr_o3) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,jbe
      !
      CALL get_indices_c(p_patch, jb,jbs,jbe, jcs,jce, rls,rle)
      IF (jcs>jce) CYCLE 
      !
      pfull(:,:)          =  pres (:,:,jb)
      avi%pres            => pfull
      !
      latc(:)             =  p_patch% cells% center(:,jb)% lat
      avi%cell_center_lat => latc
      !
      CALL lcariolle_init_o3(jcs, jce, nproma, p_patch%nlev, time_interpolation, avi, vmr_o3)
      !
      o3(jcs:jce,:,jb) = vmr_o3(jcs:jce,:)*amo3/amd
      !
    END DO
    !$OMP END PARALLEL DO
    !
    l_cariolle_initialized_o3 = .TRUE.

  END SUBROUTINE init_o3_lcariolle

END MODULE mo_aes_phy_init
