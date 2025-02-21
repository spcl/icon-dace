! Contains the setup of configuration of the
! nonhydrostatic dynamical core
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

MODULE mo_nonhydrostatic_nml

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_impl_constants,        ONLY: max_dom, TRACER_ONLY
  USE mo_io_units,              ONLY: nnml, nnml_output
  USE mo_namelist,              ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,        ONLY: use_restart_namelists
  USE mo_mpi,                   ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att,   ONLY: open_tmpfile, store_and_close_namelist,  &
                                    & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,          ONLY: temp_defaults, temp_settings

  USE mo_nonhydrostatic_config, ONLY: &
                                    ! from namelist
                                    & config_itime_scheme     => itime_scheme     , &
                                    & config_ndyn_substeps    => ndyn_substeps    , &
                                    & config_vcfl_threshold   => vcfl_threshold   , &
                                    & config_lextra_diffu     => lextra_diffu     , &
                                    & config_divdamp_fac      => divdamp_fac      , &
                                    & config_divdamp_fac2     => divdamp_fac2     , &
                                    & config_divdamp_fac3     => divdamp_fac3     , &
                                    & config_divdamp_fac4     => divdamp_fac4     , &
                                    & config_divdamp_z        => divdamp_z        , &
                                    & config_divdamp_z2       => divdamp_z2       , &
                                    & config_divdamp_z3       => divdamp_z3       , &
                                    & config_divdamp_z4       => divdamp_z4       , &
                                    & config_divdamp_fac_o2   => divdamp_fac_o2   , &
                                    & config_divdamp_order    => divdamp_order    , &
                                    & config_divdamp_type     => divdamp_type     , &
                                    & config_divdamp_trans_start => divdamp_trans_start, &
                                    & config_divdamp_trans_end => divdamp_trans_end, &
                                    & config_ivctype          => ivctype          , &
                                    & config_htop_moist_proc  => htop_moist_proc  , &
                                    & config_hbot_qvsubstep   => hbot_qvsubstep   , &
                                    & config_htop_aero_proc   => htop_aero_proc   , &
                                    & config_damp_height      => damp_height      , &
                                    & config_rayleigh_type    => rayleigh_type    , &
                                    & config_rayleigh_coeff   => rayleigh_coeff   , &
                                    & config_vwind_offctr     => vwind_offctr     , &
                                    & config_rhotheta_offctr  => rhotheta_offctr  , &
                                    & config_veladv_offctr    => veladv_offctr    , &
                                    & config_iadv_rhotheta    => iadv_rhotheta    , &
                                    & config_igradp_method    => igradp_method    , &
                                    & config_exner_expol      => exner_expol      , &
                                    & config_l_zdiffu_t       => l_zdiffu_t       , &
                                    & config_thslp_zdiffu     => thslp_zdiffu     , &
                                    & config_thhgtd_zdiffu    => thhgtd_zdiffu


  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: read_nonhydrostatic_namelist


CONTAINS
  !-------------------------------------------------------------------------
  !! Read Namelist for nonhydrostatic core. 
  !!
  !! This subroutine 
  !! - reads the Namelist for nonhydrostatic core
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  SUBROUTINE read_nonhydrostatic_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg           ! loop index
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nonhydrostatic_nml:read_nonhydrostatic_namelist'

    !-----------------------------------------------------------------------------
    ! Namelist variables
    !-----------------------------------------------------------------------------
    INTEGER  :: itime_scheme   ! parameter used to select the time stepping scheme
    ! = 1, explicit 2 time level scheme for tracer
    ! 4, 5, 6 = predictor-corrector scheme
    ! 4: Contravariant vertical velocity is computed in the predictor step only,
    !    velocity tendencies are computed in the corrector step only (most efficient option)
    ! 5: Contravariant vertical velocity is computed in both substeps (beneficial for numerical
    !    stability in very-high resolution setups with extremely steep slops, otherwise no significant impact)
    ! 6: As 5, but velocity tendencies are also computed in both substeps (no apparent benefit, but more expensive)

    INTEGER :: ndyn_substeps           ! number of dynamics substeps per fast-physics step
    REAL(wp):: vcfl_threshold          ! threshold for vertical advection CFL number at which the adaptive time step reduction
                                       ! (increase of ndyn_substeps w.r.t. the fixed fast-physics time step) is triggered
    LOGICAL :: lextra_diffu            ! if true: apply additional diffusion at grid points close
    ! to the CFL stability limit for vertical advection
    REAL(wp):: divdamp_fac             ! Scaling factor for divergence damping at height divdamp_z and below
    REAL(wp):: divdamp_fac2            ! Scaling factor for divergence damping at height divdamp_z2
    REAL(wp):: divdamp_fac3            ! Scaling factor for divergence damping at height divdamp_z3
    REAL(wp):: divdamp_fac4            ! Scaling factor for divergence damping at height divdamp_z4 and above
    REAL(wp):: divdamp_z               ! Height up to which divdamp_fac is used, start of linear profile
    REAL(wp):: divdamp_z2              ! Height of divdamp_fac2, end of linear and start of quadratic profile
    REAL(wp):: divdamp_z3              ! Height of divdamp_fac3, to define quadratic profile
    REAL(wp):: divdamp_z4              ! Height from which divdamp_fac4, end of quadratic profile
    INTEGER :: divdamp_order           ! Order of divergence damping
    INTEGER :: divdamp_type            ! Type of divergence damping (2D or 3D divergence)
    REAL(wp):: divdamp_trans_start     ! Lower bound of transition zone between 2D and 3D div damping in case of divdamp_type = 32
    REAL(wp):: divdamp_trans_end       ! Upper bound of transition zone between 2D and 3D div damping in case of divdamp_type = 32
    INTEGER :: ivctype                 ! Type of vertical coordinate (Gal-Chen / SLEVE)
    REAL(wp):: htop_moist_proc         ! Top height (in m) of the part of the model domain
                                       ! where processes related to moist physics are computed
    REAL(wp):: hbot_qvsubstep          ! Bottom height (in m) down to which water vapor is
                                       ! advected with internal substepping (to circumvent CFL
                                       ! instability in the stratopause region).
    REAL(wp):: htop_aero_proc          ! Top height (in m) of the part of the model domain
                                       ! where processes related to (ART) aerosol tracers are computed

    INTEGER :: rayleigh_type           ! type of Rayleigh damping (1: CLASSIC, 2: Klemp (2008))
    REAL(wp):: damp_height(max_dom)    ! height at which w-damping and sponge layer start
    REAL(wp):: rayleigh_coeff(max_dom) ! Rayleigh damping coefficient in w-equation
    REAL(wp):: vwind_offctr            ! Off-centering in vertical wind solver
    REAL(wp):: rhotheta_offctr         ! Off-centering for density and potential temperature at interface levels
    REAL(wp):: veladv_offctr           ! Off-centering for velocity advection
    INTEGER :: iadv_rhotheta           ! Advection scheme used for density and pot. temperature
    INTEGER :: igradp_method           ! Method for computing the horizontal presure gradient
    REAL(wp):: exner_expol             ! Temporal extrapolation of Exner for computation of
                                       ! horizontal pressure gradient
    LOGICAL :: l_zdiffu_t              ! .true.: apply truly horizontal temperature diffusion
    !         over steep slopes
    REAL(wp):: thslp_zdiffu            ! threshold slope above which temperature diffusion is applied
    REAL(wp):: thhgtd_zdiffu           ! threshold height diff. between adjacent model grid points
    ! above which temperature diffusion is applied


    NAMELIST /nonhydrostatic_nml/ itime_scheme, ndyn_substeps, ivctype, htop_moist_proc,    &
         & hbot_qvsubstep, damp_height, rayleigh_type,               &
         & rayleigh_coeff, vwind_offctr, iadv_rhotheta,              &
         & divdamp_fac, divdamp_fac2, divdamp_fac3, divdamp_fac4,    &
         & divdamp_z, divdamp_z2, divdamp_z3, divdamp_z4,            &
         & igradp_method, exner_expol, l_zdiffu_t,                   &
         & thslp_zdiffu, thhgtd_zdiffu, divdamp_order, divdamp_type, &
         & rhotheta_offctr, lextra_diffu, veladv_offctr,             &
         & divdamp_trans_start, divdamp_trans_end, htop_aero_proc,   &
         & vcfl_threshold

    !-----------------------
    ! 1. default settings
    !-----------------------

    ! Time scheme chose for the nh-model
    itime_scheme = 4 ! Predictor-corrector scheme with averaged velocity tendency
                     ! Velocity tendency is recomputed for predictor step only after physics calls

    ! number of dynamics substeps per fast-physics timestep
    ndyn_substeps = 5

    ! threshold for vertical advection CFL number at which the adaptive time step reduction
    ! (increase of ndyn_substeps w.r.t. the fixed fast-physics time step) is triggered
    vcfl_threshold = 1.05_wp

    ! apply additional horizontal diffusion on vn and w at grid points close to the stability
    ! limit for vertical advection
    lextra_diffu = .TRUE.

    ! scaling factor for divergence damping
    divdamp_fac  = 0.0025_wp
    divdamp_fac2 = 0.0040_wp
    divdamp_fac3 = 0.0040_wp
    divdamp_fac4 = 0.0040_wp
    divdamp_z    = 32500._wp
    divdamp_z2   = 40000._wp
    divdamp_z3   = 60000._wp
    divdamp_z4   = 80000._wp

    ! Order of divergence damping
    divdamp_order = 4

    ! Type of divergence damping
    divdamp_type = 3

    ! Lower and upper bound of transition zone between 2D and 3D divergence damping in case of divdamp_type = 32
    divdamp_trans_start = 12500._wp   ! 12.5 km
    divdamp_trans_end   = 17500._wp   ! 17.5 km

    ! Type of vertical coordinate (1: Gal-Chen, 2: SLEVE)
    ivctype  = 2

    ! Turn off moist physics above 22.5 km
    htop_moist_proc = 22500._wp

    !use half the transport time step above 24 km to ensure CFL stability
    !     (requires choosing ihadv_tracer(1) = 22, 32 or 42!)
    hbot_qvsubstep  = 22500._wp

    ! Turn off art aerosol tracer related processes above 22.5 km (i.e. per default like moist physics)
    htop_aero_proc = 22500._wp

    ! type of Rayleigh damping
    rayleigh_type     = 2           ! Klemp-type Rayleigh damping
    ! Rayleigh damping of w above 45 km
    damp_height(1)    = 45000.0_wp
    ! Corresponding damping coefficient
    rayleigh_coeff(1) = 0.05_wp
    ! Off-centering of vertical wind speed in vertically implicit solver
    ! When combining coarse spatial resolutions (R2B5 or coarser) with high model tops (> 50 km),
    ! this value can be increased up to 1.0 in order to stabilize the numerical treatment of sound waves
    vwind_offctr      = 0.15_wp
    ! Off-centering for density and potential temperature at interface levels
    ! Specifying a negative value here reduces the amount of vertical wind off-centering needed for
    ! stability of sound waves. 
    rhotheta_offctr   = -0.1_wp
    ! Off-centering of velocity advection in corrector step
    veladv_offctr = 0.25_wp
    ! Use Miura scheme for advection of rho and theta
    iadv_rhotheta     = 2
    ! Use truly horizontal pressure-gradient computation to ensure numerical stability
    ! without heavy orography smoothing
    igradp_method     = 3
    ! Extrapolate Exner function by 1/3 time step for computing the horizontal pressure gradient
    ! Tests indicate that for coarse resolutions (R2B5 or coarser), optimal stability is reached
    ! for values between 1/2 and 2/3, whereas for high resolutions, where stability limitations
    ! arise from large-amplitude breaking gravity waves rather than sound wave reflections, values
    ! around 1/3 are better.



    exner_expol       = 1._wp/3._wp


    ! dummy values for nested domains; will be reset to value of domain 1 
    ! if not specified explicitly in the namelist
    damp_height(2:max_dom)    = -1.0_wp
    rayleigh_coeff(2:max_dom) = -1.0_wp

    ! truly horizontal temperature diffusion
    l_zdiffu_t     = .TRUE.   ! turned on by default
    thslp_zdiffu   = 0.025_wp ! slope threshold 0.025
    thhgtd_zdiffu  = 200._wp  ! threshold for height difference between adjacent grid points 200 m

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('nonhydrostatic_nml')
      READ(funit,NML=nonhydrostatic_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('nonhydrostatic_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, nonhydrostatic_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, nonhydrostatic_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, nonhydrostatic_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !--------------------------------------------------------------------
    ! Sanity and internal checks
    !--------------------------------------------------------------------
    ! Reset values of rayleigh_coeff and damp_height to that of the 
    ! global domain if not specified
    SELECT CASE (itime_scheme)
    CASE (TRACER_ONLY, 3,4,5,6 ) !OK
    CASE DEFAULT
      CALL finish( TRIM(routine),'wrong value for nonhydrostatic_nml:itime_scheme. '//&
                 & 'See mo_impl_constants.f90 for possible options.' )
    END SELECT

    DO jg = 2, max_dom
      IF (damp_height(jg) < 0.0_wp) THEN
        damp_height(jg) = damp_height(1)
      ENDIF
      IF (rayleigh_coeff(jg) < 0.0_wp) THEN
        rayleigh_coeff(jg) = rayleigh_coeff(1)
      ENDIF
    ENDDO

    IF ( ndyn_substeps <= 0) THEN
        CALL finish( TRIM(routine), 'Invalid number of dynamics substeps '//&
          &'Value must be positive')
    ENDIF

    ! For backward compatibility with the previous implementation of divergence damping flow control
    IF (divdamp_order == 5) THEN
      divdamp_order = 4
      divdamp_type  = 3
      CALL message(TRIM(routine), 'WARNING: divdamp_order = 5 has been reset to 4')
    ENDIF

    SELECT CASE (divdamp_order)
    CASE (2,4, 24) !OK
    CASE DEFAULT
      CALL finish( TRIM(routine),'Invalid value for divdamp_order (must be 2, 24 or 4)' )
    END SELECT

    ! Checks for divdamp parameters
    IF ( divdamp_fac  < 0.0_wp ) THEN
       CALL finish( TRIM(routine), 'divdamp_fac  < 0 not allowed')
    ENDIF
    IF ( divdamp_fac2 < 0.0_wp ) THEN
       CALL finish( TRIM(routine), 'divdamp_fac2 < 0 not allowed')
    ENDIF
    IF ( divdamp_fac3 < 0.0_wp ) THEN
       CALL finish( TRIM(routine), 'divdamp_fac3 < 0 not allowed')
    ENDIF
    IF ( divdamp_fac4 < 0.0_wp ) THEN
       CALL finish( TRIM(routine), 'divdamp_fac4 < 0 not allowed')
    ENDIF
    IF ( divdamp_z2 <= divdamp_z  ) THEN
       CALL finish( TRIM(routine), 'divdamp_z2 <= divdamp_z  not allowed')
    ENDIF
    IF ( divdamp_z2 >= divdamp_z4 ) THEN
       CALL finish( TRIM(routine), 'divdamp_z2 >= divdamp_z4 not allowed')
    ENDIF
    IF ( divdamp_z3 == divdamp_z2 ) THEN
       CALL finish( TRIM(routine), 'divdamp_z3 == divdamp_z2 not allowed')
    ENDIF
    IF ( divdamp_z3 == divdamp_z4 ) THEN
       CALL finish( TRIM(routine), 'divdamp_z3 == divdamp_z4 not allowed')
    ENDIF


    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

       config_rayleigh_type     = rayleigh_type
       config_rayleigh_coeff(:) = rayleigh_coeff(:)
       config_damp_height   (:) = damp_height   (:)
       config_iadv_rhotheta     = iadv_rhotheta
       config_vwind_offctr      = vwind_offctr
       config_rhotheta_offctr   = rhotheta_offctr
       config_veladv_offctr     = veladv_offctr
       config_igradp_method     = igradp_method
       config_exner_expol       = exner_expol
       config_ndyn_substeps     = ndyn_substeps
       config_vcfl_threshold    = vcfl_threshold
       config_lextra_diffu      = lextra_diffu
       config_divdamp_fac       = divdamp_fac
       config_divdamp_fac2      = divdamp_fac2
       config_divdamp_fac3      = divdamp_fac3
       config_divdamp_fac4      = divdamp_fac4
       config_divdamp_z         = divdamp_z
       config_divdamp_z2        = divdamp_z2
       config_divdamp_z3        = divdamp_z3
       config_divdamp_z4        = divdamp_z4
       config_divdamp_fac_o2    = divdamp_fac ! initialization - divdamp_fac_o2 is a derived variable that may change during runtime
       config_divdamp_order     = divdamp_order
       config_divdamp_type      = divdamp_type
       config_divdamp_trans_start = divdamp_trans_start
       config_divdamp_trans_end   = divdamp_trans_end
       config_itime_scheme      = itime_scheme
       config_ivctype           = ivctype
       config_l_zdiffu_t        = l_zdiffu_t
       config_thslp_zdiffu      = thslp_zdiffu
       config_thhgtd_zdiffu     = thhgtd_zdiffu
       config_htop_moist_proc   = htop_moist_proc
       config_hbot_qvsubstep    = hbot_qvsubstep
       config_htop_aero_proc    = htop_aero_proc

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=nonhydrostatic_nml)
      CALL store_and_close_namelist(funit, 'nonhydrostatic_nml')
    ENDIF
    ! 6. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nonhydrostatic_nml)

  END SUBROUTINE read_nonhydrostatic_namelist


END MODULE mo_nonhydrostatic_nml
