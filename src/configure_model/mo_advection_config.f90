! configuration setup for atmospheric tracer transport
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

MODULE mo_advection_config

  USE mo_kind,                      ONLY: wp, dp
  USE mo_impl_constants,            ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom,   &
    &                                     MIURA, MIURA3, FFSL, FFSL_HYB, MCYCL,    &
    &                                     MIURA_MCYCL, MIURA3_MCYCL, FFSL_MCYCL,   &
    &                                     FFSL_HYB_MCYCL, ippm_v, ipsm_v,          &
    &                                     ino_flx, iparent_flx, inwp,              &
    &                                     iaes, SUCCESS, VNAME_LEN, NO_HADV,       &
    &                                     NO_VADV, vlname_len
  USE mo_exception,                 ONLY: message, message_text, finish
  USE mo_master_control,            ONLY: get_my_process_name
  USE mo_mpi,                       ONLY: my_process_is_stdio
  USE mo_run_config,                ONLY: msg_level
  USE mo_expression,                ONLY: expression, parse_expression_string
  USE mo_var_list,                  ONLY: t_var_list_ptr, find_list_element
  USE mo_var_list_register,         ONLY: vlr_add
  USE mo_var_list_register_utils,   ONLY: vlr_add_vref
  USE mo_var,                       ONLY: t_var
  USE mo_var_metadata_types,        ONLY: t_var_metadata, t_var_metadata_dynamic
  USE mo_var_groups,                ONLY: groups
  USE mo_var_metadata,              ONLY: get_timelevel_string
  USE mo_tracer_metadata_types,     ONLY: t_tracer_meta, t_hydro_meta
  USE mo_util_table,                ONLY: t_table, initialize_table, add_table_column, &
    &                                     set_table_entry, print_table, finalize_table

  IMPLICIT NONE

  PRIVATE

  ! types
  PUBLIC :: t_advection_config
  PUBLIC :: t_trList

  ! variables
  PUBLIC :: advection_config
  PUBLIC :: shape_func, shape_func_l
  PUBLIC :: eta, eta_l, zeta, zeta_l
  PUBLIC :: wgt_zeta, wgt_zeta_l, wgt_eta, wgt_eta_l
  PUBLIC :: lcompute, lcleanup

  ! subroutines
  PUBLIC :: configure_advection


  CHARACTER(LEN = *), PARAMETER :: modname = "mo_advection_config"


  ! Derived type to allow for the onetime computation and cleanup
  ! of tracer independent parts
  !
  TYPE t_compute
    LOGICAL :: ppm_v     (MAX_NTRACER)
    LOGICAL :: miura3_h  (MAX_NTRACER)
    LOGICAL :: ffsl_h    (MAX_NTRACER)
    LOGICAL :: ffsl_hyb_h(MAX_NTRACER)
  END TYPE t_compute


  TYPE t_scheme
    INTEGER :: iadv_min_slev     !< scheme dependent minimum vertical start level
                                 !< required for tracer-independent computations
    !INTEGER :: iadv_max_elev
    !LOGICAL :: lcompute (MAX_NTRACER)
    !LOGICAL :: lcleanup (MAX_NTRACER)
  END TYPE t_scheme



  ! Tracer ID lists
  ! Extracted from the full tracer var_list based on
  ! a specific selection criterium
  TYPE :: t_trList
    INTEGER, ALLOCATABLE :: list(:)
    INTEGER              :: len
  CONTAINS
!!$    PRIVATE
!!$    !
!!$    FINAL                :: destruct_trList
  END TYPE t_trList




  !!--------------------------------------------------------------------------
  !! Basic configuration setup for tracer advection
  !!--------------------------------------------------------------------------
  TYPE :: t_advection_config

    ! namelist variables
    CHARACTER(len=VNAME_LEN) ::  &  !< tracer-specific name suffixes
      &  tracer_names(MAX_NTRACER)  !< set by namelist, e.g. 'hus' for specific humidity
                                    !< default: 'q<tracer index>'.

    CHARACTER(len=VNAME_LEN) ::  &  !< CF convention standard names for tracer fields
      &  cfstd_names(MAX_NTRACER)   !< set internally for known tracer_names
                                    !< default: 'tracer_<tracer_names(tracer index)>'.

    CHARACTER(len=VNAME_LEN) ::  &  !< long names for tracer fields
      &  long_names(MAX_NTRACER)    !< set internally for known tracer_names
                                    !< default: 'tracer_<tracer_names(tracer index)>'.

    INTEGER :: nname                !< number of names read from transport_nml/tracer_names
                                    !< which are stored in advection_config/tracer_names

    INTEGER :: &                    !< selects horizontal transport scheme
      &  ihadv_tracer(MAX_NTRACER)  !< 0:  no horizontal advection
                                    !< 2:  2nd order miura
                                    !< 3:  3rd order miura with quadr./cubic reconstr.
                                    !< 4:  Flux form semi lagrange (FFSL)
                                    !< 5:  hybrid FFSL Miura3
                                    !< 6:  3rd or 4th order upstream (on hexagons only)
                                    !< 20: subcycling version of miura
                                    !< 22: 2nd order miura and miura_cycl
                                    !< 32: 3rd order miura with miura_cycl
                                    !< 42: FFSL with miura_cycl
                                    !< 52: FFSL_HYB with miura_cycl

    INTEGER :: &                    !< selects vertical transport scheme
      &  ivadv_tracer(MAX_NTRACER)  !< 0 : no vertical advection
                                    !< 1 : 1st order upwind
                                    !< 2 : 3rd order PSM for CFL>1
                                    !< 3 : 3rd order PPM for CFL>1

    INTEGER :: &                    !< advection of TKE
      &  iadv_tke                   !< 0 : none
                                    !< 1 : vertical advection only
                                    !< 2 : vertical and horizontal advection

    LOGICAL :: lvadv_tracer         !< if .TRUE., calculate vertical tracer advection
    LOGICAL :: lclip_tracer         !< if .TRUE., clip negative tracer values

    LOGICAL :: llsq_svd             !< least squares reconstruction with
                                    !< singular value decomposition (TRUE) or
                                    !< QR decomposition (FALSE) of design matrix A
    INTEGER :: &                    !< parameter used to select the limiter
      &  itype_vlimit(MAX_NTRACER)  !< for vertical transport

    INTEGER :: &                    !< parameter used to select the limiter
      &  itype_hlimit(MAX_NTRACER)  !< for horizontal transport
                                    !< 0: no limiter
                                    !< 3: monotonous flux limiter
                                    !< 4: positive definite flux limiter

    INTEGER :: &                    !< additional method for identifying and avoiding
      & ivlimit_selective(MAX_NTRACER)!< spurious limiting of smooth extrema
                                    !< 1: switch on
                                    !< 0: switch off

    REAL(wp):: beta_fct             !< global boost factor for range of permissible values in
                                    !< (semi-) monotonous flux limiter. A value larger than
                                    !< 1 allows for (small) over and undershoots, while a value
                                    !< of 1 gives strict monotonicity (at the price of increased
                                    !< diffusivity).

    INTEGER :: igrad_c_miura        !< parameter used to select the gradient
                                    !< reconstruction method at cell center
                                    !< for second order miura scheme

    INTEGER :: ivcfl_max            !< determines stability range of vertical
                                    !< ppm-scheme (approximate allowable maximum
                                    !< CFL-number)

    INTEGER :: npassive_tracer      !< number of additional passive tracers, in addition to
                                    !< microphysical- and ART tracers.

    INTEGER :: nadv_substeps        !< number of substeps per fast physics time step
                                    !< for the Miura-type substepping schemes 20, 22, 32, 42, 52

    CHARACTER(len=MAX_CHAR_LENGTH) :: &!< Comma separated list of initialization formulae
      &  init_formula                  !< for passive tracers.


    ! derived variables

    INTEGER  :: iubc_adv         !< selects upper boundary condition
                                 !< for tracer transport
                                 !< 0: no flux
                                 !< 1: zero gradient
                                 !< 2: interpolated flux from parent grid

    INTEGER ::  &                !< selects vertical start level for each patch
      &  iadv_slev(MAX_NTRACER)  !< and each tracer.

    INTEGER :: kstart_aero(2), & !< start and end levels for vertical flux averaging
               kend_aero(2)      !< for advection of 2D (climatological) aerosol fields

    INTEGER ::  &                !< vertical end level down to which qv is
      &  iadv_qvsubstep_elev     !< advected with internal substepping (to
                                 !< circumvent CFL instability in the
                                 !< stratopause region).

    TYPE(t_trList) ::       &    !< tracer sublist containing all tracer fields of
      &  trHydroMass             !< type hydroMass.
    TYPE(t_trList) ::       &    !< tracer sublist containing all tracer fields
      &  trAdvect                !< which are advected.
    TYPE(t_trList) ::       &    !< tracer sublist containing all tracer fields
      &  trNotAdvect             !< which are not advected.
    TYPE(t_trList) ::       &    !< tracer sublist containing all tracer fields
      &  trFeedback              !< for which child-to-parent feedback is applied.

    LOGICAL :: isAnyTypeMiura    !< TRUE if any tracer is to be advected with MIURA scheme
    LOGICAL :: isAnyTypeMcycl    !< TRUE if any tracer is to be advected with subcycling

    ! scheme specific derived variables
    !
    TYPE(t_scheme) :: ppm_v      !< vertical PPM scheme
    TYPE(t_scheme) :: miura_h    !< horizontal miura scheme (linear reconstr.)
    TYPE(t_scheme) :: miura3_h   !< horizontal miura scheme (higher order reconstr.)
    TYPE(t_scheme) :: ffsl_h     !< horizontal FFSL scheme
    TYPE(t_scheme) :: ffsl_hyb_h !< horizontal hybrid FFSL/miura3 scheme
    TYPE(t_scheme) :: mcycl_h    !< horizontal miura scheme (linear) with subcycling


  CONTAINS
    PRIVATE
    PROCEDURE :: print_setup   => advection_print_setup
  END TYPE t_advection_config

  !>
  !!
  TYPE(t_advection_config), TARGET :: advection_config(0:max_dom)

  TYPE(t_compute)  :: lcompute
  TYPE(t_compute)  :: lcleanup


  ! for first order Gauss-Legendre quadrature
  !
  REAL(wp) :: shape_func_l(4)  !< shape functions for mapping the FFSL departure
                               !< region onto the standard rectangle

  REAL(wp) :: zeta_l, eta_l    !< Gauss quadrature point in \zeta-\eta space

  REAL(wp) :: wgt_zeta_l       !< Gauss quadrature weights for zeta and eta
  REAL(wp) :: wgt_eta_l        !< points


  ! for second order Gauss-Legendre quadrature
  !
  REAL(wp) :: shape_func(4,4)  !< shape functions for mapping the FFSL departure
                               !< region onto the standard rectangle (miura3 only)

  REAL(wp) :: zeta(4), eta(4)  !< Gauss quadrature points in \zeta-\eta space
                               !< (miura3 only)

  REAL(wp) :: wgt_zeta(4)      !< Gauss quadrature weights for zeta and eta
  REAL(wp) :: wgt_eta(4)       !< points (miura3 only)


CONTAINS

  !>
  !! setup components of the transport scheme depending on this namelist
  !!
  !! Setup of additional transport control variables depending on the
  !! transport-NAMELIST and potentially other namelists. This routine is
  !! called, after all namelists have been read and a synoptic consistency
  !! check has been done.
  !!
  SUBROUTINE configure_advection( jg, nlev, nlev_1, iforcing, iqc, iqt,                &
    &                             kstart_moist, kend_qvsubstep, lvert_nest,            &
    &                             ntracer, prog_list, tracer_list, kstart_tracer )
    !
    INTEGER,              INTENT(IN)   :: jg               !< patch
    INTEGER,              INTENT(IN)   :: nlev             !< number of vertical levels
    INTEGER,              INTENT(IN)   :: nlev_1           !< vertical levels of global patch
    INTEGER,              INTENT(IN)   :: iforcing
    INTEGER,              INTENT(IN)   :: iqc, iqt         !< hydrometeor indices
    INTEGER,              INTENT(IN)   :: kstart_moist
    INTEGER,              INTENT(IN)   :: kend_qvsubstep
    INTEGER,              INTENT(IN)   :: ntracer          !< total number of tracers
    LOGICAL,              INTENT(IN)   :: lvert_nest
    TYPE(t_var_list_ptr), INTENT(IN)   :: prog_list(:)     !< list of all prognostic variables
    TYPE(t_var_list_ptr), INTENT(INOUT):: tracer_list(:)   !< list of tracers
    INTEGER,    OPTIONAL, INTENT(IN)   :: kstart_tracer(:) !< start index for (art-)tracer related processes

    !
    CHARACTER(*), PARAMETER :: routine = modname//"::configure_advection"
    INTEGER :: jt                        !< tracer loop index
    INTEGER :: jm                        !< loop index for shape functions

    INTEGER, POINTER :: ivadv_tracer(:)  !< convenience pointer
    INTEGER, POINTER :: ihadv_tracer(:)  !< convenience pointer
    INTEGER, PARAMETER :: n_timelevels = 2
    INTEGER, PARAMETER :: itime = 1      !< tracer_list time level
                                         !< here it does not matter if we use 1 or 2
    INTEGER :: z_go_tri(10)              !< for crosscheck
    CHARACTER(len=vlname_len) :: listname


    ! Build tracer list from the prognostic state for time levels `now` and `new`
    !
    ! Saftey check
    IF ( (UBOUND(prog_list,1) < n_timelevels) .OR. (UBOUND(tracer_list,1) < n_timelevels) ) THEN
      WRITE(message_text,'(a,i2,a)') 'Upper bound of prog_list or tracer_list below n_timelevels=', &
        &                        n_timelevels, '. Please check allocation.'
      CALL finish(routine, message_text)
    ENDIF
    !
    DO jt = 1, n_timelevels
      ! Build prog state tracer list
      ! no memory allocation (only references to prog list)
      !
      WRITE(listname,'(a,i2.2,a,i2.2)') 'nh_state_tracer_of_domain_',jg, &
        &                               '_and_timelev_',jt
      CALL new_nh_state_tracer_list(jg, prog_list(jt), tracer_list(jt), listname )
    ENDDO ! jt


    !--------------------------------------------------------------------
    ! Consistency checks
    !--------------------------------------------------------------------

    ! Flux computation methods - consistency check
    z_go_tri(1:10)=(/NO_HADV,MIURA,MIURA3,FFSL,FFSL_HYB,MCYCL,       &
      &              MIURA_MCYCL,MIURA3_MCYCL,FFSL_MCYCL,FFSL_HYB_MCYCL/)
    DO jt=1,ntracer
      IF ( ALL(z_go_tri /= advection_config(jg)%ihadv_tracer(jt)) ) THEN
        CALL finish( routine,                                       &
          &  'incorrect settings for TRI-C grid ihadv_tracer. Must be '// &
          &  '0,2,3,4,5,6,20,22,32,42 or 52 ')
      ENDIF
    ENDDO

    !-----------------------------------------------------------------------

    !
    ! set transport variables/model components, which depend on
    ! the transport namelist and potentially other namelsists.
    !

    !
    ! set vertical start level for each patch and each tracer
    !
    advection_config(jg)%iadv_slev(:) = 1
    advection_config(jg)%iadv_qvsubstep_elev = 1
    IF (iforcing == inwp .OR. iforcing == iaes) THEN
      ! Set iadv_slev to kstart_moist for all moisture fields but QV
      ! note: iqt denotes the first tracer index not related to moisture
      advection_config(jg)%iadv_slev(iqc:iqt-1) = kstart_moist
      advection_config(jg)%iadv_qvsubstep_elev = kend_qvsubstep
      IF ( PRESENT(kstart_tracer) ) THEN
         advection_config(jg)%iadv_slev(iqt:ntracer) = kstart_tracer(iqt:ntracer)
      ENDIF
    ENDIF


    ! set boundary condition for vertical transport
    !
    IF (.NOT. lvert_nest ) THEN ! no vertical nesting
      advection_config(jg)%iubc_adv = ino_flx    ! no flux ubc
    ELSE ! vertical nesting
      IF (nlev < nlev_1) THEN
        advection_config(jg)%iubc_adv = iparent_flx
      ELSE
        advection_config(jg)%iubc_adv = ino_flx
      ENDIF
    ENDIF


    ! dummy initialization of index fields for transport of 2D aerosol fields
    DO jt = 1, 2
      advection_config(:)%kstart_aero(jt) = 0
      advection_config(:)%kend_aero(jt) = 0
    ENDDO

    ! to save some paperwork
    ivadv_tracer => advection_config(1)%ivadv_tracer(:)
    ihadv_tracer => advection_config(1)%ihadv_tracer(:)


    ! PPM_V specific settings (vertical transport)
    !
    lcompute%ppm_v(:)   = .FALSE.
    lcleanup%ppm_v(:)   = .FALSE.

    advection_config(jg)%ppm_v%iadv_min_slev = HUGE(1)


    IF ( ANY(ivadv_tracer == ippm_v) .OR. ANY(ivadv_tracer == ipsm_v) ) THEN

      ! compute minimum required slev for this group of tracers
      DO jt=1,ntracer
        IF ( ANY( (/ippm_v, ipsm_v/) == ivadv_tracer(jt) ) ) THEN
          advection_config(jg)%ppm_v%iadv_min_slev =                           &
            &                  MIN( advection_config(jg)%ppm_v%iadv_min_slev,  &
            &                        advection_config(jg)%iadv_slev(jt) )
        ENDIF
      ENDDO

      ! Search for the first tracer jt for which vertical advection of
      ! type PPM/PSM has been selected.
      DO jt=1,ntracer
        IF ( ANY( (/ippm_v, ipsm_v/) == ivadv_tracer(jt) ) ) THEN
          lcompute%ppm_v(jt) = .TRUE.
          exit
        ENDIF
      ENDDO

      ! Search for the last tracer jt for which vertical advection of
      ! type PPM/PSM has been selected.
      DO jt=ntracer,1,-1
        IF ( ANY( (/ippm_v, ipsm_v/) == ivadv_tracer(jt) ) ) THEN
          lcleanup%ppm_v(jt) = .TRUE.
          exit
        ENDIF
      ENDDO
    END IF


    !
    ! MIURA specific settings (horizontal transport)
    !

    advection_config(jg)%miura_h%iadv_min_slev = HUGE(1)


    ! compute minimum required slev for this group of tracers
    DO jt=1,ntracer
      IF ( ANY( (/MIURA, MIURA_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        advection_config(jg)%miura_h%iadv_min_slev =  &
          &                  MIN( advection_config(jg)%miura_h%iadv_min_slev, &
          &                       advection_config(jg)%iadv_slev(jt) )
      ENDIF
    ENDDO


    !
    ! MIURA3 specific settings (horizontal transport)
    !
    lcompute%miura3_h(:) = .FALSE.
    lcleanup%miura3_h(:) = .FALSE.

    advection_config(jg)%miura3_h%iadv_min_slev = HUGE(1)

    ! compute minimum required slev for this group of tracers
    DO jt=1,ntracer
      IF ( ANY( (/MIURA3, MIURA3_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        advection_config(jg)%miura3_h%iadv_min_slev =   &
          &                  MIN( advection_config(jg)%miura3_h%iadv_min_slev, &
          &                       advection_config(jg)%iadv_slev(jt) )
      ENDIF
    ENDDO

    ! Search for the first tracer jt for which horizontal advection of
    ! type MIURA3 has been selected.
    DO jt=1,ntracer
      IF ( ANY( (/MIURA3, MIURA3_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcompute%miura3_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO

    ! Search for the last tracer jt for which horizontal advection of
    ! type MIURA3 has been selected.
    DO jt=ntracer,1,-1
      IF ( ANY( (/MIURA3, MIURA3_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcleanup%miura3_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO


    !
    ! FFSL specific settings (horizontal transport)
    !
    lcompute%ffsl_h(:) = .FALSE.
    lcleanup%ffsl_h(:) = .FALSE.

    advection_config(jg)%ffsl_h%iadv_min_slev = HUGE(1)

    ! compute minimum required slev for this group of tracers
    DO jt=1,ntracer
      IF ( ANY( (/FFSL, FFSL_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        advection_config(jg)%ffsl_h%iadv_min_slev =   &
          &                  MIN( advection_config(jg)%ffsl_h%iadv_min_slev, &
          &                       advection_config(jg)%iadv_slev(jt) )
      ENDIF
    ENDDO

    ! Search for the first tracer jt for which horizontal advection of
    ! type FFSL has been selected.
    DO jt=1,ntracer
      IF ( ANY( (/FFSL, FFSL_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcompute%ffsl_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO

    ! Search for the last tracer jt for which horizontal advection of
    ! type FFSL has been selected.
    DO jt=ntracer,1,-1
      IF ( ANY( (/FFSL, FFSL_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcleanup%ffsl_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO



    !
    ! FFSL_HYB specific settings (horizontal transport)
    !
    lcompute%ffsl_hyb_h(:) = .FALSE.
    lcleanup%ffsl_hyb_h(:) = .FALSE.

    advection_config(jg)%ffsl_hyb_h%iadv_min_slev = HUGE(1)

    ! compute minimum required slev for this group of tracers
    DO jt=1,ntracer
      IF ( ANY( (/FFSL_HYB, FFSL_HYB_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        advection_config(jg)%ffsl_hyb_h%iadv_min_slev =   &
          &                  MIN( advection_config(jg)%ffsl_hyb_h%iadv_min_slev, &
          &                       advection_config(jg)%iadv_slev(jt) )
      ENDIF
    ENDDO

    ! Search for the first tracer jt for which horizontal advection of
    ! type FFSL_HYB has been selected.
    DO jt=1,ntracer
      IF ( ANY( (/FFSL_HYB, FFSL_HYB_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcompute%ffsl_hyb_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO

    ! Search for the last tracer jt for which horizontal advection of
    ! type FFSL_HYB has been selected.
    DO jt=ntracer,1,-1
      IF ( ANY( (/FFSL_HYB, FFSL_HYB_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcleanup%ffsl_hyb_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO


    !
    ! MCYCL specific settings (horizontal transport)
    !

    advection_config(jg)%mcycl_h%iadv_min_slev = HUGE(1)

    ! compute minimum required slev for this group of tracers
    DO jt=1,ntracer
      IF ( ANY( (/MCYCL, MIURA_MCYCL, MIURA3_MCYCL, FFSL_MCYCL, FFSL_HYB_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        advection_config(jg)%mcycl_h%iadv_min_slev =   &
          &                  MIN( advection_config(jg)%mcycl_h%iadv_min_slev, &
          &                       advection_config(jg)%iadv_slev(jt) )
      ENDIF
    ENDDO


    ! Check, whether any of the tracers is supposed to be transported horizontally
    ! with a scheme of Miura type (i.e. 2nd order scheme with linear reconstruction)
    advection_config(jg)%isAnyTypeMiura=.FALSE.
    DO jt=1,ntracer
      IF (ANY((/MIURA, MIURA_MCYCL/) == ihadv_tracer(jt))) THEN
        advection_config(jg)%isAnyTypeMiura=.TRUE.
        exit
      ENDIF
    ENDDO
    !
    ! Check, whether any of the tracers is supposed to be transported horizontally
    ! with substepping (i.e. 2nd order scheme with substepping)
    advection_config(jg)%isAnyTypeMcycl=.FALSE.
    DO jt=1,ntracer
      IF (ANY((/MCYCL, MIURA_MCYCL, MIURA3_MCYCL, FFSL_MCYCL, FFSL_HYB_MCYCL/) == ihadv_tracer(jt))) THEN
        advection_config(jg)%isAnyTypeMcycl=.TRUE.
        exit
      ENDIF
    ENDDO


    !
    ! Compute shape functions for mapping the departure region onto the
    ! standard rectangle. Integration points, shape functions and quadrature
    ! weights are provided for a first and second order Gauss-Legendre
    ! quadrature.
    !

    IF (jg == 1) THEN

      !
      ! First order
      !

      ! Coordinates of integration points (in \zeta,\eta-System)
      !
      zeta_l = 0._wp
      eta_l  = 0._wp

      ! shape function for mapping
      shape_func_l(1) = 0.25_wp * (1._wp-zeta_l)*(1._wp-eta_l)
      shape_func_l(2) = 0.25_wp * (1._wp+zeta_l)*(1._wp-eta_l)
      shape_func_l(3) = 0.25_wp * (1._wp+zeta_l)*(1._wp+eta_l)
      shape_func_l(4) = 0.25_wp * (1._wp-zeta_l)*(1._wp+eta_l)
      !$ACC ENTER DATA COPYIN(shape_func_l)


      ! Gauss quadrature weights
      !
      wgt_zeta_l = 2._wp
      wgt_eta_l  = 2._wp


      !
      ! Second order
      !

      ! Coordinates of integration points (in \zeta,\eta-System)
      !
      zeta(1) = -1._wp/SQRT(3._wp)
      zeta(2) =  1._wp/SQRT(3._wp)
      zeta(3) =  1._wp/SQRT(3._wp)
      zeta(4) = -1._wp/SQRT(3._wp)

      eta(1)  = -1._wp/SQRT(3._wp)
      eta(2)  = -1._wp/SQRT(3._wp)
      eta(3)  =  1._wp/SQRT(3._wp)
      eta(4)  =  1._wp/SQRT(3._wp)

      ! shape functions for mapping
      !
      DO jm = 1,4
        shape_func(1,jm) = 0.25_wp * (1._wp-zeta(jm))*(1._wp-eta(jm))
        shape_func(2,jm) = 0.25_wp * (1._wp+zeta(jm))*(1._wp-eta(jm))
        shape_func(3,jm) = 0.25_wp * (1._wp+zeta(jm))*(1._wp+eta(jm))
        shape_func(4,jm) = 0.25_wp * (1._wp-zeta(jm))*(1._wp+eta(jm))
      END DO
      !$ACC ENTER DATA COPYIN(shape_func)

      ! Gauss quadrature weights
      !
      wgt_zeta(1) = 1._wp
      wgt_zeta(2) = 1._wp
      wgt_zeta(3) = 1._wp
      wgt_zeta(4) = 1._wp

      wgt_eta(1)  = 1._wp
      wgt_eta(2)  = 1._wp
      wgt_eta(3)  = 1._wp
      wgt_eta(4)  = 1._wp
    END IF


    !******************************************
    ! Tracer-Sublist extraction
    !******************************************

    ! create list of tracers which are advected
    ! list is allowed to have zero size.
    !
    advection_config(jg)%trAdvect = subListExtract(from_list       = tracer_list(itime), &
      &                                            extraction_rule = extraction_rule_advect)


    ! create list of tracers which are not advected
    ! list is allowed to have zero size.
    !
    advection_config(jg)%trNotAdvect = subListExtract(from_list       = tracer_list(itime),         &
      &                                               extraction_rule = extraction_rule_notAdvect)


    ! create ID list for tracer group hydroMass
    ! This list contains the IDs of all condensate fields
    ! which are required for computing the water loading term.
    !
    IF ( iforcing == inwp .OR. iforcing == iaes ) THEN
      advection_config(jg)%trHydroMass = subListExtract(from_list       = tracer_list(itime),         &
        &                                               extraction_rule = extraction_rule_hydroMass)
      !
      ! in these cases (NWP, AES) an empty group `hydroMass` is not allowed.
      IF ( advection_config(jg)%trHydroMass%len < 1 ) THEN
        CALL finish (routine, 'trHydroMass%list is empty. At least one condensate field is required.')
      ENDIF

    ENDIF


    ! create list of tracers for which child-to-parent feedback is applied
    ! list is allowed to have zero size.
    !
    advection_config(jg)%trFeedback = subListExtract(from_list       = tracer_list(itime),         &
      &                                              extraction_rule = extraction_rule_feedback)


    ! initialize passive tracers, if required
    !
    IF (advection_config(jg)%npassive_tracer > 0) THEN
      CALL init_passive_tracer (tracer_list, advection_config(jg), ntl=1)
    ENDIF

    ! print setup
    IF (msg_level >= 10) THEN
      IF(my_process_is_stdio()) THEN
        CALL advection_config(jg)%print_setup(tracer_list(itime),nlev)
      ENDIF
    ENDIF

  END SUBROUTINE configure_advection



  !-----------------------------------------------------------------------------
  !>
  !! Extract a sublist from var_list. The sublist will only contain the
  !! meta information <ncontained> from the info-state. This routine can be used,
  !! e.g. for creating a tracer sublist from the full tracer list. The
  !! extraction-rule(s) must be passed in terms of a function, which returns
  !! -999 in case that the field does not match the extraction-rule(s) and
  !! <ncontained> otherwise.
  !
  TYPE(t_trList) FUNCTION subListExtract (from_list, extraction_rule) RESULT(obj)
    TYPE(t_var_list_ptr), INTENT(IN) :: from_list         !< variable list (metadata)
    INTERFACE
      INTEGER FUNCTION extraction_rule(info, tracer_info) RESULT(id)
        IMPORT                            :: t_var_metadata, t_tracer_meta
        TYPE (t_var_metadata), INTENT(IN) :: info             ! static info state
        CLASS(t_tracer_meta) , INTENT(IN) :: tracer_info      ! dynamic (tracer) info state
      END FUNCTION extraction_rule
    END INTERFACE
    ! local vars
    CHARACTER(*), PARAMETER :: routine = modname//"::subListExtract"
    TYPE(t_var), POINTER :: element
    TYPE(t_var_metadata) , POINTER :: info             ! static info state
    CLASS(t_tracer_meta) , POINTER :: tracer_info      ! dynamic (tracer) info state
    INTEGER, ALLOCATABLE :: tmp(:)                     ! temporary array
    INTEGER :: ist, id, i

    ! allocate list with maximum size
    ALLOCATE(obj%list(from_list%p%nvars), stat=ist)
    IF(ist/=SUCCESS) CALL finish (TRIM(routine), 'alloc of obj%list failed')
    ! initialize
    obj%len = 0
    ! Sub-list extraction (IDs only)
    DO i = 1, from_list%p%nvars
      element => from_list%p%vl(i)%p
      ! retrieve information from actual linked list element
      info          => element%info
      tracer_info   => element%info_dyn%tracer
      ! extract sublist member
      id = extraction_rule(info, tracer_info)
      IF (id /= -999) THEN
        obj%len           = obj%len+1
        obj%list(obj%len) = id
      ENDIF
    ENDDO
    ! contract list
    ALLOCATE(tmp(obj%len), stat=ist)
    IF(ist/=SUCCESS) CALL finish (TRIM(routine), 'alloc of tmp failed')
    tmp(1:obj%len) = obj%list(1:obj%len)
    CALL MOVE_ALLOC(tmp,obj%list)

    ! sanity check
    IF (.NOT. ALLOCATED(obj%list)) THEN
      CALL finish (routine, 'tracer subList is not ALLOCATED')
    ENDIF
  END FUNCTION subListExtract

  !-----------------------------------------------------------------------------
  !>
  !! If the tracer at hand is a member of the hydroMass ID
  !! list, this function gives back its respective ID.
  !! Otherwise, it gives back -999
  !!
  INTEGER FUNCTION extraction_rule_hydroMass(info,tracer_info) RESULT(id)
    !
    TYPE(t_var_metadata), INTENT(IN)  :: info
    CLASS(t_tracer_meta), INTENT(IN)  :: tracer_info

    SELECT TYPE(tracer_info)
    TYPE IS (t_hydro_meta)
      id = info%ncontained
    CLASS default
      id = -999
    !
    END SELECT
  END FUNCTION extraction_rule_hydroMass


  !-----------------------------------------------------------------------------
  !>
  !! If the tracer at hand is advected (either horizontally or vertically)
  !! this function gives back its respective ID.
  !! Otherwise, it gives back -999
  !!
  INTEGER FUNCTION extraction_rule_advect(info,tracer_info) RESULT(id)
    !
    TYPE(t_var_metadata), INTENT(IN)  :: info
    CLASS(t_tracer_meta), INTENT(IN)  :: tracer_info
    !
    IF ((tracer_info%ihadv_tracer/=NO_HADV) .OR. (tracer_info%ivadv_tracer/=NO_VADV)) THEN
      id = info%ncontained
    ELSE
      id = -999
    ENDIF
  END FUNCTION extraction_rule_advect


  !-----------------------------------------------------------------------------
  !>
  !! If the tracer at hand is not advected (neither horizontally nor vertically)
  !! this function gives back its respective ID.
  !! Otherwise, it gives back -999
  !!
  INTEGER FUNCTION extraction_rule_notAdvect(info,tracer_info) RESULT(id)
    !
    TYPE(t_var_metadata), INTENT(IN)  :: info
    CLASS(t_tracer_meta), INTENT(IN)  :: tracer_info
    !
    IF ((tracer_info%ihadv_tracer==NO_HADV) .AND. (tracer_info%ivadv_tracer==NO_VADV)) THEN
      id = info%ncontained
    ELSE
      id = -999
    ENDIF
  END FUNCTION extraction_rule_notAdvect


  !-----------------------------------------------------------------------------
  !>
  !! If child-to-parent feedback should be applied
  !! this function gives back its respective ID.
  !! Otherwise, it gives back -999
  !!
  INTEGER FUNCTION extraction_rule_feedback(info,tracer_info) RESULT(id)
    !
    TYPE(t_var_metadata), INTENT(IN)  :: info
    CLASS(t_tracer_meta), INTENT(IN)  :: tracer_info
    !
    IF ( (tracer_info%lfeedback) ) THEN
      id = info%ncontained
    ELSE
      id = -999
    ENDIF
  END FUNCTION extraction_rule_feedback


  ! ATTENTION: currently not used (see FINAL statement above)
  !
  !-------------------------------------------------------------------------
  !>
  !! Deallocate object components
  !!
  !! Deallocates all components of a class t_trList object
  !!
  SUBROUTINE destruct_trList(obj)
    TYPE(t_trList) :: obj
    !
    ! local
    INTEGER :: ist
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_config: destruct_trList'

    IF (ALLOCATED(obj%list)) THEN
      DEALLOCATE(obj%list, STAT=ist)
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine), 'list deallocation failed' )
      ELSE
        write(0,*) "deallocated obj%list"
      ENDIF
    ENDIF
  END SUBROUTINE destruct_trList



  !>
  !! Initialize passive tracers
  !!
  !! Additional passive tracers are initialized by applying
  !! the initialization formulae provided via Namelist parameter
  !! 'init_formula'.
  !!
  SUBROUTINE init_passive_tracer (tracer_list, advection_config, ntl)

    TYPE(t_var_list_ptr)        , INTENT(IN) :: tracer_list(:)
    TYPE(t_advection_config), INTENT(IN) :: advection_config ! config state
    INTEGER                 , INTENT(IN) :: ntl              ! time level

    ! local variables
    !
    INTEGER :: ipassive                           ! Loop counter
    INTEGER :: start_pos
    INTEGER :: end_pos
    INTEGER :: pos
    TYPE(expression) :: formula
    CHARACTER(LEN=4) :: passive_tracer_id         ! tracer ID string
    CHARACTER(LEN=4) :: str_ntl                   ! time level string
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: tracer_name ! tracer name string

    !-------------------------------------------------------------------------------

    ! init
    start_pos= 1
    end_pos  = 0

    ! loop over all additional passive tracers
    !
    DO ipassive=1,advection_config%npassive_tracer
      pos = INDEX(advection_config%init_formula(start_pos:), ",")
      IF (pos == 0) THEN
        end_pos = MAX_CHAR_LENGTH
      ELSE
        end_pos = end_pos + pos
      ENDIF
      CALL parse_expression_string(formula, &
           advection_config%init_formula(start_pos:end_pos-1))
      ! generate tracer name
      WRITE(passive_tracer_id,'(I2)') ipassive
      str_ntl = get_timelevel_string(ntl)
      tracer_name = 'Qpassive_'//TRIM(ADJUSTL(passive_tracer_id))//TRIM(str_ntl)


      WRITE(message_text,'(2a)') 'Initialize additional passive tracer: ',TRIM(tracer_name)
      CALL message('',message_text)
      !NOTE (HB): if wp /= dp the following is not correct, since r_ptr is of type REAL(dp)
      CALL formula%evaluate( fget_var_list_element_r3d (tracer_list(ntl), &
        &                    TRIM(tracer_name)))
      CALL formula%finalize()


      start_pos=end_pos+1

    ENDDO

  CONTAINS
    FUNCTION fget_var_list_element_r3d (this_list, vname) RESULT(ptr)
      TYPE(t_var_list_ptr), INTENT(in) :: this_list    ! list
      CHARACTER(*), INTENT(in) :: vname         ! name of variable
      REAL(dp), POINTER    :: ptr(:,:,:)   ! reference to allocated field
      TYPE(t_var), POINTER :: element

      element => find_list_element(this_list, vname)
      NULLIFY (ptr)
      IF (element%info%lcontained) THEN
        IF (ASSOCIATED(element)) ptr => element%r_ptr(:,:,:,element%info%ncontained,1)
      ELSE
        IF (ASSOCIATED(element)) ptr => element%r_ptr(:,:,:,1,1)
      ENDIF
    END FUNCTION fget_var_list_element_r3d

  END SUBROUTINE init_passive_tracer



  !>
  !! Screen print out of advection setup
  !!
  SUBROUTINE advection_print_setup (config_obj, var_list_tracer, nlev)
    !
    CLASS(t_advection_config)             :: config_obj        !< object for which the setup will be printed
    TYPE(t_var_list_ptr)     , INTENT(IN) :: var_list_tracer   !< variable list (metadata)
    INTEGER                  , INTENT(IN) :: nlev              !< number of vertical levels
    ! local variables
    TYPE(t_var_metadata), POINTER :: info
    CLASS(t_tracer_meta), POINTER :: tracer_info
    TYPE(t_table) :: table
    INTEGER       :: irow, tracer_id, i
    INTEGER       :: slev
    CHARACTER(LEN=3) :: str_tracer_id
    CHARACTER(LEN=3) :: str_startlev
    CHARACTER(LEN=7) :: str_substep_range
    CHARACTER(LEN=3) :: str_nadv_substeps
    CHARACTER(LEN=3) :: str_flag

    ! could this be transformed into a table header?
    write(0,*) "Tracer meta-information for patch ", var_list_tracer%p%patch_id

    ! table-based output
    CALL initialize_table(table)
    ! the latter is no longer mandatory
    CALL add_table_column(table, "VarName")
    CALL add_table_column(table, "Tracer ID")
    CALL add_table_column(table, "feedback")
    CALL add_table_column(table, "in list trAdvect")
    CALL add_table_column(table, "in list trNotAdvect")
    CALL add_table_column(table, "in list trHydroMass")
    CALL add_table_column(table, "slev")
    CALL add_table_column(table, "substep range")
    CALL add_table_column(table, "nadv_substeps")

    irow = 0
    ! print tracer meta-information
    DO i = 1, var_list_tracer%p%nvars
      info => var_list_tracer%p%vl(i)%p%info
      tracer_info => var_list_tracer%p%vl(i)%p%info_dyn%tracer

      tracer_id = info%ncontained
      irow = irow + 1
      !
      CALL set_table_entry(table,irow,"VarName", TRIM(tracer_info%name))
      !
      write(str_tracer_id,'(i3)')  tracer_id
      CALL set_table_entry(table,irow,"Tracer ID", str_tracer_id)
      !
      str_flag = MERGE('X',' ',ANY(config_obj%trFeedback%list==tracer_id))
      CALL set_table_entry(table,irow,"feedback", TRIM(str_flag))
      !
      str_flag = MERGE('X',' ',ANY(config_obj%trAdvect%list==tracer_id))
      CALL set_table_entry(table,irow,"in list trAdvect", TRIM(str_flag))
      !
      str_flag = MERGE('X',' ',ANY(config_obj%trNotAdvect%list==tracer_id))
      CALL set_table_entry(table,irow,"in list trNotAdvect", TRIM(str_flag))
      !
      ! iforcing == inwp/iaes
      IF (ALLOCATED(config_obj%trHydroMass%list)) THEN
        str_flag = MERGE('X',' ',ANY(config_obj%trHydroMass%list==tracer_id))
        CALL set_table_entry(table,irow,"in list trHydroMass", TRIM(str_flag))
      ENDIF
      !
      ! print start level for transport and
      ! range of levels for which substepping is applied
      IF (ANY(config_obj%trAdvect%list==tracer_id)) THEN
        !
        slev = config_obj%iadv_slev(tracer_id)
        write(str_startlev,'(i3)') slev
        !
        IF (ANY((/MIURA_MCYCL, MIURA3_MCYCL, FFSL_MCYCL, FFSL_HYB_MCYCL/) &
           &     == config_obj%ihadv_tracer(tracer_id))) THEN
          write(str_substep_range,'(i3,a,i3)')  slev,'/',config_obj%iadv_qvsubstep_elev
          write(str_nadv_substeps,'(i3)')  config_obj%nadv_substeps
        ELSE IF (config_obj%ihadv_tracer(tracer_id) == MCYCL) THEN
          write(str_substep_range,'(i3,a,i3)')  slev,'/',nlev
          write(str_nadv_substeps,'(i3)')  config_obj%nadv_substeps
        ELSE
          write(str_substep_range,'(a)') '-- / --'
          str_nadv_substeps = '--'
        ENDIF
      ELSE
        !
        write(str_startlev,'(a)') '--'
        write(str_substep_range,'(a)') '-- / --'
        write(str_nadv_substeps,'(a)') '--'
      ENDIF
      CALL set_table_entry(table,irow,"slev", TRIM(str_startlev))
      CALL set_table_entry(table,irow,"substep range", TRIM(str_substep_range))
      CALL set_table_entry(table,irow,"nadv_substeps", TRIM(str_nadv_substeps))
    ENDDO

    CALL print_table(table, opt_delimiter=' | ')
    CALL finalize_table(table)

    WRITE (0,*) " " ! newline
  END SUBROUTINE advection_print_setup


  !-------------------------------------------------------------------------
  !> Creates tracer var list.
  !!
  !! Creates tracer var list containing references to all prognostic tracer
  !! fields.
  !!
  SUBROUTINE new_nh_state_tracer_list (patch_id, from_var_list, p_tracer_list, listname)
    INTEGER, INTENT(IN) :: patch_id ! current patch ID
    TYPE(t_var_list_ptr), INTENT(IN) :: from_var_list ! source list to be referenced
    TYPE(t_var_list_ptr), INTENT(INOUT) :: p_tracer_list ! new tracer list (containing all tracers)
    CHARACTER(*), INTENT(IN) :: listname
    TYPE (t_var_metadata), POINTER :: from_info
    TYPE (t_var_metadata_dynamic), POINTER :: from_info_dyn
    INTEGER :: iv

    ! Register a field list and apply default settings
    CALL vlr_add(p_tracer_list, TRIM(listname), patch_id=patch_id, &
      &          lrestart=.FALSE., loutput =.FALSE., model_type=get_my_process_name())
    ! add references to all tracer fields of the source list (prognostic state)
    DO iv = 1, from_var_list%p%nvars
      ! retrieve information from actual linked list element
      from_info => from_var_list%p%vl(iv)%p%info
      from_info_dyn => from_var_list%p%vl(iv)%p%info_dyn
      ! Only add tracer fields to the tracer list
      IF (from_info_dyn%tracer%lis_tracer .AND. .NOT.from_info%lcontainer) &
        & CALL vlr_add_vref(p_tracer_list, from_info%name, from_var_list)
    END DO
  END SUBROUTINE new_nh_state_tracer_list

END MODULE mo_advection_config
