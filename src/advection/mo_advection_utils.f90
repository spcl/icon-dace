! Some utilities which are specific to the transport algorithm.
!
! Module contains some functions and procedures which are specifically related
! to the transport schemes. These subroutines or functions are needed at
! various places within the transport scheme. Therefore outsourcing these
! routines protects from possible circular dependencies.
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
#include "omp_definitions.inc"
!----------------------------
MODULE mo_advection_utils

  USE mo_kind,                  ONLY: wp
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH, inwp, ICOSMO,     &
    &                                 IPROG, VNAME_LEN,                  &
    &                                 MAX_CHAR_LENGTH
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_fortran_tools,         ONLY: t_ptr_2d3d
  USE mo_cf_convention,         ONLY: t_cf_var
  USE mo_grib2,                 ONLY: t_grib2_var
  USE mo_var_list, ONLY: add_ref, find_list_element, t_var_list_ptr
  USE mo_var, ONLY: t_var
  USE mo_var_metadata_types,    ONLY: t_var_metadata,                    &
    &                                 t_vert_interp_meta,                &
    &                                 t_hor_interp_meta,                 &
    &                                 t_post_op_meta
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta
  USE mo_var_groups,            ONLY: MAX_GROUPS
  USE mo_advection_config,      ONLY: t_advection_config
  USE mo_comin_config,          ONLY: comin_config
#ifndef __NO_ICON_COMIN__
  USE comin_host_interface,     ONLY: comin_request_get_list_head,       &
    &                                 t_var_request_list_item
#endif

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_advection_utils"

  ! functions
  PUBLIC :: laxfr_upflux
  PUBLIC :: laxfr_upflux_v

  ! subroutines
  PUBLIC :: add_tracer_ref            ! add new tracer component
  PUBLIC :: init_tracer_settings

  ! types
  PUBLIC :: t_list2D

  INTERFACE add_tracer_ref
    MODULE PROCEDURE add_var_list_reference_tracer
  END INTERFACE add_tracer_ref

  TYPE t_list2D
    INTEGER, POINTER :: eidx(:,:)
    INTEGER, POINTER :: elev(:,:)
    INTEGER, POINTER :: len(:)
    INTEGER          :: npoints
  END TYPE t_list2D


CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux,.
  !!
  !! Lax Friedrichs first order upwind flux,
  !! used in conservative advection routines.
  !! For passive advection, equivalent to
  !! any other first order upwind flux.
  !!
  FUNCTION laxfr_upflux( p_vn, p_psi1, p_psi2 )  RESULT(p_upflux)
    !

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: p_vn
    REAL(wp), INTENT(in) :: p_psi1, p_psi2

    REAL(wp) :: p_upflux

    !$ACC ROUTINE SEQ

    !-----------------------------------------------------------------------
    p_upflux = 0.5_wp * (        p_vn  *( p_psi1 + p_psi2 )    &
      &                   - ABS( p_vn )*( p_psi2 - p_psi1 ) )

  END FUNCTION laxfr_upflux



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux for vertical advection,.
  !!
  !! Generalized Lax Friedrichs first order upwind flux,
  !! used in conservative vertical advection routines.
  !! For passive advection, equivalent to any other first
  !! order upwind flux.
  !! Applicable to height based vertical coordinate systems. 
  !!
  FUNCTION laxfr_upflux_v( p_w, p_psi1, p_psi2 )  RESULT(p_upflux)
    !

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: p_w
    REAL(wp), INTENT(in) :: p_psi1, p_psi2

    REAL(wp) :: p_upflux

    !$ACC ROUTINE SEQ

    !-----------------------------------------------------------------------
    p_upflux = 0.5_wp * (        p_w  *( p_psi1 + p_psi2 )    &
      &                   + ABS( p_w )*( p_psi2 - p_psi1 ) )

  END FUNCTION laxfr_upflux_v


  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! reference to an existing pointer to a 3D tracer field
  ! optionally overwrite some default meta data
  !
  SUBROUTINE add_var_list_reference_tracer(this_list, target_name, tracer_name,    &
    &        tracer_idx, ptr_arr, cf, grib2, advconf, ldims, loutput, lrestart,    &
    &        isteptype, tlev_source, vert_interp, hor_interp, in_group, post_op,   &
    &        tracer_info, new_element)

    TYPE(t_var_list_ptr)    , INTENT(inout)        :: this_list
    CHARACTER(len=*)    , INTENT(in)           :: target_name
    CHARACTER(len=*)    , INTENT(in)           :: tracer_name
    INTEGER             , INTENT(out)          :: tracer_idx       ! index in 4D tracer container
    TYPE(t_ptr_2d3d)    , INTENT(inout)        :: ptr_arr(:)
    TYPE(t_cf_var)      , INTENT(in)           :: cf               ! CF related metadata
    TYPE(t_grib2_var)   , INTENT(in)           :: grib2            ! GRIB2 related metadata
    TYPE(t_advection_config), INTENT(inout)    :: advconf          ! adv configure state
    INTEGER             , INTENT(in), OPTIONAL :: ldims(3)         ! local dimensions, for checking
    LOGICAL             , INTENT(in), OPTIONAL :: loutput          ! output flag
    LOGICAL             , INTENT(in), OPTIONAL :: lrestart         ! restart flag
    INTEGER,              INTENT(in), OPTIONAL :: isteptype        ! type of statistical processing
    INTEGER             , INTENT(in), OPTIONAL :: tlev_source      ! actual TL for TL dependent vars
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp   ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp    ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(MAX_GROUPS)          ! groups to which a variable belongs
    TYPE(t_post_op_meta), INTENT(in), OPTIONAL :: post_op          ! post operation (e.g. scale with const. factor or rho)
    CLASS(t_tracer_meta), INTENT(in), OPTIONAL :: tracer_info      ! tracer meta data
    TYPE(t_var), POINTER, OPTIONAL :: new_element

    ! Local variables:
    TYPE(t_var), POINTER :: target_var
    TYPE(t_var_metadata), POINTER :: target_info

    INTEGER :: zihadv_tracer, zivadv_tracer

    CHARACTER(*), PARAMETER :: routine = "add_tracer_ref"

  !------------------------------------------------------------------

    ! get pointer to target element (in this case 4D tracer container)
    target_var => find_list_element (this_list, target_name)
    ! get tracer field metadata
    target_info => target_var%info

    ! get index of current field in 4D container and set
    ! tracer index accordingly.
    ! Note that this will be repeated for each patch. Otherwise it may happen that
    ! some MPI-processes miss this assignment.
    !
    tracer_idx = target_info%ncontained+1  ! index in 4D tracer container

    WRITE (message_text,'(a,i3,a,a)')                                            &
      & "tracer index ", tracer_idx," assigned to tracer field ", TRIM(tracer_name)
    CALL message(TRIM(routine),message_text)



    ! set default values
    ! If ihadv_tracer or ivadv_tracer are not present, take respective value from the
    ! advection_config state.
    ! If ihadv_tracer or ivadv_tracer are present, take those values, and overwrite
    ! the respective values of the advection_config state.
    !
    IF( PRESENT(tracer_info) ) THEN
      IF ( tracer_info%ihadv_tracer /= 2) THEN
        zihadv_tracer = tracer_info%ihadv_tracer
        ! BE AWARE, THAT ihadv_tracer IS NOT SANITY-CHECKED. THIS OVERWRITES THE
        ! SANITY CHECKED NAMELIST SETTINGS.
        advconf%ihadv_tracer(tracer_idx) = tracer_info%ihadv_tracer
      ELSE
        zihadv_tracer = advconf%ihadv_tracer(tracer_idx)
      ENDIF
    ENDIF

    IF( PRESENT(tracer_info) ) THEN
      IF ( tracer_info%ivadv_tracer /= 3) THEN
        zivadv_tracer = tracer_info%ivadv_tracer
        ! BE AWARE, THAT ivadv_tracer IS NOT SANITY-CHECKED. THIS OVERWRITES THE
        ! SANITY CHECKED NAMELIST SETTINGS.
        advconf%ivadv_tracer(tracer_idx) = tracer_info%ivadv_tracer
      ELSE
        zivadv_tracer = advconf%ivadv_tracer(tracer_idx)
      ENDIF
    ENDIF

    ! create new table entry reference including additional tracer metadata
    CALL add_ref( this_list, target_name, tracer_name, ptr_arr(tracer_idx)%p_3d, &
       &          target_info%hgrid, target_info%vgrid, cf, grib2,               &
       &          ref_idx=tracer_idx,                                            &
       &          ldims=ldims, loutput=loutput, lrestart=lrestart,               &
       &          isteptype=isteptype, tlev_source=tlev_source,                  &
       &          vert_interp=vert_interp, hor_interp=hor_interp,                &
       &          tracer_info=tracer_info, in_group=in_group, post_op=post_op,   &
       &          new_element=new_element)

  END SUBROUTINE add_var_list_reference_tracer


  !> Determine the number (and names) of tracer variables which are
  !  later assigned automatically in mo_nonhydro_state.
  !
  SUBROUTINE init_tracer_settings(iforcing, n_dom, ltransport,                 &
    &                             inwp_turb, inwp_gscp, inwp_convection,       &
    &                             lart, iart_ntracer, ctracer_art,             &
    &                             advection_config,                            &
    &                             iqv, iqc, iqi, iqr, iqs, iqt, iqg, iqni,     &
    &                             iqh, iqnr, iqns, iqng, iqnh, iqnc,           &
    &                             iqgl, iqhl, inccn, ininact, ininpot,         &
    &                             iqtke, iqm_max, ntracer, nqtendphy, nclass_gscp, &
    &                             iqbin, iqb_i, iqb_e, iqb_s)
    INTEGER,                  INTENT(IN)    :: iforcing, n_dom
    LOGICAL,                  INTENT(IN)    :: ltransport
    INTEGER,                  INTENT(IN)    :: inwp_turb(:)
    INTEGER,                  INTENT(IN)    :: inwp_gscp(:)
    INTEGER,                  INTENT(IN)    :: inwp_convection(:)
    LOGICAL,                  INTENT(IN)    :: lart
    INTEGER,                  INTENT(IN)    :: iart_ntracer
    CHARACTER(len=MAX_CHAR_LENGTH), ALLOCATABLE, INTENT(IN) :: ctracer_art(:)
    TYPE(t_advection_config), INTENT(INOUT) :: advection_config(:)
    INTEGER,                  INTENT(INOUT) :: iqv, iqc, iqi, iqr, iqs, iqt, iqg, &
      &                                        iqni, iqh, iqnr, iqns, iqng, iqnh, &
      &                                        iqnc, iqgl, iqhl, inccn,           &
      &                                        ininact, ininpot, iqtke,           &
      &                                        iqb_i, iqb_e, iqb_s
    INTEGER, DIMENSION(:),    INTENT(INOUT) :: iqbin
    INTEGER,                  INTENT(INOUT) :: iqm_max
    INTEGER,                  INTENT(INOUT) :: ntracer
    INTEGER,                  INTENT(INOUT) :: nqtendphy
    INTEGER,                  INTENT(INOUT) :: nclass_gscp(:)
    !
    CHARACTER(len=*), PARAMETER :: routine =  modname//'::init_tracer_settings'
    INTEGER  :: jg, name_len, itracer
    CHARACTER(len=MAX_CHAR_LENGTH) :: src_name_str
#ifndef __NO_ICON_COMIN__
    TYPE(t_var_request_list_item), POINTER :: ptr
#endif

    INTEGER  :: iqb
    CHARACTER(len=7) :: qbinname(SIZE(iqbin))
    qbinname=['qbin001','qbin002','qbin003','qbin004','qbin005',&
           &  'qbin006','qbin007','qbin008','qbin009','qbin010',&
           &  'qbin011','qbin012','qbin013','qbin014','qbin015',&
           &  'qbin016','qbin017','qbin018','qbin019','qbin020',&
           &  'qbin021','qbin022','qbin023','qbin024','qbin025',&
           &  'qbin026','qbin027','qbin028','qbin029','qbin030',&
           &  'qbin031','qbin032','qbin033',                    &
           &  'qbin034','qbin035','qbin036','qbin037','qbin038',&
           &  'qbin039','qbin040','qbin041','qbin042','qbin043',&
           &  'qbin044','qbin045','qbin046','qbin047','qbin048',&
           &  'qbin049','qbin050','qbin051','qbin052','qbin053',&
           &  'qbin054','qbin055','qbin056','qbin057','qbin058',&
           &  'qbin059','qbin060','qbin061','qbin062','qbin063',&
           &  'qbin064','qbin065','qbin066']

    ! Check settings of ntracer
    !
    ! provisional number of water species for which convective
    ! and turbulent tendencies of NWP physics are stored
    nqtendphy = 0

    SELECT CASE(iforcing)
    CASE (INWP) ! iforcing

      ! ** NWP physics section **
      !
      ! IMPORTANT: For NWP physics, five microphysics tracers (QV, QC,
      !            QI, QR and QS) must always be defined because their
      !            presence is implicitly assumed in the
      !            physics-dynamics interface when converting between
      !            temperature and virtual temperature.
      !
      !            Any additional mass-related microphysics tracers
      !            must be numbered in consecutive order after iqs =
      !            5. The parameter "iqm_max" must signify the highest
      !            tracer index carrying a moisture mixing
      !            ratio. Additional tracers for hydrometeor number
      !            concentrations (in case of a two-moment scheme) or
      !            other purposes (aerosols, qt variance or anything
      !            else) can be numbered freely after iqm_max. The
      !            parameter "iqt", denoting the start index of
      !            tracers not related at all to moisture, is used in
      !            configure_advection to specify the index range of
      !            tracers for which advection is turned off in the
      !            stratosphere (i.e. all cloud and precipitation
      !            variables including number concentrations)
      !
      !            The indices for specific tracers (iqv, iqc, iqi,
      !            ...) have initial values 0, as valid for unused
      !            tracers. Below the indices are properly defined for
      !            the tracers to be used for the selected physics
      !            configuration. Accordingly also the names for
      !            theses tracers are defined.
      !
      !            Note also that the namelist parameter "ntracer" is
      !            reset automatically to the correct value when NWP
      !            physics is used in order to avoid multiple namelist
      !            changes when playing around with different physics
      !            schemes.
      !
      ! Default settings valid for all microphysics options
      !
      iqv       = 1 ; advection_config(:)%tracer_names(iqv) = 'qv' !> water vapour
      iqc       = 2 ; advection_config(:)%tracer_names(iqc) = 'qc' !! cloud water
      iqi       = 3 ; advection_config(:)%tracer_names(iqi) = 'qi' !! ice
      iqr       = 4 ; advection_config(:)%tracer_names(iqr) = 'qr' !! rain water
      iqs       = 5 ; advection_config(:)%tracer_names(iqs) = 'qs' !! snow
      !! number of water species for which convective and turbulent
      !! tendencies are stored
      nqtendphy = 3
      !
      ! The following parameters may be reset depending on the
      ! selected physics scheme
      !
      iqm_max   = 5     !! end index of water species mass mixing ratios
      iqt       = 6     !! start index of other tracers not related at all to moisture
      !
      ntracer   = 5     !! total number of tracers

      ! Taking the 'highest' microphysics option in some cases allows
      ! using more complex microphysics schemes in nested domains than
      ! in the global domain.
      ! However, a clean implementation would require a
      ! domain-dependent 'ntracer' dimension
      SELECT CASE (MAXVAL(inwp_gscp(1:n_dom)))


      CASE(2)  ! COSMO-DE (3-cat ice: snow, cloud ice, graupel)

        ! CALL finish('mo_atm_nml_crosscheck', 'Graupel scheme not implemented.')

        iqg     = 6 ; advection_config(:)%tracer_names(iqg) = 'qg' !! graupel
        iqm_max = iqg
        iqt     = iqt + 1

        ntracer = ntracer + 1  !! increase total number of tracers by 1


      CASE(3)  ! improved ice nucleation scheme C. Koehler (note:
               ! iqm_max does not change!)

        iqni     = 6 ; advection_config(:)%tracer_names(iqni)     = 'qni'     !! cloud ice number
        ininact  = 7 ; advection_config(:)%tracer_names(ininact)  = 'ninact'  !! activated ice nuclei
        iqt      = iqt + 2

        ntracer = ntracer + 2  !! increase total number of tracers by 2

      CASE(4)  ! two-moment scheme

        iqg     = 6  ; advection_config(:)%tracer_names(iqg)     = 'qg'
        iqh     = 7  ; advection_config(:)%tracer_names(iqh)     = 'qh'
        iqni    = 8  ; advection_config(:)%tracer_names(iqni)    = 'qni'
        iqnr    = 9  ; advection_config(:)%tracer_names(iqnr)    = 'qnr'
        iqns    = 10 ; advection_config(:)%tracer_names(iqns)    = 'qns'
        iqng    = 11 ; advection_config(:)%tracer_names(iqng)    = 'qng'
        iqnh    = 12 ; advection_config(:)%tracer_names(iqnh)    = 'qnh'
        iqnc    = 13 ; advection_config(:)%tracer_names(iqnc)    = 'qnc'
        ininact = 14 ; advection_config(:)%tracer_names(ininact) = 'ninact'

        nqtendphy = 3     !! number of water species for which
                          !! convective and turbulent tendencies are
                          !! stored
        iqm_max   = 7     !! end index of water species mass mixing
                          !! ratios
        iqt       = 15    !! start index of other tracers not related
                          !! at all to moisture

        ntracer = 14

      CASE(7)  ! two-moment scheme with additional prognostic liquid
               ! water (melting) variables for graupel and hail

        iqg     = 6  ; advection_config(:)%tracer_names(iqg)     = 'qg'
        iqh     = 7  ; advection_config(:)%tracer_names(iqh)     = 'qh'
        iqgl    = 8  ; advection_config(:)%tracer_names(iqgl)    = 'qgl'
        iqhl    = 9  ; advection_config(:)%tracer_names(iqhl)    = 'qhl'
        iqni    = 10 ; advection_config(:)%tracer_names(iqni)    = 'qni'
        iqnr    = 11 ; advection_config(:)%tracer_names(iqnr)    = 'qnr'
        iqns    = 12 ; advection_config(:)%tracer_names(iqns)    = 'qns'
        iqng    = 13 ; advection_config(:)%tracer_names(iqng)    = 'qng'
        iqnh    = 14 ; advection_config(:)%tracer_names(iqnh)    = 'qnh'
        iqnc    = 15 ; advection_config(:)%tracer_names(iqnc)    = 'qnc'
        ininact = 16 ; advection_config(:)%tracer_names(ininact) = 'ninact'

        nqtendphy = 3     !! number of water species for which
                          !! convective and turbulent tendencies are
                          !! stored
        iqm_max   = 9     !! end index of water species mass mixing
                          !! ratios
        iqt       = 17    !! start index of other tracers not related
                          !! at all to moisture

        ntracer = 16

      CASE(5)  ! two-moment scheme with CCN and IN budgets

        iqg     = 6  ; advection_config(:)%tracer_names(iqg)     = 'qg'
        iqh     = 7  ; advection_config(:)%tracer_names(iqh)     = 'qh'
        iqni    = 8  ; advection_config(:)%tracer_names(iqni)    = 'qni'
        iqnr    = 9  ; advection_config(:)%tracer_names(iqnr)    = 'qnr'
        iqns    = 10 ; advection_config(:)%tracer_names(iqns)    = 'qns'
        iqng    = 11 ; advection_config(:)%tracer_names(iqng)    = 'qng'
        iqnh    = 12 ; advection_config(:)%tracer_names(iqnh)    = 'qnh'
        iqnc    = 13 ; advection_config(:)%tracer_names(iqnc)    = 'qnc'
        ininact = 14 ; advection_config(:)%tracer_names(ininact) = 'ninact'
        inccn   = 15 ; advection_config(:)%tracer_names(inccn)   = 'nccn'
        ininpot = 16 ; advection_config(:)%tracer_names(ininpot) = 'ninpot'

        nqtendphy = 3     !! number of water species for which
                          !! convective and turbulent tendencies are
                          !! stored
        iqm_max   = 7     !! end index of water species mass mixing
                          !! ratios
        iqt       = 17    !! start index of other tracers not related
                          !! at all to moisture

        ntracer = 16

      CASE(6)

        iqg     = 6  ; advection_config(:)%tracer_names(iqg)     = 'qg'
        iqh     = 7  ; advection_config(:)%tracer_names(iqh)     = 'qh'
        iqni    = 8  ; advection_config(:)%tracer_names(iqni)    = 'qni'
        iqnr    = 9  ; advection_config(:)%tracer_names(iqnr)    = 'qnr'
        iqns    = 10 ; advection_config(:)%tracer_names(iqns)    = 'qns'
        iqng    = 11 ; advection_config(:)%tracer_names(iqng)    = 'qng'
        iqnh    = 12 ; advection_config(:)%tracer_names(iqnh)    = 'qnh'
        iqnc    = 13 ; advection_config(:)%tracer_names(iqnc)    = 'qnc'
        ininact = 14 ; advection_config(:)%tracer_names(ininact) = 'ninact'

        nqtendphy = 3     !! number of water species for which
                          !! convective and turbulent tendencies are
                          !! stored
        iqm_max   = 7     !! end index of water species mass mixing
                          !! ratios
        iqt       = 15    !! start index of other tracers not related
                          !! at all to moisture

        ntracer = 14

      CASE(8)

        iqg     = 6  ; advection_config(:)%tracer_names(iqg)     = 'qg'
        iqh     = 7  ; advection_config(:)%tracer_names(iqh)     = 'qh'


        DO iqb = iqb_i, iqb_e
          iqbin(iqb) = 7+iqb
          DO jg=1,SIZE(advection_config)
            advection_config(jg)%tracer_names(iqbin(iqb)) = qbinname(iqb)
          ENDDO
        END DO
        !$ACC UPDATE DEVICE(iqbin) ASYNC(1)

        iqni    = 7+iqb_e+1 ; advection_config(:)%tracer_names(iqni)    = 'qni'
        iqnr    = 7+iqb_e+2 ; advection_config(:)%tracer_names(iqnr)    = 'qnr'
        iqns    = 7+iqb_e+3 ; advection_config(:)%tracer_names(iqns)    = 'qns'
        iqng    = 7+iqb_e+4 ; advection_config(:)%tracer_names(iqng)    = 'qng'
        iqnh    = 7+iqb_e+5 ; advection_config(:)%tracer_names(iqnh)    = 'qnh'
        iqnc    = 7+iqb_e+6 ; advection_config(:)%tracer_names(iqnc)    = 'qnc'
        ininact = 7+iqb_e+7 ; advection_config(:)%tracer_names(ininact) = 'ninact'

        nqtendphy = 3     !! number of water species for which
                          !! convective and turbulent tendencies are
                          !! stored
        iqm_max = 7+iqb_s    !! stands for qv,qc,qr,qi,qs and defined before
        iqt     = 7+iqb_e+8  !! start index of other tracers not related
                             !! at all to moisture

        ntracer = 7+iqb_e+7  !! total number of tracers. Theis order is the following:
                             !! qv, qc, qi, qr, qs, qg, qh, 
                             !! 33 mass bins for cloud droplets, 
                             !! 33 mass bins for aerosols,
                             !! qni, qnr, qns, qng, qnh, qnc, ninact

      END SELECT ! microphysics schemes


      IF ( (advection_config(1)%iadv_tke) > 0 ) THEN
        IF ( ANY( (/icosmo,iprog/) == inwp_turb(1) ) ) THEN
          iqtke = iqt ; advection_config(:)%tracer_names(iqtke) = 'tke_mc' !! TKE

          ! Note that iqt is not increased, since TKE does not belong
          ! to the hydrometeor group.

          ntracer = ntracer + 1  !! increase total number of tracers by 1

          WRITE(message_text,'(a,i3)') 'Attention: TKE is advected, '//&
            'ntracer is increased by 1 to ',ntracer
          CALL message(routine,message_text)
        ELSE
          WRITE(message_text,'(a,i2)') 'TKE advection not supported for inwp_turb= ', &
            &                          inwp_turb(1)

          CALL finish(routine, message_text )
        ENDIF
      ENDIF

      ! Note: Indices for additional tracers are assigned automatically
      ! via add_tracer_ref in mo_nonhydro_state.

      WRITE(message_text,'(a,i3)') 'Attention: NWP physics is used, '//&
        'ntracer is automatically reset to ',ntracer
      CALL message(routine,message_text)

      ! set the nclass_gscp variable for land-surface scheme to number
      ! of hydrometeor mixing ratios
      DO jg = 1, n_dom
        nclass_gscp(jg) = iqm_max
      ENDDO

    CASE default ! iforcing

      ! set indices for iqv, iqc and iqi dependent on the number of
      ! specified tracers
      !
      iqm_max = 0  ! end index of water species mixing ratios
      !
      iqv = MERGE(1,0,ntracer>=1) ; IF (iqv/=0) iqm_max=1
      iqc = MERGE(2,0,ntracer>=2) ; IF (iqc/=0) iqm_max=2
      iqi = MERGE(3,0,ntracer>=3) ; IF (iqi/=0) iqm_max=3
      !
      iqt = iqm_max+1 ! starting index of non-water species

    END SELECT ! iforcing

    IF (lart) THEN

      IF(iart_ntracer > 0) THEN

        ! assign tracer name strings:
        DO itracer = 1,iart_ntracer
          ! the following code takes care of different character
          ! string lengths:
          src_name_str = ctracer_art(itracer)
          name_len = MIN(LEN_TRIM(src_name_str), VNAME_LEN)
          advection_config(1)%tracer_names(ntracer+itracer) = src_name_str(1:name_len)
        END DO
        ntracer = ntracer + iart_ntracer

      END IF

      WRITE(message_text,'(a,i3,a,i3)') &
        & 'Attention: transport of ART tracers is active, '//&
        'ntracer is increased by ', iart_ntracer, ' to ',ntracer
      CALL message(routine,message_text)
    ENDIF

#ifndef __NO_ICON_COMIN__
    ! loop over the total list of additional requested tracer
    ! variables and add them to the `advection_config` data structure.
    !
    ! Since ICON does not accept a domain-specifc number of tracers, we do this for domain jg=1 only.
    jg=1
    ptr => comin_request_get_list_head()
    VAR_LOOP : DO WHILE (ASSOCIATED(ptr))
      ASSOCIATE (comin_request_item => ptr%item_value)
        ! skip if variable was not requested for this domain
        IF (comin_request_item%descriptor%id /= jg) THEN
          ptr => ptr%next()
          CYCLE VAR_LOOP
        END IF

        IF (.NOT. comin_request_item%metadata%tracer) THEN
          ptr => ptr%next()
          CYCLE VAR_LOOP
        END IF

        ntracer = ntracer + 1
        advection_config(jg)%tracer_names(ntracer) = comin_request_item%descriptor%name

        WRITE(message_text,*) 'Attention: ComIn active, adding tracer ', &
          &                   TRIM(comin_request_item%descriptor%name)
        CALL message(routine,message_text)

      END ASSOCIATE
      ptr => ptr%next()
    END DO VAR_LOOP

    WRITE(message_text,'(a,i3)') 'Attention: ComIn active, '//&
      &   'ntracer is increased to ',ntracer
    CALL message(routine,message_text)

    ! Seperate loop for adding the tracer to turbulence/convection
    ! Generally allowing a domain-dependent specification of lturb and/or lconv
    DOM_LOOP : DO jg = 1, n_dom
      ptr => comin_request_get_list_head()
      VAR_LOOP_TURB : DO WHILE (ASSOCIATED(ptr))
        ASSOCIATE (comin_request_item => ptr%item_value)
          ! skip if variable was not requested for this domain
          IF (comin_request_item%descriptor%id /= jg) THEN
            ptr => ptr%next()
            CYCLE VAR_LOOP_TURB
          END IF

          IF (.NOT. comin_request_item%metadata%tracer) THEN
            ptr => ptr%next()
            CYCLE VAR_LOOP_TURB
          END IF

          IF (comin_request_item%metadata%tracer_turb) THEN
            comin_config%comin_icon_domain_config(jg)%nturb_tracer = &
              &  comin_config%comin_icon_domain_config(jg)%nturb_tracer + 1
          END IF

          IF (comin_request_item%metadata%tracer_conv) THEN
            comin_config%comin_icon_domain_config(jg)%nconv_tracer = &
              &  comin_config%comin_icon_domain_config(jg)%nconv_tracer + 1
          END IF

        END ASSOCIATE
        ptr => ptr%next()
      END DO VAR_LOOP_TURB
      WRITE (message_text,'(A,I2,A,I2,A)') "Domain ",jg," contains ",                              &
        &                                  comin_config%comin_icon_domain_config(jg)%nturb_tracer, &
        &                                  " comin tracer registered for turbulent transport"
      CALL message(routine, message_text)
      WRITE (message_text,'(A,I2,A,I2,A)') "Domain ",jg," contains ",                              &
        &                                  comin_config%comin_icon_domain_config(jg)%nconv_tracer, &
        &                                  " comin tracer registered for convective transport"
      CALL message(routine, message_text)

      IF ( (inwp_turb(jg) /= icosmo) .AND. &
        & comin_config%comin_icon_domain_config(jg)%nturb_tracer > 0) THEN
        WRITE (message_text,'(A,A,I2)') "ComIn tracer-turbulence interaction only valid ", &
          &                           "for inwp_turb = ",icosmo
        CALL finish(routine, message_text)
      END IF
      IF ( (inwp_convection(jg) /= 1) .AND. &
        & comin_config%comin_icon_domain_config(jg)%nconv_tracer > 0) THEN
        WRITE (message_text,'(A,A,I2)') "ComIn tracer-convection interaction only valid ", &
          &                           "for inwp_convection = 1"
        CALL finish(routine, message_text)
      END IF
    END DO DOM_LOOP
#endif

    ! take into account additional passive tracers, if present
    ! no need to update iqt, since passive tracers do not belong to
    ! the hydrometeor group.
    ! ATTENTION: as long as ntracer is not domain specific, we set jg=1
    IF ( advection_config(1)%npassive_tracer > 0) THEN
      ntracer = ntracer + advection_config(1)%npassive_tracer
      WRITE(message_text,'(a,i3,a,i3)') 'Attention: passive tracers have been added, '//&
        'ntracer is increased by ',advection_config(1)%npassive_tracer, &
        ' to ',ntracer
      CALL message(routine,message_text)
    ENDIF


    IF (ltransport) THEN

      DO jg = 1,n_dom

        SELECT CASE ( iforcing )
        CASE ( INWP )
          !...........................................................
          ! in NWP physics
          !...........................................................

          ! Force settings for tracer iqtke, if TKE advection is
          ! performed
          !
          IF ( advection_config(jg)%iadv_tke > 0 ) THEN

            ! force monotonous slope limiter for vertical advection
            advection_config(jg)%itype_vlimit(iqtke) = 2

            ! force positive definite flux limiter for horizontal advection
            advection_config(jg)%itype_hlimit(iqtke) = 4

            SELECT CASE (advection_config(jg)%iadv_tke)
            CASE (1)
              ! switch off horizontal advection
              advection_config(jg)%ihadv_tracer(iqtke) = 0
            CASE (2)
              ! check whether horizontal substepping is switched on
              IF (ALL( (/22,32,42,52/) /= advection_config(jg)%ihadv_tracer(iqtke)) ) THEN
                ! choose Miura with substepping
                advection_config(jg)%ihadv_tracer(iqtke) = 22
              ENDIF
            END SELECT

          ENDIF

        END SELECT ! iforcing

      END DO ! jg = 1,n_dom
    END IF ! ltransport

  END SUBROUTINE init_tracer_settings


END MODULE mo_advection_utils
