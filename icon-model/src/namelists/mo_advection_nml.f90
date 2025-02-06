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

! Namelist for scalar transport
!
! These subroutines are called by  read_atmo_namelists and do the transport
! setup.

MODULE mo_advection_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message_text, message, print_value
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_run_config,          ONLY: ntracer
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_ntracer, max_dom,      &
    &                               MIURA, FFSL_HYB_MCYCL, ippm_v,              &
    &                               inol, ifluxl_m, ifluxl_sm, inol_v,          &
    &                               islopel_vsm, ifluxl_vpd, VNAME_LEN, NO_VADV
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_advection_config,    ONLY: advection_config 
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_transport_namelist




  !----------------------------------!
  ! transport_nml namelist variables !
  !----------------------------------!

  CHARACTER(len=VNAME_LEN) ::  &   !< tracer-specific name suffixes  
    &  tracer_names(MAX_NTRACER)   !< these are only required for 
                                   !< idealized runs without NWP or ECHAM forcing.

  CHARACTER(len=MAX_CHAR_LENGTH) :: &!< list of tracers to initialize
    &  ctracer_list


  INTEGER :: &                     !< selects horizontal transport scheme
    &  ihadv_tracer(max_ntracer)   !< 0:  no horizontal advection
                                   !< 2:  2nd order miura
                                   !< 3:  3rd order miura with quadr./cubic reconstr.
                                   !< 4:  Flux form semi lagrange (FFSL)
                                   !< 5:  hybrid FFSL Miura3
                                   !< 20: subcycling version of miura
                                   !< 22: 2nd order miura and miura_cycl
                                   !< 32: 3rd order miura with miura_cycl
                                   !< 42: FFSL with miura_cycl
                                   !< 52: FFSL_HYB with miura_cycl

  INTEGER :: &                     !< selects vertical transport scheme
    &  ivadv_tracer(max_ntracer)   !< 0 : no vertical advection
                                   !< 1 : 1st order upwind
                                   !< 2 : 3rd order PSM for CFL>                            
                                   !< 3 : 3rd order PPM for CFL>               

  INTEGER :: &                     !< advection of TKE
    &  iadv_tke                    !< 0 : none
                                   !< 1 : vertical advection only
                                   !< 2 : vertical and horizontal advection


  LOGICAL :: lvadv_tracer          !< if .TRUE., calculate vertical tracer advection
  LOGICAL :: lclip_tracer          !< if .TRUE., clip negative tracer values

  LOGICAL :: llsq_svd              !< least squares reconstruction with 
                                   !< singular value decomposition (TRUE) or 
                                   !< QR decomposition (FALSE) of design matrix A

  INTEGER :: &                     !< parameter used to select the limiter
    &  itype_vlimit(max_ntracer)   !< for vertical transport
                               

  INTEGER :: &                     !< parameter used to select the limiter
    &  itype_hlimit(max_ntracer)   !< for horizontal transport
                                   !< 0: no limiter
                                   !< 3: monotonous flux limiter
                                   !< 4: positive definite flux limiter

  INTEGER :: &                     !< additional method for identifying and avoiding 
    & ivlimit_selective(max_ntracer)!< spurious limiting of smooth extrema
                                   !< 1: switch on
                                   !< 0: switch off

  REAL(wp):: beta_fct              !< global boost factor for range of permissible values in 
                                   !< (semi-) monotonous flux limiter. A value larger than 
                                   !< 1 allows for (small) over and undershoots, while a value 
                                   !< of 1 gives strict monotonicity (at the price of increased 
                                   !< diffusivity).

  INTEGER :: igrad_c_miura         !< parameter used to select the gradient
                                   !< reconstruction method at cell center
                                   !< for second order miura scheme

  INTEGER :: ivcfl_max             !< determines stability range of vertical 
                                   !< ppm-scheme (approximate allowable maximum 
                                   !< CFL-number)

  INTEGER :: npassive_tracer       !< number of additional passive tracers, in addition to
                                   !< microphysical- and ART tracers. 

  INTEGER :: nadv_substeps(max_dom)!< domain specific number of substeps per fast physics time step 
                                   !< for the Miura-type substepping schemes 20, 22, 32, 42, 52

  CHARACTER(len=MAX_CHAR_LENGTH) :: &!< Comma separated list of initialization formulae 
    &  init_formula                  !< for passive tracers.

  NAMELIST/transport_nml/ ihadv_tracer, ivadv_tracer, lvadv_tracer,       &
    &                     itype_vlimit, ivlimit_selective,                &
    &                     ivcfl_max, itype_hlimit,                        &
    &                     iadv_tke, beta_fct,                             &
    &                     lclip_tracer, tracer_names,                     &
    &                     ctracer_list, igrad_c_miura, llsq_svd,          &
    &                     npassive_tracer, init_formula, nadv_substeps



CONTAINS


  !-------------------------------------------------------------------------
  !
  !! Read Namelist for tracer transport. 
  !!
  !! This subroutine 
  !! - reads the Namelist for tracer transport
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  SUBROUTINE read_transport_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: iunit
    INTEGER :: it          !< loop counter
    INTEGER :: nname       !< number of names given in transport_nml/tracer_names
    CHARACTER(len=VNAME_LEN) :: tname

    INTEGER, PARAMETER :: default_nadv_substeps = 3 ! default value for nadv_substeps

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_advection_nml: read_transport_nml'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    ctracer_list         = ''
    ihadv_tracer(:)      = MIURA     ! miura horizontal advection scheme
    itype_hlimit(:)      = ifluxl_sm ! positive definite flux limiter
    ivadv_tracer(:)      = ippm_v    ! PPM vertical advection scheme
    itype_vlimit(:)      = islopel_vsm ! semi-monotonous slope limiter
    ivlimit_selective(:) = 0         ! do not use method for preserving smooth extrema
    iadv_tke             = 0         ! no TKE advection
    beta_fct             = 1.005_wp  ! factor of allowed over-/undershooting in monotonous limiter
    ivcfl_max            = 5         ! CFL-stability range for vertical advection
    lvadv_tracer         = .TRUE.    ! vertical advection yes/no
    lclip_tracer         = .FALSE.   ! clipping of negative values yes/no

    igrad_c_miura        = 1         ! MIURA linear least squares reconstruction

    llsq_svd             = .TRUE.    ! apply singular-value-decomposition (FALSE: use QR-decomposition)
    npassive_tracer      = 0         ! no additional passive tracers
    init_formula         = ''        ! no explizit initialization of passive tracers

    tracer_names(:) = '...'          ! default name for undefined tracers

    nadv_substeps(:)     = default_nadv_substeps  ! domain specific number of substeps per fast physics time step 
                                                  ! for the Miura-type substepping schemes 20, 22, 32, 42, 52

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('transport_nml')
      READ(funit,NML=transport_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, transport_nml)   ! write defaults to temporary text file
    END IF

    CALL open_nml(TRIM(filename))
    CALL position_nml ('transport_nml', STATUS=istat)
    !
    SELECT CASE (istat)
    CASE (POSITIONED)

      ! Set array parameters to dummy values to determine which ones are actively set in the namelist
      !
      nadv_substeps(:) = -1

      READ (nnml, transport_nml)                                        ! overwrite default settings

      ! Restore default values for global domain where nothing at all has been specified
      ! 
      IF (nadv_substeps(1) < 0) nadv_substeps(1) = default_nadv_substeps

      ! Copy values of parent domain (in case of linear nesting) to nested domains where nothing has been specified
      !
      DO jg = 2, max_dom
        IF (nadv_substeps(jg) < 0) nadv_substeps(jg) = nadv_substeps(jg-1)
      ENDDO

      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, transport_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------

    ! check consistency of run_nml/ntracer and the number of tracer
    ! names given in transport
    !
    ! Count the number of names given in tracer_names(:) before the first
    ! default name '...' is found.
    CALL message(' ',' ')
    CALL message(routine,'Evaluate the number of names in transport_nml/tracer_names')
    nname=0
    DO
      IF (TRIM(ADJUSTL(tracer_names(nname+1))) /= '...') THEN
        nname=nname+1
        CALL print_value('tracer name '//TRIM(ADJUSTL(tracer_names(nname)))//' in position',nname)
      ELSE
        CALL print_value('number of named tracers',nname)
        EXIT
      END IF
    END DO
    !
    ! Set ntracer if zero so far but tracer names are given
    IF (ntracer==0 .AND. nname > 0) THEN
      ntracer=nname
      CALL print_value('ntracer changed from 0 to',nname)
    END IF
    !
    ! Stop if ntracer has been set to non-zero but less than the number of tracer names
    IF (0<ntracer .AND. ntracer<nname) THEN
      CALL finish(routine,'Found inconsistent setup: 0 < run_nml/ntracer < ' //&
           &      'number of names in transport_nml/tracer_names')
    END IF
    !
    ! If ntracer has been set and ntracer > nname, then fill in default names 'q<no>'
    ! for tracers with indices nname+1 to ntracer.
    IF (nname<ntracer) THEN
      CALL message(routine,'Setting default names for tracers without names.')
      CALL print_value('number of tracers without names',ntracer-nname)
      DO it=nname+1,ntracer
        WRITE(tname,'(i3)') it
        tracer_names(it) = 'q'//TRIM(ADJUSTL(tname))
        CALL print_value('tracer name '//TRIM(ADJUSTL(tracer_names(it)))//' in position',it)
      END DO
    END IF
    CALL message(' ',' ')

    ! flux computation methods - sanity check
    !
    IF ( ANY(ihadv_tracer(1:max_ntracer) > FFSL_HYB_MCYCL) .OR.       &
      &  ANY(ihadv_tracer(1:max_ntracer) < 0) )    THEN
      CALL finish(routine,                                              &
        &  'incorrect settings for ihadv_tracer. Must be 0,1,2,3,4,5,'//&
        &  '20,22,32,42 or 52 ')
    ENDIF

    IF ( ANY(ivadv_tracer(1:max_ntracer) > ippm_v) .OR.               &
      &  ANY(ivadv_tracer(1:max_ntracer) < 0) ) THEN
      CALL finish(routine,                                     &
        &  'incorrect settings for ivadv_tracer. Must be 0,1,2, or 3 ')
    ENDIF

    IF (.NOT. lvadv_tracer) THEN
      IF ( ANY(ivadv_tracer(1:max_ntracer) /= NO_VADV) ) THEN
        ivadv_tracer(1:max_ntracer) = NO_VADV
        WRITE(message_text,'(a,I2,a,a)') &
          &  'Resetting ivadv_tracer(1:max_ntracer)=',NO_VADV,' for all tracers, ', &
          &  'as vertical tracer transport is disabled (lvadv_tracer=.FALSE.).'
        CALL message(TRIM(routine), message_text)
      ENDIF
    ENDIF



    ! limiter - sanity check
    !
    IF ( ANY(itype_vlimit(1:max_ntracer) < inol_v ) .OR.              &
      &  ANY(itype_vlimit(1:max_ntracer) > ifluxl_vpd)) THEN
      CALL finish(routine,                                     &
        &  'incorrect settings for itype_vlimit. Permissible choices [0,..,3]')
    ENDIF
    DO it = 1, max_ntracer
      IF ( ALL((/inol,ifluxl_m,ifluxl_sm/) /= itype_hlimit(it)) ) THEN
        CALL finish(routine,                                     &
          &  'incorrect settings for itype_hlimit. Permissible choices [0,3,4]')
      ENDIF
    ENDDO
    IF ( ANY(ivlimit_selective(1:max_ntracer) < 0 ) .OR.              &
      &  ANY(ivlimit_selective(1:max_ntracer) > 1 )) THEN
      CALL finish(routine,                                     &
        &  'incorrect settings for ivlimit_selective. Permissible values 0 or 1')
    ENDIF

    ! FCT multiplicative spreading - sanity check
    ! 
    IF ( (beta_fct < 1.0_wp) .OR. (beta_fct >= 2.0_wp) ) THEN
      CALL finish(routine,                                     & 
        & 'permissible range of values for beta_fct:  [1,2[ ')
    ENDIF


    ! TKE advection - sanity check
    !
    IF ( ALL((/0,1,2/) /= iadv_tke) ) THEN
      CALL finish(routine,                                     &
        &  'incorrect settings for iadv_tke. Must be 0, 1 or 2')
    ENDIF

    IF (ctracer_list /= '') THEN
      WRITE(message_text,'(a,a)') &
        &  'ctracer_list is obsolescent and inactive. ', &
        &  'Please use variables tracer_names and/or tracer_inidist_list, instead'
      CALL message('WARNING (transport_nml):', message_text)
    ENDIF

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 0,max_dom
      advection_config(jg)%tracer_names(:)     = ADJUSTL(tracer_names(:))
      advection_config(jg)%nname               = nname
      advection_config(jg)%ihadv_tracer(:)     = ihadv_tracer(:)
      advection_config(jg)%ivadv_tracer(:)     = ivadv_tracer(:)
      advection_config(jg)%lvadv_tracer        = lvadv_tracer
      advection_config(jg)%lclip_tracer        = lclip_tracer
      advection_config(jg)%llsq_svd            = llsq_svd
      advection_config(jg)%itype_vlimit(:)     = itype_vlimit(:)
      advection_config(jg)%itype_hlimit(:)     = itype_hlimit(:)
      advection_config(jg)%ivlimit_selective(:)= ivlimit_selective(:)
      advection_config(jg)%beta_fct            = beta_fct
      advection_config(jg)%igrad_c_miura       = igrad_c_miura
      advection_config(jg)%iadv_tke            = iadv_tke
      advection_config(jg)%ivcfl_max           = ivcfl_max
      advection_config(jg)%npassive_tracer     = npassive_tracer 
      IF (jg==0) THEN
        advection_config(jg)%nadv_substeps     = -1   ! dummy value, as nadv_substeps is only defined from dom01 onwards.
      ELSE
        advection_config(jg)%nadv_substeps     = nadv_substeps(jg) 
      ENDIF
      advection_config(jg)%init_formula        = init_formula 
    ENDDO



    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=transport_nml)                    
      CALL store_and_close_namelist(funit, 'transport_nml')             
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=transport_nml)


  END SUBROUTINE read_transport_namelist

END MODULE mo_advection_nml
