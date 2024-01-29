! Namelist for configuration of grid refinement algorithms.
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

MODULE mo_gridref_nml

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile

  USE mo_gridref_config,      ONLY:   &
    &                            config_rbf_vec_kern_grf_e => rbf_vec_kern_grf_e,& 
    &                            config_rbf_scale_grf_e    => rbf_scale_grf_e,&
    &                            config_grf_velfbk         => grf_velfbk,&
    &                            config_grf_scalfbk        => grf_scalfbk,&
    &                            config_grf_tracfbk        => grf_tracfbk,&
    &                            config_grf_intmethod_c    => grf_intmethod_c,&
    &                            config_grf_intmethod_e    => grf_intmethod_e,&
    &                            config_grf_intmethod_ct   => grf_intmethod_ct,&
    &                            config_l_density_nudging  => l_density_nudging,&
    &                            config_denom_diffu_v      => denom_diffu_v,&
    &                            config_denom_diffu_t      => denom_diffu_t,&
    &                            config_fbk_relax_timescale => fbk_relax_timescale
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_gridref_namelist

  !---------------------
  ! namelist variables
  !---------------------

  INTEGER  :: rbf_vec_kern_grf_e ! rbf kernel for vector interpolation

  ! scale factors for rbf grid refinement interpolation
  REAL(wp) :: rbf_scale_grf_e(max_dom)

  INTEGER  :: grf_intmethod_c,  &  ! switch for type of grid refinement interpolation
    &         grf_intmethod_ct, &  ! (see below for explanation of options)
    &         grf_intmethod_e

  INTEGER  :: grf_velfbk     ! switch for velocity feedback method
                             ! 1 = averaging over child edges 1 and 2;
                             ! 2 = 2nd-order method using RBF reconstruction to child vertices
  
  INTEGER  :: grf_scalfbk    ! switch for feedback method of scalar dynamical variables
                             ! 1 = area-weighted averaging
                             ! 2 = bilinear interpolation

  INTEGER  :: grf_tracfbk    ! switch for feedback method of passive tracer variables
                             ! 1 = area-weighted averaging
                             ! 2 = bilinear interpolation

  LOGICAL  :: l_density_nudging ! .true.: apply density nudging near lateral nest boundaries if feedback is turned on
                                ! (in case of one-way nesting, all prognostic variables are nudged irrespective of this switch)

  ! Denominators of normalized diffusion coefficients for boundary diffusion
  REAL(wp) :: denom_diffu_v, denom_diffu_t

  ! Relaxation time scale for feedback in case of ifeedback_type = 2
  REAL(wp) :: fbk_relax_timescale

  NAMELIST/gridref_nml/  rbf_vec_kern_grf_e, rbf_scale_grf_e,                   &
    &                    grf_velfbk, grf_scalfbk, grf_tracfbk,                  &
    &                    grf_intmethod_c, grf_intmethod_e,                      &
    &                    grf_intmethod_ct, denom_diffu_v, denom_diffu_t,        &
    &                    l_density_nudging, fbk_relax_timescale

CONTAINS
  !-------------------------------------------------------------------------
  !! This subroutine 
  !! - reads the Namelist for local grid refinement 
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  SUBROUTINE read_gridref_namelist( filename )

    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER :: istat, funit, iunit
    CHARACTER(len=*),PARAMETER :: routine = 'mo_gridref_nml:read_gridref_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    ! Switch for interpolation method used for cell-based dynamical
    grf_intmethod_c   = 2     ! 1: copying, 2: gradient-based interpolation
    ! Switch for interpolation method used for tracer variables
    grf_intmethod_ct  = 2     ! 1: copying, 2: gradient-based interpolation
    ! Currently, grf_intmethod_c is used for temperature only; other variables are copied
    grf_intmethod_e   = 6     ! 1: (removed),
                              ! 2: RBF,
                              ! 3: (removed),
                              ! 4: RBF/gradient-based,
                              ! 5: (removed),
                              ! 6: same as 4, but direct interpolation of mass fluxes

    ! Switch for velocity feedback method.
    grf_velfbk      = 1       ! 1: average over child edges 1 and 2
                              ! 2: 2nd-order method using RBF reconstruction to child vertices
    ! Switch for feedback method for scalar dynamical variables
    grf_scalfbk     = 2       ! 1: area-weighted averaging
                              ! 2: bilinear interpolation
    ! Switch for feedback method for passive tracer variables
    grf_tracfbk     = 2       ! 1: area-weighted averaging
                              ! 2: bilinear interpolation

    ! RBF kernels for grid refinement interpolation
    rbf_vec_kern_grf_e = 1    ! 1: Gaussian, 2: 1/(1+r**2), 3: inverse multiquadric

    ! Initialize namelist fields for scaling factors (dimension 1:max_dom)
    rbf_scale_grf_e(1:max_dom) = -1.0_wp  ! dummy value; resolution-dependent defaults are set in gridref_config

    ! Denominator for temperature boundary diffusion
    denom_diffu_t = 135._wp

    ! Denominator for velocity boundary diffusion
    denom_diffu_v = 200._wp

    ! Density nudging near nest boundaries turned off by default
    ! only applicable for grf_intmethod_e == 2 or 4
    l_density_nudging = .FALSE.

    ! Relaxation time scale for feedback in case of ifeedback_type = 2
    fbk_relax_timescale = 10800._wp ! 3 hours

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('gridref_nml')
      READ(funit,NML=gridref_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('gridref_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, gridref_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, gridref_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, gridref_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml


    !--------------------------------------------------------------------
    ! Sanity and internal checks
    !--------------------------------------------------------------------

    IF (ALL((/2,4,6/) /= grf_intmethod_e)) THEN
      WRITE(message_text,'(a,i2,a)') 'Invalid value grf_intmethod_e=',grf_intmethod_e, &
        &                       ' (must be 2, 4, or 6)'
      CALL finish(TRIM(routine), message_text)
    ENDIF


    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
      config_rbf_vec_kern_grf_e = rbf_vec_kern_grf_e
      config_rbf_scale_grf_e = rbf_scale_grf_e
      config_grf_velfbk = grf_velfbk
      config_grf_scalfbk = grf_scalfbk
      config_grf_tracfbk = grf_tracfbk
      config_grf_intmethod_c = grf_intmethod_c
      config_grf_intmethod_e = grf_intmethod_e
      config_grf_intmethod_ct = grf_intmethod_ct
      config_denom_diffu_v = denom_diffu_v
      config_denom_diffu_t = denom_diffu_t
      config_l_density_nudging = l_density_nudging
      config_fbk_relax_timescale = fbk_relax_timescale

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=gridref_nml)
      CALL store_and_close_namelist(funit, 'gridref_nml')
    ENDIF
    ! 6. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=gridref_nml)

  END SUBROUTINE read_gridref_namelist

END MODULE mo_gridref_nml
