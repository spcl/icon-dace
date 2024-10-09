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

! Contains the setup of variables related to horizontal diffusion

MODULE mo_diffusion_nml

  USE mo_diffusion_config,    ONLY: diffusion_config
  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_exception,           ONLY: message, finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_diffusion_namelist

  !-------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary parameters setting up the
  !     configuration of the dynamical core
  !-------------------------------------------------------------------------
  INTEGER :: hdiff_order  ! order of horizontal diffusion
                          ! -1: no diffusion
                          ! 2: 2nd order linear diffusion on all vertical levels
                          ! 4: 4th order linear diffusion on all vertical levels
                          ! 5: Smagorinsky diffusion with optional fourth-order background diffusion

  REAL(wp) :: hdiff_efdt_ratio      ! ratio of e-folding time to (2*)time step
  REAL(wp) :: hdiff_w_efdt_ratio    ! ratio of e-folding time to time step for w diffusion (NH only)
  REAL(wp) :: hdiff_min_efdt_ratio  ! minimum value of hdiff_efdt_ratio 
                                    ! (for upper sponge layer)
  REAL(wp) :: hdiff_tv_ratio        ! the ratio of diffusion coefficient: temp:mom
  REAL(wp) :: hdiff_smag_fac        ! scaling factor for Smagorinsky diffusion at height hdiff_smag_z and below
  REAL(wp) :: hdiff_smag_fac2       ! scaling factor for Smagorinsky diffusion at height hdiff_smag_z2
  REAL(wp) :: hdiff_smag_fac3       ! scaling factor for Smagorinsky diffusion at height hdiff_smag_z3
  REAL(wp) :: hdiff_smag_fac4       ! scaling factor for Smagorinsky diffusion at height hdiff_smag_z4 and above
  REAL(wp) :: hdiff_smag_z          ! height up to which hdiff_smag_fac is used, start of linear profile
  REAL(wp) :: hdiff_smag_z2         ! height of hdiff_smag_fac2, end of linear and start of quadratic profile
  REAL(wp) :: hdiff_smag_z3         ! height of hdiff_smag_fac3, to define quadratic profile
  REAL(wp) :: hdiff_smag_z4         ! height from which hdiff_smag_fac4, end of quadratic profile
  REAL(wp) :: hdiff_multfac         ! multiplication factor of normalized diffusion
                                    ! coefficient for nested domains
  INTEGER  :: itype_vn_diffu        ! options for discretizing the Smagorinsky momentum diffusion
  INTEGER  :: itype_t_diffu         ! options for discretizing the Smagorinsky temperature diffusion
  LOGICAL :: lhdiff_temp   ! if .TRUE., apply horizontal diffusion to temp.
  LOGICAL :: lhdiff_vn     ! if .TRUE., apply horizontal diffusion to horizontal momentum.
  LOGICAL :: lhdiff_w      ! if .TRUE., apply horizontal diffusion to vertical momentum.
  LOGICAL :: lhdiff_q      ! if .TRUE., apply horizontal diffusion to QV and QC.
  LOGICAL :: lsmag_3d(max_dom)      ! if .TRUE., compute 3D Smagorinsky diffusion coefficient.
  LOGICAL :: lhdiff_smag_w(max_dom) ! if .TRUE., apply additional Smagorinsky diffusion to vertical momentum.

  NAMELIST/diffusion_nml/ hdiff_order,                                                       &
                          hdiff_efdt_ratio, hdiff_min_efdt_ratio, hdiff_tv_ratio,            &
                          hdiff_smag_fac, hdiff_smag_fac2, hdiff_smag_fac3, hdiff_smag_fac4, &
                          hdiff_smag_z,   hdiff_smag_z2,   hdiff_smag_z3,   hdiff_smag_z4,   &
                          hdiff_multfac, lhdiff_temp, lhdiff_vn, itype_vn_diffu,             &
                          itype_t_diffu, hdiff_w_efdt_ratio, lhdiff_w, lsmag_3d, lhdiff_smag_w, lhdiff_q

CONTAINS
  !-------------------------------------------------------------------------
  !! Read Namelist for horizontal diffusion. 
  !!
  !! This subroutine 
  !! - reads the Namelist for diffusion
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)  
  !!
  SUBROUTINE read_diffusion_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename 
    INTEGER :: istat, funit, iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_diffusion_nml: read_diffusion_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------
    lhdiff_temp          = .TRUE.
    lhdiff_vn            = .TRUE.
    lhdiff_w             = .TRUE.
    lhdiff_q             = .FALSE.
    lsmag_3d(:)          = .FALSE.
    lhdiff_smag_w(:)     = .FALSE.

    hdiff_order          = 5
    hdiff_efdt_ratio     = 36.0_wp
    hdiff_smag_fac       = 0.015_wp
    hdiff_smag_fac2      = REAL(2.e-6*(1600. + 25000. + SQRT(1600.*(1600.+50000.))),wp) ! (1)
    hdiff_smag_fac3      = 0.000_wp
    hdiff_smag_fac4      = 1.000_wp
    hdiff_smag_z         = 32500._wp
    hdiff_smag_z2        = REAL(       1600. + 50000. + SQRT(1600.*(1600.+50000.)) ,wp) ! (1)
    hdiff_smag_z3        = 50000._wp
    hdiff_smag_z4        = 90000._wp
    !
    ! (1) These irrational values are computed in single precision and then converted
    !     to working precision, meaning double precision, to prevent that these values
    !     differ after writing to and reading from the restart file. If this happens,
    !     then a restart changes results compared to a reference simulation without restart.
    !     This problem was found in restart tests with executables built with Intel compilers
    !     when these values were computed with all numbers in working precision.
    !     Using single precision computations solved the problem.
    !

    hdiff_min_efdt_ratio = 1.0_wp
    hdiff_w_efdt_ratio   = 15.0_wp
    hdiff_multfac        = 1.0_wp
    hdiff_tv_ratio       = 1.0_wp
    itype_vn_diffu       = 1
    itype_t_diffu        = 2


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('diffusion_nml')
      READ(funit,NML=diffusion_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('diffusion_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, diffusion_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, diffusion_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, diffusion_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------
    SELECT CASE( hdiff_order)
    CASE(-1)
      CALL message(TRIM(routine),'Horizontal diffusion switched off.')
      lhdiff_temp = .FALSE.
      lhdiff_vn   = .FALSE.
      lhdiff_w    = .FALSE.

    CASE(2,4,5)

      IF ((.NOT.lhdiff_temp).AND.(.NOT.lhdiff_vn)) THEN
        CALL message('','')
        CALL message('','lhdiff_temp and lhdiff_vn both set to .FALSE. by user.')
        CALL message('','Horizontal diffusion is thus switched off and '// &
                        'hdiff_order reset to -1')
        CALL message('','')
        hdiff_order = -1
      END IF

    CASE DEFAULT
      CALL finish(TRIM(routine),                     &
        & 'Error: Invalid choice of hdiff_order. '// &
        & 'Choose from -1, 2, 4, and 5.')
    END SELECT

    IF ( hdiff_efdt_ratio<=0._wp) THEN
      CALL message(TRIM(routine),'No horizontal background diffusion is used')
    ENDIF

    ! Checks for hdiff_smag parameters
    IF ( hdiff_smag_fac  < 0.0_wp ) THEN
       CALL finish( TRIM(routine), 'hdiff_smag_fac  < 0 not allowed')
    ENDIF
    IF ( hdiff_smag_fac2 < 0.0_wp ) THEN
       CALL finish( TRIM(routine), 'hdiff_smag_fac2 < 0 not allowed')
    ENDIF
    IF ( hdiff_smag_fac3 < 0.0_wp ) THEN
       CALL finish( TRIM(routine), 'hdiff_smag_fac3 < 0 not allowed')
    ENDIF
    IF ( hdiff_smag_fac4 < 0.0_wp ) THEN
       CALL finish( TRIM(routine), 'hdiff_smag_fac4 < 0 not allowed')
    ENDIF
    IF ( hdiff_smag_z2 <= hdiff_smag_z  ) THEN
       CALL finish( TRIM(routine), 'hdiff_smag_z2 <= hdiff_smag_z  not allowed')
    ENDIF
    IF ( hdiff_smag_z2 >= hdiff_smag_z4 ) THEN
       CALL finish( TRIM(routine), 'hdiff_smag_z2 >= hdiff_smag_z4 not allowed')
    ENDIF
    IF ( hdiff_smag_z3 == hdiff_smag_z2 ) THEN
       CALL finish( TRIM(routine), 'hdiff_smag_z3 == hdiff_smag_z2 not allowed')
    ENDIF
    IF ( hdiff_smag_z3 == hdiff_smag_z4 ) THEN
       CALL finish( TRIM(routine), 'hdiff_smag_z3 == hdiff_smag_z4 not allowed')
    ENDIF

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------
    diffusion_config(:)% lhdiff_temp          =  lhdiff_temp
    diffusion_config(:)% lhdiff_vn            =  lhdiff_vn
    diffusion_config(:)% lhdiff_w             =  lhdiff_w
    diffusion_config(:)% lhdiff_q             =  lhdiff_q
    diffusion_config(:)% lsmag_3d             =  lsmag_3d(:)
    diffusion_config(:)% lhdiff_smag_w        =  lhdiff_smag_w(:)
    diffusion_config(:)% hdiff_order          =  hdiff_order
    diffusion_config(:)% hdiff_efdt_ratio     =  hdiff_efdt_ratio
    diffusion_config(:)% hdiff_w_efdt_ratio   =  hdiff_w_efdt_ratio
    diffusion_config(:)% hdiff_min_efdt_ratio =  hdiff_min_efdt_ratio
    diffusion_config(:)% hdiff_smag_fac       =  hdiff_smag_fac
    diffusion_config(:)% hdiff_smag_fac2      =  hdiff_smag_fac2
    diffusion_config(:)% hdiff_smag_fac3      =  hdiff_smag_fac3
    diffusion_config(:)% hdiff_smag_fac4      =  hdiff_smag_fac4
    diffusion_config(:)% hdiff_smag_z         =  hdiff_smag_z
    diffusion_config(:)% hdiff_smag_z2        =  hdiff_smag_z2
    diffusion_config(:)% hdiff_smag_z3        =  hdiff_smag_z3
    diffusion_config(:)% hdiff_smag_z4        =  hdiff_smag_z4
    diffusion_config(:)% hdiff_multfac        =  hdiff_multfac
    diffusion_config(:)% hdiff_tv_ratio       =  hdiff_tv_ratio 
    diffusion_config(:)%itype_vn_diffu        =  itype_vn_diffu
    diffusion_config(:)%itype_t_diffu         =  itype_t_diffu 

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=diffusion_nml)                    
      CALL store_and_close_namelist(funit,'diffusion_nml') 
    ENDIF
    ! 7. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=diffusion_nml)

  END SUBROUTINE read_diffusion_namelist

END MODULE mo_diffusion_nml
