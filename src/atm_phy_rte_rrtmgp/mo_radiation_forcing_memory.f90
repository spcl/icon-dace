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

! Data types and variables used by the rte-rrtmgp radiation package.
!
! This module contains
!  definition of data types for organising rte-rrtmgp diagnostics,
!  the actual variables that are declared of these types, and
!  subroutines for (de-)allocating memory for the variables.

MODULE mo_radiation_forcing_memory

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_exception,           ONLY: message, finish
  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_parallel_config,     ONLY: nproma
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_model_domain,        ONLY: t_patch
  USE mo_var_list,            ONLY: t_var_list_ptr, add_var
  USE mo_var_list_register,   ONLY: vlr_add, vlr_del
  USE mo_var_metadata,        ONLY: create_vert_interp_metadata, vintp_types
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_PACK16, &
    &                               DATATYPE_FLT32,  DATATYPE_FLT64,   &
    &                               GRID_UNSTRUCTURED
  USE mo_zaxis_type,          ONLY: ZA_REFERENCE, ZA_REFERENCE_HALF,         &
    &                               ZA_SURFACE
  USE mo_radiation_config,    ONLY: lradforcing
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: prm_radiation_forcing                           !< variables
  PUBLIC :: prm_radiation_forcing_list                      !< variable lists
  PUBLIC :: construct_radiation_forcing_list                !< subroutine
  PUBLIC :: destruct_radiation_forcing_list                 !< subroutines
  PUBLIC :: t_radiation_forcing                             !< derived types

  !!--------------------------------------------------------------------------
  !!                               DATA TYPES
  !!--------------------------------------------------------------------------

  !>
  !! Derived data type: t_radiation_forcing
  !!
  !! This structure contains components:
  !! <ol>
  !! <li> related to the instantaneous aerosol radiative forcing
  !! </ol>
  !!
  !! All components are arrays of one of the following shapes:
  !! <ol>
  !! <li> (nproma,           nblks_phy)
  !! <li> (nproma, nlev_phy(+1), nblks_phy)
  !! </ol>
  !! Currently the physics grid has the same spatial resolution as the
  !! dynamics grid, but is unstaggered. This means
  !!
  !!    nlev_phy = nlev
  !!   nblks_phy = patch%nblks_c
  !!
  !! In the long run, the physics grid and dynamics grid may differ in
  !! horizontal and/or vertical resolution, or even in shape.

  TYPE t_radiation_forcing
      !--- Auxiliary radiative fluxes:

  REAL(wp),  POINTER :: emter_for(:,:,:) !<Net dwnwrd LW flux [Wm2]
  REAL(wp),  POINTER :: emtef_for(:,:,:) !<Net dwnnrd LW flux (clear sky) [Wm2]
  REAL(wp),  POINTER :: trsol_for(:,:,:) !<Net dwnwrd SW flux [Wm2] for forcing
  REAL(wp),  POINTER :: trsof_for(:,:,:) !<Net dwnwrd SW flux (clear sky) [Wm2]

  !--- Forcing fields TOA, SUR:

  REAL(wp),  POINTER :: fsw_clear_top(:,:) !<instantaneous sw forcing clear sky top of atmosphere
  REAL(wp),  POINTER :: fsw_total_top(:,:) !<instantaneous sw forcing all sky top of atmosphere
  REAL(wp),  POINTER :: fsw_clear_sur(:,:) !<instantaneous sw forcing clear sky surface
  REAL(wp),  POINTER :: fsw_total_sur(:,:) !<instantaneous sw forcing all sky surface

  REAL(wp),  POINTER :: flw_clear_top(:,:) !<instantaneous lw forcing clear sky top of atmosphere
  REAL(wp),  POINTER :: flw_total_top(:,:) !<instantaneous lw forcing all sky top of atmosphere
  REAL(wp),  POINTER :: flw_clear_sur(:,:) !<instantaneous lw forcing clear sky surface
  REAL(wp),  POINTER :: flw_total_sur(:,:) !<instantaneous lw forcing all sky surface

  !--- Radiation flux forcing:
  REAL(wp),  POINTER :: d_aflx_sw(:,:,:)  !<3d instantaneous sw forcing all sky
  REAL(wp),  POINTER :: d_aflx_lw(:,:,:)  !<3d instantaneous lw forcing all sky
  REAL(wp),  POINTER :: d_aflx_swc(:,:,:) !<3d instantaneous sw forcing clear sky
  REAL(wp),  POINTER :: d_aflx_lwc(:,:,:) !<3d instantaneous lw forcing clear sky

  !--- Heating rate forcing:
  REAL(wp), POINTER  :: netht_lw(:,:,:)   !<3d forcing of net lw heating rate (K/day) 
  REAL(wp), POINTER  :: netht_sw(:,:,:)   !<3d forcing of net sw heating rate (K/day)

  END TYPE t_radiation_forcing

  !!--------------------------------------------------------------------------
  !!                          diagnostic variables 
  !!--------------------------------------------------------------------------
  !! The variable names have the prefix "prm_" in order to emphasize that they
  !! are defined for and used in parameterisations of radiation.

  TYPE(t_radiation_forcing),ALLOCATABLE,TARGET :: prm_radiation_forcing(:)  !< shape: (n_dom)

  !!--------------------------------------------------------------------------
  !!                          variable lists
  !!--------------------------------------------------------------------------
  TYPE(t_var_list_ptr),ALLOCATABLE :: prm_radiation_forcing_list(:)  !< shape: (n_dom)

CONTAINS

  !!--------------------------------------------------------------------------
  !!                SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
  !!--------------------------------------------------------------------------
  !>
  !! Top-level procedure for building the physics state
  !!
  SUBROUTINE construct_radiation_forcing_list( patch_array )

    TYPE(t_patch),INTENT(IN) :: patch_array(:)
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname
    INTEGER :: ndomain, jg, ist, nblks, nlev
    CHARACTER(len=*), PARAMETER :: thissubprog='construct_radiation_forcing_list of mo_radiation_forcing_memory'
    
    !---

    IF (.NOT.(lradforcing(1).OR.lradforcing(2))) RETURN
    CALL message(TRIM(thissubprog),'Construction of radiation_forcing_list started.')

    ! Allocate pointer arrays prm_field and prm_tend, 
    ! as well as the corresponding list arrays.

    ndomain = SIZE(patch_array)

    ALLOCATE( prm_radiation_forcing(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(TRIM(thissubprog), &
      &'allocation of prm_radiation_forcing array failed')

    ALLOCATE( prm_radiation_forcing_list(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(TRIM(thissubprog), &
      &'allocation of prm_radiation_forcing list array failed')

    !$ACC ENTER DATA CREATE(prm_radiation_forcing)

    ! Build a rte-rrtmgp forcing list for each grid level.
    ! This includes memory allocation. 
    
    DO jg = 1,ndomain

      nblks = patch_array(jg)%nblks_c
      nlev  = patch_array(jg)%nlev

      WRITE(listname,'(a,i2.2)') 'prm_radiation_forcing_D',jg
      CALL new_radiation_forcing_list( jg, nproma, nlev, nblks,          &
                                   & TRIM(listname), '',             &
                                   & prm_radiation_forcing_list(jg),     &
                                   & prm_radiation_forcing(jg)          )

      CALL message(TRIM(thissubprog),'Construction of rte-rrtmgp forcing list finished.')

   END DO

  END SUBROUTINE construct_radiation_forcing_list
  !--------------------------------------------------------------------
  !>
  !! Release memory used by the radiation_forcing arrays and list arrays
  !!
  SUBROUTINE destruct_radiation_forcing_list

    INTEGER :: ndomain  !< total # of grid levels/domains
    INTEGER :: jg       !< grid level/domain index
    INTEGER :: ist      !< system status code
    CHARACTER(len=*), PARAMETER :: thissubprog='destruct_radiation_forcing_list of mo_radiation_forcing_memory'
    !---
    IF (.NOT.(lradforcing(1).OR.lradforcing(2))) RETURN
    CALL message(TRIM(thissubprog),'Destruction of radiation_forcing_list started.')

    ndomain = SIZE(prm_radiation_forcing)

    DO jg = 1,ndomain
      CALL vlr_del( prm_radiation_forcing_list(jg) )
    ENDDO

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(prm_radiation_forcing)

    DEALLOCATE( prm_radiation_forcing_list, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thissubprog), &
      & 'deallocation of radiation_forcing list array failed')

    DEALLOCATE( prm_radiation_forcing, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thissubprog), &
      & 'deallocation of radiation_forcing array failed')

    CALL message(TRIM(thissubprog),'Destruction of radiation_forcing_list finished.')

  END SUBROUTINE destruct_radiation_forcing_list
  !--------------------------------------------------------------------

  SUBROUTINE new_radiation_forcing_list( k_jg, kproma, klev, kblks,          &
                                   & listname, prefix,             &
                                   & field_list,     &
                                   & field          )

    INTEGER,INTENT(IN) :: k_jg !> patch ID
    INTEGER,INTENT(IN) :: kproma, klev, kblks !< dimension sizes

    CHARACTER(len=*)              ,INTENT(IN) :: listname, prefix

    TYPE(t_var_list_ptr),     INTENT(INOUT)   :: field_list
    TYPE(t_radiation_forcing),INTENT(INOUT)   :: field

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shape3d_layer_interfaces(3)
!!$    INTEGER :: shape4d(4)
    INTEGER :: ibits
    INTEGER :: datatype_flt
!!$    INTEGER :: jsfc, jtrc

    ibits = DATATYPE_PACK16

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    shape2d  = (/kproma,       kblks/)
    shape3d  = (/kproma, klev, kblks/)
    shape3d_layer_interfaces = (/kproma,klev+1,kblks/)


    ! Register a field list and apply default settings

    CALL vlr_add( field_list, TRIM(listname), patch_id=k_jg, lrestart=.TRUE., &
      &           model_type=get_my_process_name() )

    ! Auxiliary flux variables
    IF (lradforcing(2)) THEN
    cf_desc    = t_cf_var('emter_for', 'W m-2', 'thermal radiation flux', datatype_flt)
    grib2_desc = grib2_var(0, 5, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'emter_for', field%emter_for,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,                &
                & ldims=shape3d_layer_interfaces,                                             &
                & vert_interp =                                                               &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ) )

    cf_desc    = t_cf_var('emter_for', 'W m-2', 'thermal radiation flux', datatype_flt)
    grib2_desc = grib2_var(0, 5, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'emtef_for', field%emtef_for,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,                &
                & ldims=shape3d_layer_interfaces,                                             &
                & vert_interp =                                                               &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ) )
    END IF

    IF (lradforcing(1)) THEN
    cf_desc    = t_cf_var('trsol_for', 'W m-2', 'solar radiation flux', datatype_flt)
    grib2_desc = grib2_var(0, 4, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'trsol_for', field%trsol_for,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,                &
                & ldims=shape3d_layer_interfaces,                                             &
                & vert_interp =                                                               &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ) )

    cf_desc    = t_cf_var('trsof_for', 'W m-2', 'solar radiation flux', datatype_flt)
    grib2_desc = grib2_var(0, 4, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'trsof_for', field%trsof_for,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,                &
                & ldims=shape3d_layer_interfaces,                                             &
                & vert_interp =                                                               &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ) )
    END IF

    ! diagnostic fields

    IF (lradforcing(1)) THEN
    cf_desc    = t_cf_var('fsw_clear_top', 'W m-2',                                        &
                & 'instantaneous sw forcing clear sky top of atmosphere', datatype_flt)
    grib2_desc = grib2_var(0, 4, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'fsw_clear_top', field%fsw_clear_top,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,  &
                & lrestart=.FALSE. )

    cf_desc    = t_cf_var('fsw_total_top', 'W m-2',                                        &
                & 'instantaneous sw forcing all sky top of atmosphere', datatype_flt)
    grib2_desc = grib2_var(0, 4, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'fsw_total_top', field%fsw_total_top,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,  &
                & lrestart=.FALSE. )

    cf_desc    = t_cf_var('fsw_clear_sur', 'W m-2',                                        &
                & 'instantaneous sw forcing clear sky surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'fsw_clear_sur', field%fsw_clear_sur,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,  &
                & lrestart=.FALSE. )

    cf_desc    = t_cf_var('fsw_total_sur', 'W m-2',                                        &
                & 'instantaneous sw forcing all sky surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'fsw_total_sur', field%fsw_total_sur,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,  &
                & lrestart=.FALSE. )
    END IF

    IF (lradforcing(2)) THEN
    cf_desc    = t_cf_var('flw_clear_top', 'W m-2',                                        &
                & 'instantaneous lw forcing clear sky top of atmosphere', datatype_flt)
    grib2_desc = grib2_var(0, 5, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'flw_clear_top', field%flw_clear_top,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,  &
                & lrestart=.FALSE. )

    cf_desc    = t_cf_var('flw_total_top', 'W m-2',                                        &
                & 'instantaneous lw forcing all sky top of atmosphere', datatype_flt)
    grib2_desc = grib2_var(0, 5, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'flw_total_top', field%flw_total_top,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,  &
                & lrestart=.FALSE. )

    cf_desc    = t_cf_var('flw_clear_sur', 'W m-2',                                        &
                & 'instantaneous lw forcing clear sky surface', datatype_flt)
    grib2_desc = grib2_var(0, 5, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'flw_clear_sur', field%flw_clear_sur,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,  &
                & lrestart=.FALSE. )

    cf_desc    = t_cf_var('flw_total_sur', 'W m-2',                                        &
                & 'instantaneous lw forcing all sky surface', datatype_flt)
    grib2_desc = grib2_var(0, 5, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'flw_total_sur', field%flw_total_sur,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,  &
                & lrestart=.FALSE. )
    END IF
    
    IF (lradforcing(1)) THEN
    cf_desc    = t_cf_var('d_aflx_sw', 'W m-2', '3d instantaneous sw forcing all sky',       &
               & datatype_flt)
    grib2_desc = grib2_var(0, 4, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'d_aflx_sw', field%d_aflx_sw,                          &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,               &
                ldims=shape3d_layer_interfaces,                                              &
                & vert_interp =                                                              &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ),  &
                & lrestart=.FALSE. )

    cf_desc    = t_cf_var('d_aflx_swc', 'W m-2', '3d instantaneous sw forcing clear sky',    &
               & datatype_flt)
    grib2_desc = grib2_var(0, 4, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'d_aflx_swc', field%d_aflx_swc,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,               &
                & ldims=shape3d_layer_interfaces,                                            &
                & vert_interp =                                                              &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ),  &
                & lrestart=.FALSE. )
    END IF

    IF (lradforcing(2)) THEN
    cf_desc    = t_cf_var('d_aflx_lw', 'W m-2', '3d instantaneous lw forcing all sky',       &
               & datatype_flt)
    grib2_desc = grib2_var(0, 5, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'d_aflx_lw', field%d_aflx_lw,                          &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,               &
                & ldims=shape3d_layer_interfaces,                                            &
                & vert_interp =                                                              &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ),  &
                & lrestart=.FALSE. )

    cf_desc    = t_cf_var('d_aflx_lwc', 'W m-2', '3d instantaneous lw forcing clear sky',    &
               & datatype_flt)
    grib2_desc = grib2_var(0, 5, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'d_aflx_lwc', field%d_aflx_lwc,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,               &
                & ldims=shape3d_layer_interfaces,                                            &
                & vert_interp =                                                              &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ),  &
                & lrestart=.FALSE. )

    cf_desc    = t_cf_var('netht_lw', 'K d-1', '3d forcing of net lw heating rate',          &
               & datatype_flt)
    grib2_desc = grib2_var(0, 0, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'netht_lw', field%netht_lw,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d,     & 
                & vert_interp =                                                              &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ),  &
                & lrestart=.FALSE. )
    END IF

    IF (lradforcing(1)) THEN
    cf_desc    = t_cf_var('netht_sw', 'K d-1', '3d forcing of net sw heating rate',          &
               & datatype_flt)
    grib2_desc = grib2_var(0, 0, 23, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'netht_sw', field%netht_sw,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d,     &
                & vert_interp =                                                              &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ),  &
                & lrestart=.FALSE. )
    END IF

  END SUBROUTINE new_radiation_forcing_list


END MODULE mo_radiation_forcing_memory
