! Contains basic diagnostics for ICON ocean model.
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
MODULE mo_ocean_check_total_content
  USE mo_master_control,     ONLY: get_my_process_name
  USE mo_kind,               ONLY: wp, dp, i8
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range, t_subset_indexed
  USE mo_grid_tools,         ONLY: get_oriented_edges_from_global_vertices, check_global_indexes
  USE mo_mpi,                ONLY: my_process_is_stdio, p_field_sum, get_my_mpi_work_id, &
    & p_comm_work_test, p_comm_work, p_io, p_bcast, my_process_is_mpi_workroot, p_sum, &
    & my_process_is_mpi_parallel
  USE mo_sync,               ONLY: global_sum_array, disable_sync_checks, enable_sync_checks, &
    &                              sync_c, sync_e, sync_patch_array
  USE mo_math_types,         ONLY: t_cartesian_coordinates
  USE mo_math_utilities,     ONLY: cvec2gvec
  USE mo_advection_utils,    ONLY: laxfr_upflux
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_dbg_nml,            ONLY: idbg_val
  USE mo_math_constants,     ONLY: rad2deg, dbl_eps
  USE mo_impl_constants,     ONLY: sea_boundary,sea, &
    & min_rlcell, min_rledge, min_rlcell, &
    & max_char_length, min_dolic
   USE mo_cdi_constants,      ONLY: GRID_EDGE, GRID_CELL, GRID_UNSTRUCTURED_EDGE, &
    & GRID_UNSTRUCTURED_CELL
   USE mo_sea_ice_nml,        ONLY: kice, sice
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_parallel_config,    ONLY: nproma, p_test_run
  USE mo_run_config,         ONLY: dtime, nsteps
  USE mo_physical_constants, ONLY: grav, rhos, rhoi, rho_ref, clw, alf
  USE mo_model_domain,       ONLY: t_patch, t_patch_3d,t_patch_vert, t_grid_edges
  USE mo_ocean_types,        ONLY: t_hydro_ocean_state, t_hydro_ocean_diag
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_exception,          ONLY: message, finish, message_text, warning
  USE mo_sea_ice_types,      ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_ocean_surface_types,ONLY: t_ocean_surface
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_scalar_product,     ONLY: map_edges2cell_3d
  USE mo_io_units,           ONLY: find_next_free_unit
  USE mo_util_file,          ONLY: util_symlink, util_rename, util_islink, util_unlink
  USE mo_statistics,         ONLY: subset_sum, levels_horizontal_mean, total_mean, gather_sums, &
    & verticallyIntegrated_field
  USE mo_var_groups,         ONLY: groups
  USE mo_cf_convention
  USE mo_grib2,              ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, GRID_UNSTRUCTURED
  USE mo_cdi_constants,      ONLY: GRID_EDGE, GRID_CELL, GRID_UNSTRUCTURED_EDGE, &
    &                              GRID_UNSTRUCTURED_CELL
  USE mo_zaxis_type,         ONLY: ZA_DEPTH_BELOW_SEA
  USE mo_io_config,          ONLY: lnetcdf_flt64_output
  USE mo_name_list_output_init, ONLY: isRegistered

  USE mtime,                 ONLY: datetime, MAX_DATETIME_STR_LEN, datetimeToPosixString

  USE mo_ocean_nml,          ONLY: n_zlev, no_tracer, vert_cor_type
  USE mo_fortran_tools,      ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_total_salt_content, check_total_salt_content, check_total_si_volume
  PUBLIC :: check_total_salt_content_zstar, calc_total_salt_content_zstar
  PUBLIC :: check_accumulated_volume_difference


  REAL(wp) :: total_salt_old=1.0_wp
  REAL(wp) :: total_saltinseaice_old=1.0_wp
  REAL(wp) :: total_saltinliquidwater_old=1.0_wp
  REAL(wp) :: initial_total_salt = 0.0_wp
  REAL(wp) :: accumulated_run_error(500)

  REAL(wp) :: total_watvolume_old = 1.0_wp
  REAL(wp) :: total_icevolume_old = 1.0_wp
  REAL(wp) :: total_snowvolume_old = 1.0_wp
  REAL(wp) :: total_icearea_old = 1.0_wp

  REAL(wp) :: accumulated_volume_difference(256)  = 0.0_wp
  REAL(wp) :: accumulated_h_difference(256)       = 0.0_wp


CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE check_accumulated_volume_difference(id, patch_2d, h_new, h_old)

    INTEGER, INTENT(IN)                                              :: id
    TYPE(t_patch), TARGET, INTENT(IN)                                :: patch_2d
    REAL(wp),DIMENSION(nproma,patch_2d%alloc_cell_blocks),INTENT(IN) :: h_new, h_old

    REAL(wp)                                               :: volume_diff, h_mean_diff
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: h_diff
    TYPE(t_subset_range), POINTER :: owned_cells
    INTEGER :: blockNo
 
    owned_cells => patch_2D%cells%owned

!ICON_OMP_PARALLEL_DO ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = owned_cells%start_block, owned_cells%end_block
      h_diff(:,blockNo) = h_new(:,blockNo) - h_old(:,blockNo)
    ENDDO
    
    volume_diff = subset_sum(h_diff, patch_2d%cells%area, owned_cells, h_mean_diff)
    
    accumulated_volume_difference(id) = accumulated_volume_difference(id) + volume_diff
    accumulated_h_difference(id)      = accumulated_h_difference(id)      + h_mean_diff
    
    IF (my_process_is_stdio()) THEN
      WRITE(0,*) '(',id,')', ' -- mean h error, accumulated:', h_mean_diff, accumulated_h_difference(id)
      WRITE(0,*) '(',id,')', ' -- volume error, accumulated:', volume_diff, accumulated_volume_difference(id)
    ENDIF
  END SUBROUTINE check_accumulated_volume_difference
  !-------------------------------------------------------------------------
   
  !-------------------------------------------------------------------------
  SUBROUTINE check_total_salt_content(id, so, patch_2d, h, thickness, ice, im)

    TYPE(t_patch), TARGET, INTENT(IN)                                 :: patch_2d
    REAL(wp),DIMENSION(nproma,patch_2d%alloc_cell_blocks),INTENT(IN) :: h
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: thickness
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: so !salinity
    TYPE (t_sea_ice),       INTENT(IN)                    :: ice

    INTEGER :: id
    INTEGER :: im

    REAL(wp)                                              :: total_salt
    REAL(wp)                                              :: total_saltinseaice
    REAL(wp)                                              :: total_saltinliquidwater

    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: salt
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinseaice
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinliquidwater


    CALL calc_salt_content(so, patch_2d, h, thickness, ice, im,  &
                           salt, saltinseaice, saltinliquidwater )

    total_salt = global_sum_array(salt)
    total_saltinseaice = global_sum_array(saltinseaice)
    total_saltinliquidwater = global_sum_array(saltinliquidwater)
    IF (initial_total_salt == 0.0_wp) THEN 
      initial_total_salt = total_salt
      total_salt_old = total_salt
      accumulated_run_error(:) = 0.0_wp
    ELSE
      accumulated_run_error(id) = accumulated_run_error(id) + total_salt - total_salt_old
    ENDIF
!     IF (my_process_is_stdio()) THEN
!       WRITE(0,*) '(',id,')', ' salt   :' , total_salt, (total_salt / total_salt_old) - 1.0_wp
!       WRITE(0,*) '(',id,')', ' saltice:' , total_saltinseaice, ( total_saltinseaice / max(1.0e-19_wp,total_saltinseaice_old)) - 1.0_wp
!       WRITE(0,*) '(',id,')', ' saltwat:' , total_saltinliquidwater, (total_saltinliquidwater / total_saltinliquidwater_old) - 1.0_wp
!     ENDIF

     IF (my_process_is_stdio()) THEN
       WRITE(0,*) '(',id,')', ' salt   :' , total_salt, total_salt - total_salt_old
!       WRITE(0,*) '(',id,')', ' saltice:' , total_saltinseaice, total_saltinseaice - total_saltinseaice_old
!       WRITE(0,*) '(',id,')', ' saltwat:' , total_saltinliquidwater, total_saltinliquidwater - total_saltinliquidwater_old
     ENDIF
   
!    IF (my_process_is_stdio()) THEN
!       WRITE(0,*) '(',id,')', ' -- total init error:', (total_salt - initial_total_salt) / initial_total_salt 
!       WRITE(0,*) '(',id,')', ' -- accumulated run error:', accumulated_run_error(id) /  total_salt
!    ENDIF
   
    total_salt_old = total_salt
    total_saltinseaice_old = total_saltinseaice
    total_saltinliquidwater_old = total_saltinliquidwater

  END SUBROUTINE check_total_salt_content
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
 SUBROUTINE calc_total_salt_content(so, patch_2d, h, thickness, ice, im, &
                  total_salt, total_saltinseaice, total_saltinliquidwater, lacc)

    TYPE(t_patch), TARGET, INTENT(IN)                                 :: patch_2d
    REAL(wp),DIMENSION(nproma,patch_2d%alloc_cell_blocks),INTENT(IN) :: h
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: thickness
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: so !salinity
    TYPE (t_sea_ice),       INTENT(IN)                    :: ice

    REAL(wp)                                              :: total_salt
    REAL(wp)                                              :: total_saltinseaice
    REAL(wp)                                              :: total_saltinliquidwater
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: salt
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinseaice
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinliquidwater

    INTEGER :: im 
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPY(salt, saltinseaice, saltinliquidwater) IF(lzacc)

    CALL calc_salt_content(so, patch_2d, h, thickness, ice, im , &
                           salt, saltinseaice, saltinliquidwater, lacc=lzacc )

    !$ACC END DATA

    total_salt = global_sum_array(salt)
    total_saltinseaice = global_sum_array(saltinseaice)
    total_saltinliquidwater = global_sum_array(saltinliquidwater)

  END SUBROUTINE calc_total_salt_content
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE calc_salt_content(so, patch_2d, h, thickness, ice, im, &
            salt, saltinseaice, saltinliquidwater, lacc)
    TYPE(t_patch), TARGET, INTENT(IN)                                 :: patch_2d
    REAL(wp),DIMENSION(nproma,patch_2d%alloc_cell_blocks),INTENT(IN) :: h
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: thickness
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: so !salinity
    TYPE (t_sea_ice),       INTENT(IN)                    :: ice
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! locals
    REAL(wp), INTENT(INOUT), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: salt
    REAL(wp), INTENT(INOUT), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinseaice
    REAL(wp), INTENT(INOUT), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinliquidwater

    REAL(wp) :: rhoicwa,rhosnwa,draftave
    LOGICAL :: lzacc

 
    TYPE(t_subset_range), POINTER                         :: subset
    INTEGER                                               :: block, cell, cellStart,cellEnd, level
    INTEGER                                               :: im

    CALL set_acc_host_or_device(lzacc, lacc)

    IF(no_tracer<=1)RETURN

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    salt(:,:)              = 0.0_wp
    saltinseaice(:,:)      = 0.0_wp
    saltinliquidwater(:,:) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)

    rhoicwa = rhoi / rho_ref
    rhosnwa = rhos / rho_ref

    subset => patch_2d%cells%owned

    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, cellStart, cellEnd)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO cell = cellStart, cellEnd
        IF (subset%vertical_levels(cell,block) < 1) CYCLE

        ! surface:
        saltInSeaice(cell,BLOCK)    = sice*rhoicwa &
             &                         * SUM(ice%hi(cell,:,BLOCK)*ice%conc(cell,:,BLOCK)) &
             &                         * patch_2d%cells%area(cell,BLOCK)

        IF (im .eq. 0) then 

          draftave=sum((rhoicwa*ice%hi(cell,:,BLOCK)+rhosnwa*ice%hs(cell,:,BLOCK))&
                    *ice%conc(cell,:,BLOCK))

        else

          draftave = ice%draftave_old(cell,BLOCK)

        endif

        IF ( vert_cor_type .EQ. 1) THEN
          draftave = 0.0_wp
        END IF

        saltInLiquidWater(cell,BLOCK) = so(cell,1,BLOCK) &
             &                    * (thickness(cell,1,BLOCK)+h(cell,BLOCK)-draftave) &
             &                    * patch_2d%cells%area(cell,BLOCK)

        !$ACC LOOP SEQ
        DO level=2,subset%vertical_levels(cell,BLOCK)
          saltinliquidwater(cell,BLOCK) = saltinliquidwater(cell,BLOCK) + so(cell,level,BLOCK) &
            &                    * thickness(cell,level,block) &
            &                    * patch_2d%cells%area(cell,block)
        END DO

        salt(cell,block) = saltInSeaice(cell,block) + saltInLiquidWater(cell,block)

        ! rest of the underwater world
      END DO ! cell
      !$ACC END PARALLEL LOOP
    END DO !block
    !$ACC WAIT(1)

  END SUBROUTINE calc_salt_content
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE calc_si_volume(patch_2d, ice, icevolume, snowvolume, icearea)

    TYPE(t_patch), TARGET, INTENT(IN)                      :: patch_2d
    TYPE (t_sea_ice),      INTENT(IN)                      :: ice

    ! locals
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: icevolume
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: snowvolume
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: icearea

 
    TYPE(t_subset_range), POINTER                          :: subset
    INTEGER                                                :: block, cell, cellStart,cellEnd, level

    IF(no_tracer<=1)RETURN

    icevolume  = 0.0_wp
    snowvolume = 0.0_wp
    icearea = 0.0_wp

    subset => patch_2d%cells%owned
    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, cellStart, cellEnd)
      DO cell = cellStart, cellEnd
        IF (subset%vertical_levels(cell,block) < 1) CYCLE

        ! icevolume
        icevolume(cell,BLOCK)    = SUM(ice%hi(cell,:,BLOCK)*ice%conc(cell,:,BLOCK)) &
             * patch_2d%cells%area(cell,BLOCK)

        ! snowvolume
        snowvolume(cell,BLOCK)    = SUM(ice%hs(cell,:,BLOCK)*ice%conc(cell,:,BLOCK)) &
             * patch_2d%cells%area(cell,BLOCK)

        ! icearea
        icearea(cell,BLOCK)    = SUM(ice%conc(cell,:,BLOCK)) &
             * patch_2d%cells%area(cell,BLOCK)


      END DO ! cell
    END DO !block
  END SUBROUTINE calc_si_volume
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_total_si_volume(id,  patch_2d, ice, h)

    TYPE(t_patch), TARGET, INTENT(IN)                                 :: patch_2d
    TYPE (t_sea_ice),       INTENT(IN)                    :: ice

    REAL(wp),DIMENSION(nproma,patch_2d%alloc_cell_blocks),INTENT(IN) :: h

    INTEGER :: id

    REAL(wp)                                              :: total_watvolume
    REAL(wp)                                              :: total_icevolume
    REAL(wp)                                              :: total_snowvolume
    REAL(wp)                                              :: total_icearea

    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: icevolume
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: snowvolume
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: icearea

    CALL calc_si_volume( patch_2d, ice,  &
                           icevolume, snowvolume, icearea )
    total_watvolume = global_sum_array(h)
    total_icevolume = global_sum_array(icevolume)
    total_snowvolume = global_sum_array(snowvolume)
    total_icearea = global_sum_array(icearea)

     IF (my_process_is_stdio()) THEN
       WRITE(0,*) '(',id,')', ' watvol:' , total_watvolume, (total_watvolume / MAX(1.0e-19_wp,total_watvolume_old)) - 1.0_wp
       WRITE(0,*) '(',id,')', ' icevol:' , total_icevolume, (total_icevolume / MAX(1.0e-19_wp,total_icevolume_old)) - 1.0_wp
       WRITE(0,*) '(',id,')', ' snovol:' , total_snowvolume, ( total_snowvolume / MAX(1.0e-19_wp,total_snowvolume_old)) - 1.0_wp
!       WRITE(0,*) '(',id,')', ' icearea:' , total_icearea, (total_icearea / MAX(1.0e-19_wp, total_icearea_old)) - 1.0_wp

     ENDIF

    total_watvolume_old = total_watvolume
    total_icevolume_old = total_icevolume
    total_snowvolume_old = total_snowvolume
    total_icearea_old = total_icearea

  END SUBROUTINE check_total_si_volume
  !-------------------------------------------------------------------------

   
   SUBROUTINE check_total_salt_content_zstar(id, so, patch_2d, stretch, thickness, ice, p_oce_sfc)

    TYPE(t_patch), TARGET, INTENT(IN)                                 :: patch_2d
    REAL(wp),DIMENSION(nproma,patch_2d%alloc_cell_blocks),INTENT(IN)  :: stretch
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: thickness
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: so !salinity
    TYPE (t_sea_ice),       INTENT(IN)                    :: ice

    INTEGER :: id

    REAL(wp)                                              :: total_salt
    REAL(wp)                                              :: total_saltinseaice
    REAL(wp)                                              :: total_saltinliquidwater

    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: salt
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinseaice
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinliquidwater

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! REMOVE
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    TYPE(t_subset_range), POINTER                :: subset
    INTEGER                                      :: block, cell, cellStart,cellEnd, level
    REAL(wp)  :: flux, flux_tot 

    flux     = 0.0_wp
    flux_tot = 0.0_wp

    subset => patch_2d%cells%owned
    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, cellStart, cellEnd)
      DO cell = cellStart, cellEnd
        flux = flux + patch_2d%cells%area(cell,BLOCK) *p_oce_sfc%FrshFlux_IceSalt(cell, block) * dtime
      END DO ! cell
    END DO !block
    
    flux_tot = global_sum_array(flux)
!    IF (my_process_is_stdio()) THEN
!       WRITE(0,*) '(',id,')', ' flux   :', flux_tot 
!     ENDIF
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL calc_salt_content_zstar(so, patch_2d, stretch, thickness, ice, &
                           salt, saltinseaice, saltinliquidwater )

    total_salt = global_sum_array(salt)
    total_saltinseaice = global_sum_array(saltinseaice)
    total_saltinliquidwater = global_sum_array(saltinliquidwater)
    IF (initial_total_salt == 0.0_wp) THEN 
      initial_total_salt = total_salt
      total_salt_old = total_salt
      accumulated_run_error(:) = 0.0_wp
    ELSE
      accumulated_run_error(id) = accumulated_run_error(id) + total_salt - total_salt_old
    ENDIF

     IF (my_process_is_stdio()) THEN
       WRITE(0,*) '(',id,')', ' salt   :' , total_salt, total_salt - total_salt_old
!       WRITE(0,*) '(',id,')', ' saltice:' , total_saltinseaice, total_saltinseaice - total_saltinseaice_old
!       WRITE(0,*) '(',id,')', ' saltwat:' , total_saltinliquidwater, total_saltinliquidwater - total_saltinliquidwater_old
!       WRITE(0,*) '(',id,')', ' saltwat:' , total_saltinliquidwater, total_saltinseaice
     ENDIF
   
!    IF (my_process_is_stdio()) THEN
!       WRITE(0,*) '(',id,')', ' -- total init error:', (total_salt - initial_total_salt) / initial_total_salt 
!       WRITE(0,*) '(',id,')', ' -- accumulated run error:', accumulated_run_error(id) /  total_salt
!    ENDIF
   
    total_salt_old = total_salt
    total_saltinseaice_old = total_saltinseaice
    total_saltinliquidwater_old = total_saltinliquidwater



  END SUBROUTINE check_total_salt_content_zstar

  
  SUBROUTINE calc_total_salt_content_zstar(so, patch_2d, stretch, thickness, ice, &
      & total_salt, total_saltinseaice, total_saltinliquidwater, lacc)

    TYPE(t_patch), TARGET, INTENT(IN)                                 :: patch_2d
    REAL(wp),DIMENSION(nproma,patch_2d%alloc_cell_blocks),INTENT(IN)  :: stretch
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: thickness
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: so !salinity
    TYPE (t_sea_ice),       INTENT(IN)                    :: ice

    REAL(wp), INTENT(OUT)                                 :: total_salt
    REAL(wp), INTENT(OUT)                                 :: total_saltinseaice
    REAL(wp), INTENT(OUT)                                 :: total_saltinliquidwater
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: salt
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinseaice
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinliquidwater
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPY(salt, saltinseaice, saltinliquidwater) IF(lzacc)

    CALL calc_salt_content_zstar(so, patch_2d, stretch, thickness, ice, &
                           salt, saltinseaice, saltinliquidwater, lacc=lzacc )

    !$ACC END DATA

    total_salt = global_sum_array(salt)
    total_saltinseaice = global_sum_array(saltinseaice)
    total_saltinliquidwater = global_sum_array(saltinliquidwater)

  END SUBROUTINE calc_total_salt_content_zstar
 
  
  SUBROUTINE calc_salt_content_zstar(so, patch_2d, stretch, thickness, ice, &
            salt, saltinseaice, saltinliquidwater, lacc)
    TYPE(t_patch), TARGET, INTENT(IN)                                 :: patch_2d
    REAL(wp),DIMENSION(nproma,patch_2d%alloc_cell_blocks),INTENT(IN)  :: stretch
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: thickness
    REAL(wp),DIMENSION(nproma,n_zlev,patch_2d%alloc_cell_blocks),INTENT(IN) :: so !salinity
    TYPE (t_sea_ice),       INTENT(IN)                    :: ice
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! locals
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: salt
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinseaice
    REAL(wp), DIMENSION(nproma,patch_2d%alloc_cell_blocks) :: saltinliquidwater

    REAL(wp) :: rhoicwa,rhosnwa,draftave
    LOGICAL :: lzacc

 
    TYPE(t_subset_range), POINTER                         :: subset
    INTEGER                                               :: block, cell, cellStart,cellEnd, level

    CALL set_acc_host_or_device(lzacc, lacc)

    IF(no_tracer<=1)RETURN

    salt         = 0.0_wp
    saltinseaice = 0.0_wp
    saltinliquidwater = 0.0_wp

    rhoicwa = rhoi / rho_ref
    rhosnwa = rhos / rho_ref

    subset => patch_2d%cells%owned

    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, cellStart, cellEnd)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO cell = cellStart, cellEnd
        IF (subset%vertical_levels(cell,block) < 1) CYCLE

        ! surface:
        saltInSeaice(cell,BLOCK)    = sice*rhoicwa &
             &                         * SUM(ice%hi(cell,:,BLOCK)*ice%conc(cell,:,BLOCK)) &
             &                         * patch_2d%cells%area(cell,BLOCK)

        !$ACC LOOP SEQ
        DO level=1,subset%vertical_levels(cell,BLOCK)
          saltinliquidwater(cell,BLOCK) = saltinliquidwater(cell,BLOCK) + so(cell,level,BLOCK) &
            &   * stretch(cell, block) * thickness(cell,level,block) &
            &                          * patch_2d%cells%area(cell,block)
        END DO

        salt(cell,block) = saltInSeaice(cell,block) + saltInLiquidWater(cell,block)

        ! rest of the underwater world
      END DO ! cell
      !$ACC END PARALLEL LOOP
    END DO !block
    !$ACC WAIT(1)

  END SUBROUTINE calc_salt_content_zstar



END MODULE mo_ocean_check_total_content
