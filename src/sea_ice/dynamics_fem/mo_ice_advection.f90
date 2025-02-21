! Contains sea ice advection routines (on ICON grid).
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

MODULE mo_ice_advection
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_exception,           ONLY: message, message_text, finish  
  USE mo_run_config,          ONLY: dtime
  USE mo_parallel_config,     ONLY: nproma
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, sync_patch_array, sync_patch_array_mult
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_impl_constants,      ONLY: sea_boundary

  USE mo_advection_utils,     ONLY: laxfr_upflux
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_ocean_math_operators,  ONLY: div_oce_3D

  USE mo_sea_ice_types,       ONLY: t_sea_ice
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PUBLIC  :: ice_advection_upwind
  PUBLIC  :: ice_advection_upwind_einar

  PRIVATE :: upwind_hflux_ice

  CHARACTER(len=12)           :: str_module    = 'IceAdvect'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1            ! Level of detail for 1 line debug

CONTAINS

  !-------------------------------------------------------------------------
  !
  !> Advection of sea ice and snow on ice.
  !> Modified to advect ice%hi*ice%conc and conc instead of ice%vol.
  !! This uses the upwind_hflux_ice routine and the ocean's div_oce_3D routine to do upwind
  !! advection of the relevant variables.
  !!
  SUBROUTINE ice_advection_upwind( p_patch_3D, p_op_coeff, p_ice, lacc )
    TYPE(t_patch_3D), TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_operator_coeff),   INTENT(IN)    :: p_op_coeff
    TYPE(t_sea_ice),          INTENT(INOUT) :: p_ice
    LOGICAL, INTENT(IN), OPTIONAL           :: lacc

    ! Local variables
    ! Patch and range
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_subset_range), POINTER :: cells_in_domain

    ! Indexing
    INTEGER  :: jk, jb, jc
    INTEGER  :: i_startidx_c, i_endidx_c
    LOGICAL  :: lzacc

    ! Temporary variables/buffers
    REAL(wp) :: z_adv_flux_h (nproma,p_ice%kice,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: flux_hi  (nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: flux_conc(nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: flux_hs  (nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: tmp(SIZE(p_ice%hi,1), SIZE(p_ice%hi,2), SIZE(p_ice%hi,3))

    CALL set_acc_host_or_device(lzacc, lacc)

!--------------------------------------------------------------------------------------------------
    patch_2D => p_patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
!--------------------------------------------------------------------------------------------------

    !$ACC DATA CREATE(tmp, z_adv_flux_h, flux_hi, flux_conc, flux_hs) IF(lzacc)

    !upwind estimate of tracer flux
    !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
    tmp(:,:,:) = p_ice%hi(:,:,:) * p_ice%conc(:,:,:)
    !$ACC END KERNELS
    CALL upwind_hflux_ice( p_patch_3D, tmp,  p_ice%vn_e, z_adv_flux_h, lacc=lzacc )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), patch_2D, p_op_coeff%div_coeff, flux_hi  (:,jk,:),&
        & 1, cells_in_domain, lacc=lzacc )
    ENDDO

    CALL upwind_hflux_ice( p_patch_3D, p_ice%conc, p_ice%vn_e, z_adv_flux_h, lacc=lzacc )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), patch_2D, p_op_coeff%div_coeff, flux_conc(:,jk,:),&
        & 1, cells_in_domain, lacc=lzacc )
    ENDDO

    !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
    tmp(:,:,:) = p_ice%hs(:,:,:) * p_ice%conc(:,:,:)
    !$ACC END KERNELS
    CALL upwind_hflux_ice( p_patch_3D, tmp, p_ice%vn_e, z_adv_flux_h, lacc=lzacc )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), patch_2D, p_op_coeff%div_coeff, flux_hs  (:,jk,:),&
        & 1, cells_in_domain, lacc=lzacc )
    ENDDO

    DO jk = 1,p_ice%kice
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

            ! will be divided by p_ice%conc(jc,jk,jb) after the -dtime*flux_conc is added
            p_ice%hi(jc,jk,jb)=p_ice%hi(jc,jk,jb)*p_ice%conc(jc,jk,jb)
            p_ice%hs(jc,jk,jb)=p_ice%hs(jc,jk,jb)*p_ice%conc(jc,jk,jb)

            p_ice%conc(jc,jk,jb)= p_ice%conc(jc,jk,jb)-dtime*flux_conc(jc,jk,jb)

            IF ( p_ice%conc(jc,jk,jb) > 0.0_wp) THEN
              p_ice%hi(jc,jk,jb)= (p_ice%hi(jc,jk,jb)-dtime*flux_hi(jc,jk,jb))/p_ice%conc(jc,jk,jb)
              p_ice%hs(jc,jk,jb)= (p_ice%hs(jc,jk,jb)-dtime*flux_hs(jc,jk,jb))/p_ice%conc(jc,jk,jb)
            ENDIF
          ENDIF
          ! TODO ram - remove patch_2D%cells%area(jc,jb) and test
          IF ( p_ice%conc(jc,jk,jb) > 0.0_wp) THEN
            p_ice%vol(jc,jk,jb) = p_ice%hi (jc,jk,jb)   &
              &         *( p_ice%conc(jc,jk,jb)*patch_2D%cells%area(jc,jb) )
            p_ice%vols(jc,jk,jb) = p_ice%hs(jc,jk,jb)   &
              &         *( p_ice%conc(jc,jk,jb)*patch_2D%cells%area(jc,jb) )
          ENDIF
        END DO
        !$ACC END PARALLEL LOOP
      END DO
    END DO

!--------------------------------------------------------------------------------------------------
! Sync results
!--------------------------------------------------------------------------------------------------
    CALL sync_patch_array_mult(sync_c, patch_2D, 5, &
      & p_ice%vol(:,:,:), p_ice%vols(:,:,:), p_ice%conc(:,:,:), p_ice%hs(:,:,:), p_ice%hi(:,:,:))

!     CALL sync_patch_array(SYNC_C, patch_2D, p_ice%vol (:,:,:))
!     CALL sync_patch_array(SYNC_C, patch_2D, p_ice%vols(:,:,:))
!     CALL sync_patch_array(SYNC_C, patch_2D, p_ice%conc(:,:,:))
!     CALL sync_patch_array(SYNC_C, patch_2D, p_ice%hs  (:,:,:))
!     CALL sync_patch_array(SYNC_C, patch_2D, p_ice%hi  (:,:,:))

!     CALL finish('sea_ice:upwind advect','Test sync ok')
   
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('ice_adv: vol ice'  , p_ice%vol , str_module, 4, in_subset=patch_2D%cells%owned)
    CALL dbg_print('ice_adv: vol snow' , p_ice%vols, str_module, 4, in_subset=patch_2D%cells%owned)
    CALL dbg_print('ice_adv: hi'       , p_ice%hi  , str_module, 4, in_subset=patch_2D%cells%owned)
    CALL dbg_print('ice_adv: hs'       , p_ice%hs  , str_module, 4, in_subset=patch_2D%cells%owned)
    CALL dbg_print('ice_adv: conc'     , p_ice%conc, str_module, 4, in_subset=patch_2D%cells%owned)
    !---------------------------------------------------------------------

    !$ACC END DATA

  END SUBROUTINE ice_advection_upwind

  !-------------------------------------------------------------------------
  !
  !> Advection of sea ice and snow on ice
  !! This uses the upwind_hflux_ice routine and the ocean's div_oce_3D routine to do upwind
  !! advection of the relevant variables.
  !!
  SUBROUTINE ice_advection_upwind_einar( p_patch_3D, p_op_coeff, p_ice )
    TYPE(t_patch_3D), TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_operator_coeff),   INTENT(IN)    :: p_op_coeff
    TYPE(t_sea_ice),          INTENT(INOUT) :: p_ice

    ! Local variables
    ! Patch and range
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_subset_range), POINTER :: cells_in_domain

    ! Indexing
    INTEGER  :: jk, jb, jc
    INTEGER  :: i_startidx_c, i_endidx_c

    ! Temporary variables/buffers
    REAL(wp) :: z_adv_flux_h (nproma,p_ice%kice,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: flux_hi  (nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: flux_conc(nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: flux_hs  (nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

!--------------------------------------------------------------------------------------------------
    patch_2D => p_patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
!--------------------------------------------------------------------------------------------------

    !upwind estimate of tracer flux
    CALL upwind_hflux_ice( p_patch_3D, p_ice%vol,  p_ice%vn_e, z_adv_flux_h )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), patch_2D, p_op_coeff%div_coeff, flux_hi  (:,jk,:),&
        & 1, cells_in_domain)
    ENDDO

    CALL upwind_hflux_ice( p_patch_3D, p_ice%conc, p_ice%vn_e, z_adv_flux_h )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), patch_2D, p_op_coeff%div_coeff, flux_conc(:,jk,:),&
        & 1, cells_in_domain)
    ENDDO

    CALL upwind_hflux_ice( p_patch_3D, p_ice%vols, p_ice%vn_e, z_adv_flux_h )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), patch_2D, p_op_coeff%div_coeff, flux_hs  (:,jk,:),&
        & 1, cells_in_domain)
    ENDDO

    DO jk = 1,p_ice%kice
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN
            p_ice%vol (jc,jk,jb)= p_ice%vol (jc,jk,jb)-dtime*flux_hi  (jc,jk,jb)
            p_ice%conc(jc,jk,jb)= p_ice%conc(jc,jk,jb)-dtime*flux_conc(jc,jk,jb)
            p_ice%vols(jc,jk,jb)= p_ice%vols(jc,jk,jb)-dtime*flux_hs  (jc,jk,jb)
          ENDIF
          ! TODO ram - remove patch_2D%cells%area(jc,jb) and test
          ! See also thermodyn/mo_sea_ice.f90
          IF ( p_ice%conc(jc,jk,jb) > 0.0_wp) THEN 
            p_ice%hi(jc,jk,jb) = p_ice%vol (jc,jk,jb)   &
              &         /( p_ice%conc(jc,jk,jb)*patch_2D%cells%area(jc,jb) )
            p_ice%hs(jc,jk,jb) = p_ice%vols(jc,jk,jb)   &
              &         /( p_ice%conc(jc,jk,jb)*patch_2D%cells%area(jc,jb) )
          ENDIF
        END DO
      END DO
    END DO

!--------------------------------------------------------------------------------------------------
! Sync results
!--------------------------------------------------------------------------------------------------

    CALL sync_patch_array(SYNC_C, patch_2D, p_ice%vol (:,:,:))
    CALL sync_patch_array(SYNC_C, patch_2D, p_ice%vols(:,:,:))
    CALL sync_patch_array(SYNC_C, patch_2D, p_ice%conc(:,:,:))
    CALL sync_patch_array(SYNC_C, patch_2D, p_ice%hs  (:,:,:))
    CALL sync_patch_array(SYNC_C, patch_2D, p_ice%hi  (:,:,:))

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('ice_adv: vol ice'  , p_ice%vol , str_module, 4, in_subset=patch_2D%cells%owned)
    CALL dbg_print('ice_adv: vol snow' , p_ice%vols, str_module, 4, in_subset=patch_2D%cells%owned)
    CALL dbg_print('ice_adv: hi'       , p_ice%hi  , str_module, 4, in_subset=patch_2D%cells%owned)
    CALL dbg_print('ice_adv: hs'       , p_ice%hs  , str_module, 4, in_subset=patch_2D%cells%owned)
    CALL dbg_print('ice_adv: conc'     , p_ice%conc, str_module, 4, in_subset=patch_2D%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE ice_advection_upwind_einar

  !-------------------------------------------------------------------------------
  !>
  !! First order upwind scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using the first
  !! order Godunov method.
  !!
  SUBROUTINE upwind_hflux_ice( p_patch_3D, pvar_c, pvn_e, pupflux_e, opt_slev, opt_elev, lacc )

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(IN)              :: pvar_c   (:,:,:) !< advected cell centered variable
    REAL(wp), INTENT(IN)              :: pvn_e    (:,:)   !< normal velocity on edges
    REAL(wp), INTENT(OUT)             :: pupflux_e(:,:,:) !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL     :: opt_slev    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL     :: opt_elev    ! optional vertical end level
    LOGICAL, INTENT(IN), OPTIONAL     :: lacc

    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER  :: slev, elev
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER         :: patch_2D
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    patch_2D         => p_patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-----------------------------------------------------------------------
    !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
    pupflux_e(:,:,:) = 0.0_wp
    !$ACC END KERNELS
    
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = UBOUND(pvar_c,2)
    END IF
    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
    ! for no-slip boundary conditions, boundary treatment for tracer (zero at leteral walls)
    !is implicit done via velocity boundary conditions
    !
    ! line and block indices of two neighboring cells
    iilc => patch_2D%edges%cell_idx
    iibc => patch_2D%edges%cell_blk

    ! loop through all patch edges (and blocks)
!ICONOMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,je,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block

      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)

!CDIR UNROLL=6
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) IF(lzacc)
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
!          IF ( p_patch_3D%lsm_e(je,1,jb) <= sea_boundary ) THEN
            pupflux_e(je,jk,jb) =  &
            &  laxfr_upflux( pvn_e(je,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
            &                                 pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
!          ELSE
!            pupflux_e(je,jk,jb) = 0.0_wp
!          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
      !$ACC END PARALLEL LOOP
    END DO  ! end loop over blocks
!ICON_OMP_END_DO_NOWAIT
!ICONOMP_END_PARALLEL

  END SUBROUTINE upwind_hflux_ice

END MODULE mo_ice_advection
