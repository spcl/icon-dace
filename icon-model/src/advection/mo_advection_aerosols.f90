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

! Routines for idealized transport of 2D aerosol fields

MODULE mo_advection_aerosols

  USE mo_kind,                ONLY: wp, vp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_impl_constants,      ONLY: min_rlcell_int, min_rledge_int, nclass_aero, idu
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_advection_config,    ONLY: advection_config
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_edge
  USE mo_math_gradients,      ONLY: grad_green_gauss_cell
  USE mo_math_divrot,         ONLY: recon_lsq_cell_l_svd
  USE mo_advection_traj,      ONLY: t_back_traj, btraj_compute_o1
  USE mo_exception,           ONLY: message, message_text

  IMPLICIT NONE

  PRIVATE



  PUBLIC :: aerosol_2D_advection, setup_aerosol_advection


CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! Computes indices for vertically averaged fluxes used for transport of 2D aerosol fields
  SUBROUTINE setup_aerosol_advection(p_patch)

    TYPE(t_patch),     TARGET, INTENT(IN) :: p_patch  !< patch of current domain

    INTEGER :: nlev, jk, jk1, jg, i

    REAL(wp) :: zf, zstart(2), zend(2)

    nlev = p_patch%nlev
    jg   = p_patch%id

    ! Start and end heights for vertical averaging of mass fluxes
    ! 1. all aerosol classes except dust
    zstart(1) =  500._wp
    zend(1)   = 1500._wp
    ! 2. dust (has a larger decay scale than the others)
    zstart(2) = 2000._wp
    zend(2)   = 6000._wp

    ! Compute corresponding level indices
!$NEC novector
    DO jk = nlev, 1, -1
      jk1 = jk + p_patch%nshift_total
      zf = 0.5_wp*(vct_a(jk1)+vct_a(jk1+1))
      DO i = 1, 2
        IF (zf >= zstart(i) .AND. advection_config(jg)%kend_aero(i) == 0) THEN
          advection_config(jg)%kend_aero(i) = jk
        ENDIF
        IF (zf >= zend(i) .AND. advection_config(jg)%kstart_aero(i) == 0) THEN
          advection_config(jg)%kstart_aero(i) = jk
        ENDIF
      ENDDO
    ENDDO

    WRITE(message_text,'(a,i2)') '2D aerosol advection: start and end levels for flux averaging, domain ',jg
    CALL message('',message_text)
    WRITE(message_text,'(2(a,2i5))') 'kstart: ',advection_config(jg)%kstart_aero(1:2), &
                                     ', kend: ',advection_config(jg)%kend_aero(1:2)
    CALL message('',message_text)

  END SUBROUTINE setup_aerosol_advection

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Driver routine for idealized transport of 2D aerosol fields
  !!
  SUBROUTINE aerosol_2D_advection(p_patch, p_int, iprog_aero, dtime, aerosol, vn_traj, mflx_h, mflx_v, &
                                  deltaz_e, rhodz_now, rhodz_new)


    TYPE(t_patch),     TARGET, INTENT(IN) :: p_patch  ! patch of current domain

    TYPE(t_int_state), TARGET, INTENT(IN) :: p_int    ! interpolation state

    INTEGER,  INTENT(IN) :: iprog_aero      ! option for prognostic treatment (1 = dust only, 2 = all)
    REAL(wp), INTENT(IN) :: dtime           ! advection time step

    REAL(wp), INTENT(INOUT) :: aerosol(:,:,:)   ! 2D aerosol optical depth fields (middle index = aerosol class)

    REAL(wp), INTENT(IN)    :: vn_traj(:,:,:)   ! edge-normal velocity for back-trajectory calculation
    REAL(wp), INTENT(IN)    :: mflx_h(:,:,:)    ! horizontal mass flux at edges (vn*rho*deltaz)
    REAL(wp), INTENT(IN)    :: mflx_v(:,:,:)    ! vertical mass flux at cell interfaces (w*rho)
    REAL(vp), INTENT(IN)    :: deltaz_e(:,:,:)  ! layer thickness at edges
    REAL(wp), INTENT(IN)    :: rhodz_now(:,:,:) ! air mass per unit area in cell at time level now
    REAL(wp), INTENT(IN)    :: rhodz_new(:,:,:) ! air mass per unit area in cell at time level new


    ! Vertically averaged edge-normal and tangential velocities
    REAL(wp) :: vn_traj_avg(nproma,2,p_patch%nblks_e), vt_traj_avg(nproma,2,p_patch%nblks_e)

    REAL(wp) :: mflx_h_int(nproma,2,p_patch%nblks_e)  !< vertically integrated horizontal mass flux

    ! Vertically integrated air masses per unit area at both time levels
    REAL(wp) :: rhodz_now_int(nproma,nclass_aero,p_patch%nblks_c), rhodz_new_int(nproma,nclass_aero,p_patch%nblks_c)

    ! backward trajectory information
    TYPE(t_back_traj) :: btraj

    ! Horizontal gradient field of aerosols
    REAL(wp) :: lsq_aero (3,nproma,nclass_aero,p_patch%nblks_c)
    REAL(vp) :: grad_aero(2,nproma,nclass_aero,p_patch%nblks_c)
    LOGICAL  :: use_zlsq                       !< true if lsq_aero is used to store the gradients
    !
    REAL(wp) :: flx_aero(nproma,nclass_aero,p_patch%nblks_e)

    REAL(vp) ::  fluxdiv_c(nproma,nclass_aero), dz(nproma)
    REAL(wp) ::  dthalf

    INTEGER  :: jb, jk, jt, jc, je, jg, ilc, ibc, kst, kend
    INTEGER  :: i_startblk, i_startidx, i_endblk, i_endidx
    INTEGER  :: i_rlstart, i_rlend, jtstart, jtend, jtstep

    ! Mapping array between 5 aerosol types and 2 layer thickness classes (iss,iorg,ibc,iso4,idu)
    INTEGER, PARAMETER :: ji(nclass_aero) = (/1,1,1,1,2/)

    ! Pointer to index fields
    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

   !-----------------------------------------------------------------------

    CALL btraj%construct(nproma,p_patch%nlev,p_patch%nblks_e,2)

    jg  = p_patch%id

    dthalf = 0.5_wp*dtime

    ! line and block indices of edges as seen from cells
    iidx => p_patch%cells%edge_idx
    iblk => p_patch%cells%edge_blk

    ! Compute vertically averaged mass fluxes, back-trajectory velocities, and air masses
    IF (iprog_aero == 1) THEN
      jtstart = 2 ! dust only
    ELSE
      jtstart = 1 ! all aerosol classes
    ENDIF
    jtend = 2

!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

    i_rlstart  = 5
    i_rlend    = min_rledge_int-2
    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

!$OMP DO PRIVATE(jb,jk,jt,je,i_startidx,i_endidx,dz)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,       &
                          i_startidx, i_endidx, i_rlstart, i_rlend )


      vn_traj_avg(:,:,jb) = 0._wp
      mflx_h_int(:,:,jb)  = 0._wp

      DO jt = jtstart, jtend
        dz(:) = 0._vp

        DO jk = advection_config(jg)%kstart_aero(jt), advection_config(jg)%kend_aero(jt)
          DO je = i_startidx, i_endidx

            vn_traj_avg(je,jt,jb) = vn_traj_avg(je,jt,jb) + deltaz_e(je,jk,jb)*vn_traj(je,jk,jb)
            mflx_h_int(je,jt,jb)  = mflx_h_int(je,jt,jb)  + mflx_h(je,jk,jb)
            dz(je) = dz(je) + deltaz_e(je,jk,jb)

          ENDDO
        ENDDO
        DO je = i_startidx, i_endidx
          vn_traj_avg(je,jt,jb) = vn_traj_avg(je,jt,jb) / dz(je)
          ! ensure that averaged velocity and averaged mass flux have the same sign
          IF (vn_traj_avg(je,jt,jb)*mflx_h_int(je,jt,jb) < 0._wp) THEN
            vn_traj_avg(je,jt,jb) = SIGN(0.01_wp,mflx_h_int(je,jt,jb))
          ENDIF
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO

    i_rlstart  = grf_bdywidth_c-1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    IF (iprog_aero == 1) THEN
      jtstart = idu ! dust only
      jtend   = idu
      jtstep  = 1
    ELSE
      jtstart = 1 ! all aerosol classes
      jtend   = nclass_aero
      jtstep  = 4
    ENDIF

!$OMP DO PRIVATE(jb,jk,jt,jc,i_startidx,i_endidx,kst,kend)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
                          i_startidx, i_endidx, i_rlstart, i_rlend )

      ! Decay scales for all aerosol types but dust are the same. Therefore, the computation is done
      ! only for types 1 and 5, and the remaining ones are filled afterwards
      DO jt = jtstart, jtend, jtstep
        rhodz_now_int(:,jt,jb) = 0._wp
        rhodz_new_int(:,jt,jb) = 0._wp

        kst  = advection_config(jg)%kstart_aero(ji(jt))
        kend = advection_config(jg)%kend_aero(ji(jt))

        DO jk = kst, kend
          DO jc = i_startidx, i_endidx

            rhodz_now_int(jc,jt,jb) = rhodz_now_int(jc,jt,jb) + rhodz_now(jc,jk,jb)
            rhodz_new_int(jc,jt,jb) = rhodz_new_int(jc,jt,jb) + rhodz_new(jc,jk,jb)

          ENDDO
        ENDDO

        ! Account for vertical fluxes across top and bottom levels
        DO jc = i_startidx, i_endidx
          rhodz_now_int(jc,jt,jb) = rhodz_now_int(jc,jt,jb) + dtime * &
            (mflx_v(jc,kend+1,jb) - mflx_v(jc,kst,jb))
        ENDDO

      ENDDO

      ! Copy integrated values from type 1 to types 2-4
      IF (iprog_aero > 1) THEN
        DO jc = i_startidx, i_endidx
          rhodz_now_int(jc,2:4,jb) = rhodz_now_int(jc,1,jb)
          rhodz_new_int(jc,2:4,jb) = rhodz_new_int(jc,1,jb)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Compute tangential component of vertically averaged back-trajectory winds
    CALL rbf_vec_interpol_edge(vn_traj_avg, p_patch, p_int, vt_traj_avg,               &
                               opt_rlstart=grf_bdywidth_e-2, opt_rlend=min_rledge_int-1)

    ! Compute back trajectories
    CALL btraj_compute_o1 ( btraj    = btraj,            & !inout
      &                  ptr_p       = p_patch,          & !in
      &                  ptr_int     = p_int,            & !in
      &                  p_vn        = vn_traj_avg,      & !in
      &                  p_vt        = vt_traj_avg,      & !in
      &                  p_dthalf    = dthalf,           & !in
      &                  opt_rlstart = grf_bdywidth_e-2, & !in
      &                  opt_rlend   = min_rledge_int-1, & !in
      &                  opt_elev    = 2                 ) !in


    ! Reconstruct 2D gradient fields of aerosol
    IF (advection_config(jg)%igrad_c_miura == 1 .AND. advection_config(jg)%llsq_svd) THEN
      use_zlsq = .TRUE.
      CALL recon_lsq_cell_l_svd(aerosol, p_patch, p_int%lsq_lin, lsq_aero, opt_rlend=min_rlcell_int-1, &
                                opt_slev = jtstart, opt_elev = jtend)
    ELSE
      use_zlsq = .FALSE.
      CALL grad_green_gauss_cell(aerosol, p_patch, p_int, grad_aero, opt_rlend=min_rlcell_int-1, &
                                opt_slev = jtstart, opt_elev = jtend)
    ENDIF


!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

    i_rlstart  = grf_bdywidth_e-2
    i_rlend    = min_rledge_int-1
    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

!$OMP DO PRIVATE(jb,jt,je,i_startidx,i_endidx,ilc,ibc)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      IF ( use_zlsq ) THEN 
        DO jt = jtstart, jtend
          DO je = i_startidx, i_endidx

            ilc = btraj%cell_idx(je,ji(jt),jb)
            ibc = btraj%cell_blk(je,ji(jt),jb)
            flx_aero(je,jt,jb) = ( lsq_aero(1,ilc,jt,ibc)                                  &
              &                + btraj%distv_bary(je,ji(jt),jb,1)*lsq_aero(2,ilc,jt,ibc)   &
              &                + btraj%distv_bary(je,ji(jt),jb,2)*lsq_aero(3,ilc,jt,ibc) ) &
              &                * mflx_h_int(je,ji(jt),jb)

          ENDDO
        ENDDO
      ELSE
        DO jt = jtstart, jtend
          DO je = i_startidx, i_endidx

            ilc = btraj%cell_idx(je,ji(jt),jb)
            ibc = btraj%cell_blk(je,ji(jt),jb)
            flx_aero(je,jt,jb) = ( aerosol(ilc,jt,ibc)                                      &
              &                + btraj%distv_bary(je,ji(jt),jb,1)*grad_aero(1,ilc,jt,ibc)   &
              &                + btraj%distv_bary(je,ji(jt),jb,2)*grad_aero(2,ilc,jt,ibc) ) &
              &                * mflx_h_int(je,ji(jt),jb)

          ENDDO
        ENDDO
      ENDIF  ! use_zlsq

    ENDDO
!$OMP END DO
!$OMP END PARALLEL



!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

    i_rlstart  = grf_bdywidth_c+1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP DO PRIVATE(jb,jt,jc,i_startidx,i_endidx,fluxdiv_c)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
                          i_startidx, i_endidx, i_rlstart, i_rlend )


      DO jc = i_startidx, i_endidx
        DO jt = jtstart, jtend

          fluxdiv_c(jc,jt) =                                                     &
            flx_aero(iidx(jc,jb,1),jt,iblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
            flx_aero(iidx(jc,jb,2),jt,iblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
            flx_aero(iidx(jc,jb,3),jt,iblk(jc,jb,3))*p_int%geofac_div(jc,3,jb)

        ENDDO
      ENDDO

      DO jt = jtstart, jtend
        DO jc = i_startidx, i_endidx
          aerosol(jc,jt,jb) = MAX(0._wp, ( aerosol(jc,jt,jb)*rhodz_now_int(jc,jt,jb) - &
            dtime*fluxdiv_c(jc,jt) ) / rhodz_new_int(jc,jt,jb))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    CALL btraj%destruct()

  END SUBROUTINE aerosol_2D_advection


END MODULE mo_advection_aerosols

