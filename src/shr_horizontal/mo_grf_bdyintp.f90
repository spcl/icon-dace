! Contains the interpolation routines needed for grid refinement.
!
! These had originally been included in mo_grf_interpolation but then were
! packed into a separate module to clean up the code
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

MODULE mo_grf_bdyintp
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                ONLY: wp
USE mo_model_domain,        ONLY: t_patch
USE mo_parallel_config,     ONLY: nproma, p_test_run, cpu_min_nproma
USE mo_communication,       ONLY: exchange_data_grf

USE mo_grf_intp_data_strc

IMPLICIT NONE

PRIVATE
PUBLIC :: interpol_vec_grf, interpol2_vec_grf, interpol_scal_grf


CONTAINS



!-------------------------------------------------------------------------
!
!
!>
!! Performs interpolation from parent edges to child edges.
!!
!! Performs interpolation from parent edges to child edges using inverse
!! distance weighting (IDW) or RBF. The output field (p_vn_out) contains
!! the normal wind components along the lateral boundaries of nested
!! model domains.
!!
SUBROUTINE interpol_vec_grf (p_pp, p_pc, p_grf, p_vn_in, p_vn_out)
  !
  TYPE(t_patch), TARGET, INTENT(in) :: p_pp
  TYPE(t_patch), TARGET, INTENT(inout) :: p_pc

  ! input normal components of vectors at edge midpoints
  REAL(wp), INTENT(IN) :: p_vn_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

  ! Indices of source points and interpolation coefficients
  TYPE(t_gridref_single_state),   TARGET, INTENT(IN) ::  p_grf

  ! reconstructed edge-normal wind component
  REAL(wp),INTENT(INOUT) :: p_vn_out(:,:,:) ! dim: (nproma,nlev_c,nblks_e)

  INTEGER :: jb, jk, je                ! loop indices
  INTEGER :: js                        ! shift parameter
  INTEGER :: nproma_bdyintp, nblks_bdyintp, npromz_bdyintp, nlen, nshift

  REAL(wp) :: vn_aux(p_pc%nlev,p_grf%npoints_bdyintp_e,4)

  ! Pointers to index fields/lists
  INTEGER,  DIMENSION(:,:),   POINTER :: iidx, iblk

  INTEGER :: nlev_c       !< number of vertical full levels (child domain)

  !-----------------------------------------------------------------------

  ! Set pointers to index lists
  iidx    => p_grf%idxlist_bdyintp_e
  iblk    => p_grf%blklist_bdyintp_e

  ! Compute values for dynamic nproma blocking
  nproma_bdyintp = cpu_min_nproma(nproma,256)
  nblks_bdyintp  = INT(p_grf%npoints_bdyintp_e/nproma_bdyintp)
  npromz_bdyintp = MOD(p_grf%npoints_bdyintp_e,nproma_bdyintp)
  IF (npromz_bdyintp > 0) THEN
    nblks_bdyintp = nblks_bdyintp + 1
  ELSE
    npromz_bdyintp = nproma_bdyintp
  ENDIF

  ! number of vertical full levels (child domain)
  nlev_c = p_pc%nlev
  ! difference between upper boundary of parent domain and
  ! upper boundary of child domain (in terms of vertical levels)
  js = p_pc%nshift

  !$ACC DATA CREATE(vn_aux) IF(i_am_accel_node)

!$OMP PARALLEL
!$OMP DO PRIVATE (jb,jk,je,nlen,nshift) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_bdyintp

    IF (jb == nblks_bdyintp) THEN
      nlen = npromz_bdyintp
    ELSE
      nlen = nproma_bdyintp
    ENDIF
    nshift = (jb-1)*nproma_bdyintp

    ! PRESENT(p_pp%edges%refin_ctrl) is required, as PGI 20.9 can not find this variable otherwise.
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(p_pp%edges%refin_ctrl) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO je = nshift+1, nshift+nlen
      DO jk = 1, nlev_c

        ! child edge 1
        vn_aux(jk,je,1) = p_grf%coeff_bdyintp_e12(1,je) * &
          p_vn_in(iidx(1,je),jk+js,iblk(1,je)) +          &
          p_grf%coeff_bdyintp_e12(2,je) *                 &
          p_vn_in(iidx(2,je),jk+js,iblk(2,je)) +          &
          p_grf%coeff_bdyintp_e12(3,je)*                  &
          p_vn_in(iidx(6,je),jk+js,iblk(6,je)) +          &
          p_grf%coeff_bdyintp_e12(4,je)*                  &
          p_vn_in(iidx(10,je),jk+js,iblk(10,je)) +        &
          p_grf%coeff_bdyintp_e12(5,je) *                 &
          p_vn_in(iidx(11,je),jk+js,iblk(11,je)) +        &
          p_grf%coeff_bdyintp_e12(6,je)*                  &
          p_vn_in(iidx(12,je),jk+js,iblk(12,je))

        ! child edge 2
        vn_aux(jk,je,2) = p_grf%coeff_bdyintp_e12(7,je) * &
          p_vn_in(iidx(1,je),jk+js,iblk(1,je)) +          &
          p_grf%coeff_bdyintp_e12(8,je) *                 &
          p_vn_in(iidx(3,je),jk+js,iblk(3,je)) +          &
          p_grf%coeff_bdyintp_e12(9,je)*                  &
          p_vn_in(iidx(7,je),jk+js,iblk(7,je)) +          &
          p_grf%coeff_bdyintp_e12(10,je)*                 &
          p_vn_in(iidx(13,je),jk+js,iblk(13,je)) +        &
          p_grf%coeff_bdyintp_e12(11,je) *                &
          p_vn_in(iidx(14,je),jk+js,iblk(14,je)) +        &
          p_grf%coeff_bdyintp_e12(12,je)*                 &
          p_vn_in(iidx(15,je),jk+js,iblk(15,je))

        ! child edge 3
        vn_aux(jk,je,3) = p_grf%coeff_bdyintp_e34(1,je) * &
          p_vn_in(iidx(1,je),jk+js,iblk(1,je)) +          &
          p_grf%coeff_bdyintp_e34(2,je) *                 &
          p_vn_in(iidx(2,je),jk+js,iblk(2,je)) +          &
          p_grf%coeff_bdyintp_e34(3,je)*                  &
          p_vn_in(iidx(3,je),jk+js,iblk(3,je)) +          &
          p_grf%coeff_bdyintp_e34(4,je)*                  &
          p_vn_in(iidx(4,je),jk+js,iblk(4,je)) +          &
          p_grf%coeff_bdyintp_e34(5,je) *                 &
          p_vn_in(iidx(5,je),jk+js,iblk(5,je))

        ! child edge 4
        IF (p_pp%edges%refin_ctrl(iidx(1,je),iblk(1,je)) == -1) CYCLE
        vn_aux(jk,je,4) = p_grf%coeff_bdyintp_e34(6,je) * &
          p_vn_in(iidx(1,je),jk+js,iblk(1,je)) +          &
          p_grf%coeff_bdyintp_e34(7,je) *                 &
          p_vn_in(iidx(6,je),jk+js,iblk(6,je)) +          &
          p_grf%coeff_bdyintp_e34(8,je)*                  &
          p_vn_in(iidx(7,je),jk+js,iblk(7,je)) +          &
          p_grf%coeff_bdyintp_e34(9,je)*                  &
          p_vn_in(iidx(8,je),jk+js,iblk(8,je)) +          &
          p_grf%coeff_bdyintp_e34(10,je) *                &
          p_vn_in(iidx(9,je),jk+js,iblk(9,je))

      ENDDO
    ENDDO
    !$ACC END PARALLEL

  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  ! Store results in p_vn_out
  CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_vec_grf,1,nlev_c, &
    RECV1=p_vn_out,SEND1=vn_aux)
  !$ACC WAIT(1)
  !$ACC END DATA ! vn_aux

END SUBROUTINE interpol_vec_grf

!
!>
!! Performs gradient-based interpolation from parent edges to child edges 1 and 2
!! and RBF/IDW-based interpolation to child edges 3 and 4
!!
SUBROUTINE interpol2_vec_grf (p_pp, p_pc, p_grf, nfields, f3din1, f3dout1, &
                              f3din2, f3dout2, f3din3, f3dout3, f3din4, f3dout4)


  TYPE(t_patch), TARGET, INTENT(in) :: p_pp
  TYPE(t_patch), TARGET, INTENT(inout) :: p_pc

  ! input normal components of vectors at edge midpoints; dim: (nproma,nlev,nblks_e)
  REAL(wp), INTENT(IN), DIMENSION(:,:,:), OPTIONAL, TARGET :: f3din1, f3din2, f3din3, f3din4

  ! Indices of source points and interpolation coefficients
  TYPE(t_gridref_single_state),   TARGET, INTENT(IN)    ::  p_grf

  ! number of fields provided on input (needed for aux fields and pointer allocation)
  INTEGER, INTENT(IN) :: nfields

  ! reconstructed edge-normal vector component; dim: (nproma,nlev_c,nblks_e)
  REAL(wp), INTENT(INOUT), DIMENSION(:,:,:), OPTIONAL :: f3dout1, f3dout2, f3dout3, f3dout4

  INTEGER :: jb, jk, je, jv, jn        ! loop indices
  INTEGER :: js                        ! shift parameter

  INTEGER :: nlevtot ! for MPI communication call
  INTEGER :: nproma_bdyintp, nblks_bdyintp_e, npromz_bdyintp_e, nblks_bdyintp_v, &
             npromz_bdyintp_v, nlen, nshift

  REAL(wp) :: dvn_tang

  ! Auxiliary fields
  REAL(wp) :: vn_aux(p_pc%nlev,p_grf%npoints_bdyintp_e,4,nfields)
  REAL(wp), DIMENSION(p_pc%nlev,p_grf%npoints_bdyintp_v,nfields) :: u_vert,v_vert

  ! Pointers to index fields/lists
  INTEGER,  DIMENSION(:,:),   POINTER :: iidx, iblk, ividx, ivblk, ievidx


  ! Allocatable pointer to input and output fields
  TYPE t_fieldptr
    REAL(wp), POINTER :: fld(:,:,:)
  END TYPE t_fieldptr
  TYPE(t_fieldptr) :: p_in(nfields)

  INTEGER :: nlev_c        !< number of vertical full levels (child domain)

  !-----------------------------------------------------------------------

  IF (p_test_run) vn_aux=0._wp

  IF (PRESENT(f3din1)) THEN
    p_in(1)%fld  => f3din1
  ENDIF
  IF (PRESENT(f3din2)) THEN
    p_in(2)%fld  => f3din2
  ENDIF
  IF (PRESENT(f3din3)) THEN
    p_in(3)%fld  => f3din3
  ENDIF
  IF (PRESENT(f3din4)) THEN
    p_in(4)%fld  => f3din4
  ENDIF

  ! Set pointers to index lists
  iidx    => p_grf%idxlist_bdyintp_e
  iblk    => p_grf%blklist_bdyintp_e
  ividx   => p_grf%idxlist_rbfintp_v
  ivblk   => p_grf%blklist_rbfintp_v
  ievidx  => p_grf%edge_vert_idx

  ! Compute values for dynamic nproma blocking
  nproma_bdyintp = cpu_min_nproma(nproma,256)
  nblks_bdyintp_e  = INT(p_grf%npoints_bdyintp_e/nproma_bdyintp)
  npromz_bdyintp_e = MOD(p_grf%npoints_bdyintp_e,nproma_bdyintp)
  IF (npromz_bdyintp_e > 0) THEN
    nblks_bdyintp_e = nblks_bdyintp_e + 1
  ELSE
    npromz_bdyintp_e = nproma_bdyintp
  ENDIF
  nblks_bdyintp_v  = INT(p_grf%npoints_bdyintp_v/nproma_bdyintp)
  npromz_bdyintp_v = MOD(p_grf%npoints_bdyintp_v,nproma_bdyintp)
  IF (npromz_bdyintp_v > 0) THEN
    nblks_bdyintp_v = nblks_bdyintp_v + 1
  ELSE
    npromz_bdyintp_v = nproma_bdyintp
  ENDIF

  ! number of vertical full levels (child domain)
  nlev_c = p_pc%nlev
  ! difference between upper boundary of parent domain and
  ! upper boundary of child domain (in terms of vertical levels
  js = p_pc%nshift

  !$ACC DATA CREATE(vn_aux, u_vert, v_vert) COPYIN(p_in) IF(i_am_accel_node)

  DO jn = 1, nfields
    ! Make sure that all fields in p_in are available from gpu via p_in(jn)%fld.
    !$ACC ENTER DATA ATTACH(p_in(jn)%fld) IF(i_am_accel_node)
  ENDDO

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jn,jk,jv,nlen,nshift) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_bdyintp_v

    IF (jb == nblks_bdyintp_v) THEN
      nlen = npromz_bdyintp_v
    ELSE
      nlen = nproma_bdyintp
    ENDIF
    nshift = (jb-1)*nproma_bdyintp

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jv = nshift+1, nshift+nlen
      DO jn = 1, nfields
        DO jk = 1, nlev_c

          u_vert(jk,jv,jn) =  &
            p_grf%coeff_rbf_v(1,1,jv)*p_in(jn)%fld(ividx(1,jv),jk+js,ivblk(1,jv)) + &
            p_grf%coeff_rbf_v(2,1,jv)*p_in(jn)%fld(ividx(2,jv),jk+js,ivblk(2,jv)) + &
            p_grf%coeff_rbf_v(3,1,jv)*p_in(jn)%fld(ividx(3,jv),jk+js,ivblk(3,jv)) + &
            p_grf%coeff_rbf_v(4,1,jv)*p_in(jn)%fld(ividx(4,jv),jk+js,ivblk(4,jv)) + &
            p_grf%coeff_rbf_v(5,1,jv)*p_in(jn)%fld(ividx(5,jv),jk+js,ivblk(5,jv)) + &
            p_grf%coeff_rbf_v(6,1,jv)*p_in(jn)%fld(ividx(6,jv),jk+js,ivblk(6,jv))
          v_vert(jk,jv,jn) =  &
            p_grf%coeff_rbf_v(1,2,jv)*p_in(jn)%fld(ividx(1,jv),jk+js,ivblk(1,jv)) + &
            p_grf%coeff_rbf_v(2,2,jv)*p_in(jn)%fld(ividx(2,jv),jk+js,ivblk(2,jv)) + &
            p_grf%coeff_rbf_v(3,2,jv)*p_in(jn)%fld(ividx(3,jv),jk+js,ivblk(3,jv)) + &
            p_grf%coeff_rbf_v(4,2,jv)*p_in(jn)%fld(ividx(4,jv),jk+js,ivblk(4,jv)) + &
            p_grf%coeff_rbf_v(5,2,jv)*p_in(jn)%fld(ividx(5,jv),jk+js,ivblk(5,jv)) + &
            p_grf%coeff_rbf_v(6,2,jv)*p_in(jn)%fld(ividx(6,jv),jk+js,ivblk(6,jv))

        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL

  ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jb,jn,jk,je,nlen,nshift,dvn_tang) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_bdyintp_e

    IF (jb == nblks_bdyintp_e) THEN
      nlen = npromz_bdyintp_e
    ELSE
      nlen = nproma_bdyintp
    ENDIF
    nshift = (jb-1)*nproma_bdyintp

    ! PRESENT(p_grf, iblk, iidx) is necessary to help PGI_20_9 with address-not-mapped errors
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(p_grf, iblk, iidx) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(3) PRIVATE(dvn_tang)
    DO je = nshift+1, nshift+nlen
      DO jn = 1, nfields
!DIR$ IVDEP
        DO jk = 1, nlev_c

          ! child edges 1 and 2
          dvn_tang = u_vert(jk,ievidx(2,je),jn) * p_grf%prim_norm(2,1,je) + &
                     v_vert(jk,ievidx(2,je),jn) * p_grf%prim_norm(2,2,je) - &
                    (u_vert(jk,ievidx(1,je),jn) * p_grf%prim_norm(1,1,je) + &
                     v_vert(jk,ievidx(1,je),jn) * p_grf%prim_norm(1,2,je) )

          vn_aux(jk,je,1,jn) = p_in(jn)%fld(iidx(1,je),jk+js,iblk(1,je)) + &
                               dvn_tang*p_grf%dist_pe2ce(1,je)
          vn_aux(jk,je,2,jn) = p_in(jn)%fld(iidx(1,je),jk+js,iblk(1,je)) + &
                               dvn_tang*p_grf%dist_pe2ce(2,je)

          ! child edge 3
          vn_aux(jk,je,3,jn) = p_grf%coeff_bdyintp_e34(1,je) * &
            p_in(jn)%fld(iidx(1,je),jk+js,iblk(1,je)) +        &
            p_grf%coeff_bdyintp_e34(2,je) *                    &
            p_in(jn)%fld(iidx(2,je),jk+js,iblk(2,je)) +        &
            p_grf%coeff_bdyintp_e34(3,je)*                     &
            p_in(jn)%fld(iidx(3,je),jk+js,iblk(3,je)) +        &
            p_grf%coeff_bdyintp_e34(4,je)*                     &
            p_in(jn)%fld(iidx(4,je),jk+js,iblk(4,je)) +        &
            p_grf%coeff_bdyintp_e34(5,je) *                    &
            p_in(jn)%fld(iidx(5,je),jk+js,iblk(5,je))

          ! child edge 4
          IF (p_pp%edges%refin_ctrl(iidx(1,je),iblk(1,je)) /= -1) THEN
            vn_aux(jk,je,4,jn) = p_grf%coeff_bdyintp_e34(6,je) * &
              p_in(jn)%fld(iidx(1,je),jk+js,iblk(1,je)) +        &
              p_grf%coeff_bdyintp_e34(7,je) *                    &
              p_in(jn)%fld(iidx(6,je),jk+js,iblk(6,je)) +        &
              p_grf%coeff_bdyintp_e34(8,je)*                     &
              p_in(jn)%fld(iidx(7,je),jk+js,iblk(7,je)) +        &
              p_grf%coeff_bdyintp_e34(9,je)*                     &
              p_in(jn)%fld(iidx(8,je),jk+js,iblk(8,je)) +        &
              p_grf%coeff_bdyintp_e34(10,je) *                   &
              p_in(jn)%fld(iidx(9,je),jk+js,iblk(9,je))
          ENDIF

        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL

  ENDDO ! blocks
  !$ACC WAIT(1)
!$OMP END DO
!$OMP END PARALLEL

  DO jn = 1, nfields
    !$ACC EXIT DATA DETACH(p_in(jn)%fld) IF(i_am_accel_node)
  ENDDO



  IF (nfields == 1) THEN
    nlevtot = nlev_c
    CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_vec_grf,nfields,nlevtot, &
      RECV1=f3dout1,SEND1=vn_aux(:,:,:,1) )

  ELSE IF (nfields == 2) THEN
    nlevtot = 2*nlev_c
    CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_vec_grf,nfields,nlevtot, &
      RECV1=f3dout1,SEND1=vn_aux(:,:,:,1),RECV2=f3dout2,  &
      SEND2=vn_aux(:,:,:,2) )

  ELSE IF (nfields == 3) THEN
    nlevtot = 3*nlev_c
    CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_vec_grf,nfields,nlevtot, &
      RECV1=f3dout1,SEND1=vn_aux(:,:,:,1),RECV2=f3dout2,  &
      SEND2=vn_aux(:,:,:,2),RECV3=f3dout3,SEND3=vn_aux(:,:,:,3) )

  ELSE IF (nfields == 4) THEN
    nlevtot = 4*nlev_c
    CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_vec_grf,nfields,nlevtot, &
      RECV1=f3dout1,SEND1=vn_aux(:,:,:,1),RECV2=f3dout2, &
      SEND2=vn_aux(:,:,:,2),RECV3=f3dout3,SEND3=vn_aux(:,:,:,3),  &
      RECV4=f3dout4,SEND4=vn_aux(:,:,:,4) )
  ENDIF

  !$ACC END DATA ! CREATE( vn_aux,u_vert,v_vert ) COPYIN(p_in)

END SUBROUTINE interpol2_vec_grf


!-------------------------------------------------------------------------
!
!>
!! Performs interpolation of scalar fields from parent cells to child
!! cells using the 2D gradient at the cell center.
!!
SUBROUTINE interpol_scal_grf (p_pp, p_pc, p_grf, nfields,&
                              f3din1, f3dout1, f3din2, f3dout2, f3din3, f3dout3, &
                              f3din4, f3dout4, f3din5, f3dout5, f3din6, f3dout6, &
                              f4din1, f4dout1, f4din2, f4dout2,                  &
                              llimit_nneg, lnoshift)
  !
  TYPE(t_patch), TARGET, INTENT(in) :: p_pp
  TYPE(t_patch), TARGET, INTENT(inout) :: p_pc


  ! Indices of source points and interpolation coefficients
  TYPE(t_gridref_single_state),   TARGET, INTENT(IN)    ::  p_grf

  ! number of fields provided on input (needed for aux fields and pointer allocation)
  INTEGER, INTENT(IN) :: nfields

  ! logical switch: if present and true, limit horizontal gradient so as to avoid negative values
  LOGICAL, INTENT(IN), OPTIONAL :: llimit_nneg(nfields)

  ! logical switch: if present and true, turn off accounting for shifts related to vertical nesting
  ! (needed for boundary interpolation of 2D diagnostic variables combined into a 3D field)
  LOGICAL, INTENT(IN), OPTIONAL :: lnoshift

  ! input scalar fields at cell points (up to six 3D fields or two 4D fields)
  REAL(wp), INTENT(IN), OPTIONAL, TARGET ::  & ! dim: (nproma,nlev,nblks_c)
    f3din1(:,:,:), f3din2(:,:,:), f3din3(:,:,:), f3din4(:,:,:), f3din5(:,:,:), f3din6(:,:,:), &
    f4din1(:,:,:,:), f4din2(:,:,:,:)

  ! reconstructed scalar output fields (up to six 3D fields or two 4D fields)
  REAL(wp), INTENT(INOUT), OPTIONAL, TARGET ::  & ! dim: (nproma,nlev_c,nblks_c)
    f3dout1(:,:,:), f3dout2(:,:,:), f3dout3(:,:,:), f3dout4(:,:,:), f3dout5(:,:,:), f3dout6(:,:,:), &
    f4dout1(:,:,:,:), f4dout2(:,:,:,:)

  INTEGER :: jb, jk, jc, jn, n         ! loop indices
  INTEGER :: js                        ! shift parameter
  INTEGER :: elev                      ! end index of vertical loop
  INTEGER :: nproma_bdyintp, nblks_bdyintp, npromz_bdyintp, nlen, nshift

  INTEGER :: nlevtot  ! for MPI communication call
  INTEGER :: n4d

  LOGICAL :: l_limit_nneg(nfields)     ! local variable corresponding to llimit_nneg
  LOGICAL :: l_noshift                 ! local variable corresponding to lnoshift
  LOGICAL :: l4d                       ! 4D field is provided as input

  ! Auxiliary fields
  REAL(wp), DIMENSION(MAX(90,p_pc%nlevp1),p_grf%npoints_bdyintp_c) :: &
    grad_x, grad_y, maxval_neighb, minval_neighb, val_ctr
  REAL(wp) :: h_aux(MAX(90,p_pc%nlevp1),p_grf%npoints_bdyintp_c,4,nfields)
  REAL(wp) :: limfac1, limfac2, limfac, min_expval, max_expval, epsi, ovsht_fac, r_ovsht_fac, &
              relaxed_minval, relaxed_maxval

  ! Pointers to index fields
  INTEGER, DIMENSION(:,:),   POINTER :: iidx, iblk

  ! Allocatable pointer to input and output fields
  TYPE t_fieldptr
    REAL(wp), POINTER :: fld(:,:,:)
  END TYPE t_fieldptr
  TYPE(t_fieldptr) :: p_in(nfields), p_out(nfields)

  !-----------------------------------------------------------------------

  IF (PRESENT(f4din1) .AND. .NOT. PRESENT(f4din2)) THEN
    DO n = 1, nfields
      p_in(n)%fld  => f4din1(:,:,:,n)
      p_out(n)%fld => f4dout1(:,:,:,n)
    ENDDO
    l4d = .TRUE.
  ELSE IF (PRESENT(f4din1) .AND. PRESENT(f4din2)) THEN
    n4d = nfields/2
    DO n = 1, n4d
      p_in(n)%fld  => f4din1(:,:,:,n)
      p_out(n)%fld => f4dout1(:,:,:,n)
    ENDDO
    DO n = 1, n4d
      p_in(n4d+n)%fld  => f4din2(:,:,:,n)
      p_out(n4d+n)%fld => f4dout2(:,:,:,n)
    ENDDO
    l4d = .TRUE.
  ELSE
    IF (PRESENT(f3din1)) THEN
      p_in(1)%fld  => f3din1
      p_out(1)%fld => f3dout1
    ENDIF
    IF (PRESENT(f3din2)) THEN
      p_in(2)%fld  => f3din2
      p_out(2)%fld => f3dout2
    ENDIF
    IF (PRESENT(f3din3)) THEN
      p_in(3)%fld  => f3din3
      p_out(3)%fld => f3dout3
    ENDIF
    IF (PRESENT(f3din4)) THEN
      p_in(4)%fld  => f3din4
      p_out(4)%fld => f3dout4
    ENDIF
    IF (PRESENT(f3din5)) THEN
      p_in(5)%fld  => f3din5
      p_out(5)%fld => f3dout5
    ENDIF
    IF (PRESENT(f3din6)) THEN
      p_in(6)%fld  => f3din6
      p_out(6)%fld => f3dout6
    ENDIF
    l4d = .FALSE.
  ENDIF

  ! Check if gradient limiting is required
  IF (PRESENT(llimit_nneg)) THEN
    l_limit_nneg(:) = llimit_nneg(:)
  ELSE
    l_limit_nneg(:) = .FALSE.
  ENDIF

  ! Check if vertical shifting is to be turned off
  IF (PRESENT(lnoshift)) THEN
    l_noshift = lnoshift
  ELSE
    l_noshift = .FALSE.
  ENDIF

  ! Compute values for dynamic nproma blocking
  nproma_bdyintp = cpu_min_nproma(nproma,256)
  nblks_bdyintp  = INT(p_grf%npoints_bdyintp_c/nproma_bdyintp)
  npromz_bdyintp = MOD(p_grf%npoints_bdyintp_c,nproma_bdyintp)
  IF (npromz_bdyintp > 0) THEN
    nblks_bdyintp = nblks_bdyintp + 1
  ELSE
    npromz_bdyintp = nproma_bdyintp
  ENDIF

  epsi = 1.e3_wp*TINY(1.0_wp) ! factor 1000 is needed to safely prevent overflow in divisions

  ovsht_fac = 1.05_wp ! factor of allowed overshooting
  r_ovsht_fac = 1._wp/ovsht_fac

  ! Pointers to index lists for gradient computation
  iidx => p_grf%idxlist_bdyintp_c
  iblk => p_grf%blklist_bdyintp_c


  IF (l_noshift) THEN
    js = 0
  ELSE
    js = p_pc%nshift
  ENDIF

  !$ACC DATA CREATE(grad_x, grad_y, maxval_neighb, minval_neighb, val_ctr, h_aux) &
  !$ACC   COPYIN(p_in) IF(i_am_accel_node)

  !$ACC SERIAL IF(i_am_accel_node)
  IF (p_test_run) h_aux(:,:,:,:) = 0._wp
  !$ACC END SERIAL

  DO jn = 1, nfields
    ! Make sure that all fields in p_in are available from gpu via p_in(jn)%fld.
    !$ACC ENTER DATA ATTACH(p_in(jn)%fld) IF(i_am_accel_node)
  ENDDO

! FIXME: PGI runs into deadlock on loop exit (?), if OMP-parallelized. Compiler bug suspected
!$OMP PARALLEL DO PRIVATE (jb,nlen,nshift,jk,jc,jn,elev,limfac1,limfac2,limfac, &
!$OMP   min_expval,max_expval,relaxed_minval,relaxed_maxval, grad_x, grad_y, &
!$OMP   val_ctr, maxval_neighb, minval_neighb) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_bdyintp

      nlen = MERGE(nproma_bdyintp, npromz_bdyintp, jb /= nblks_bdyintp)
      nshift = (jb-1)*nproma_bdyintp

      DO jn = 1, nfields
        elev   = UBOUND(p_out(jn)%fld,2)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jc = nshift+1, nshift+nlen
          DO jk = 1, elev
            val_ctr(jk,jc) = p_in(jn)%fld(iidx(1,jc),jk+js,iblk(1,jc))
            grad_x(jk,jc) =  &
              p_grf%coeff_bdyintp_c(1,1,jc)*p_in(jn)%fld(iidx(1,jc),jk+js,iblk(1,jc)) + &
              p_grf%coeff_bdyintp_c(2,1,jc)*p_in(jn)%fld(iidx(2,jc),jk+js,iblk(2,jc)) + &
              p_grf%coeff_bdyintp_c(3,1,jc)*p_in(jn)%fld(iidx(3,jc),jk+js,iblk(3,jc)) + &
              p_grf%coeff_bdyintp_c(4,1,jc)*p_in(jn)%fld(iidx(4,jc),jk+js,iblk(4,jc)) + &
              p_grf%coeff_bdyintp_c(5,1,jc)*p_in(jn)%fld(iidx(5,jc),jk+js,iblk(5,jc)) + &
              p_grf%coeff_bdyintp_c(6,1,jc)*p_in(jn)%fld(iidx(6,jc),jk+js,iblk(6,jc)) + &
              p_grf%coeff_bdyintp_c(7,1,jc)*p_in(jn)%fld(iidx(7,jc),jk+js,iblk(7,jc)) + &
              p_grf%coeff_bdyintp_c(8,1,jc)*p_in(jn)%fld(iidx(8,jc),jk+js,iblk(8,jc)) + &
              p_grf%coeff_bdyintp_c(9,1,jc)*p_in(jn)%fld(iidx(9,jc),jk+js,iblk(9,jc)) + &
              p_grf%coeff_bdyintp_c(10,1,jc)*p_in(jn)%fld(iidx(10,jc),jk+js,iblk(10,jc))
            grad_y(jk,jc) =  &
              p_grf%coeff_bdyintp_c(1,2,jc)*p_in(jn)%fld(iidx(1,jc),jk+js,iblk(1,jc)) + &
              p_grf%coeff_bdyintp_c(2,2,jc)*p_in(jn)%fld(iidx(2,jc),jk+js,iblk(2,jc)) + &
              p_grf%coeff_bdyintp_c(3,2,jc)*p_in(jn)%fld(iidx(3,jc),jk+js,iblk(3,jc)) + &
              p_grf%coeff_bdyintp_c(4,2,jc)*p_in(jn)%fld(iidx(4,jc),jk+js,iblk(4,jc)) + &
              p_grf%coeff_bdyintp_c(5,2,jc)*p_in(jn)%fld(iidx(5,jc),jk+js,iblk(5,jc)) + &
              p_grf%coeff_bdyintp_c(6,2,jc)*p_in(jn)%fld(iidx(6,jc),jk+js,iblk(6,jc)) + &
              p_grf%coeff_bdyintp_c(7,2,jc)*p_in(jn)%fld(iidx(7,jc),jk+js,iblk(7,jc)) + &
              p_grf%coeff_bdyintp_c(8,2,jc)*p_in(jn)%fld(iidx(8,jc),jk+js,iblk(8,jc)) + &
              p_grf%coeff_bdyintp_c(9,2,jc)*p_in(jn)%fld(iidx(9,jc),jk+js,iblk(9,jc)) + &
              p_grf%coeff_bdyintp_c(10,2,jc)*p_in(jn)%fld(iidx(10,jc),jk+js,iblk(10,jc))
            maxval_neighb(jk,jc) =                           &
              MAX(p_in(jn)%fld(iidx(1,jc),jk+js,iblk(1,jc)), &
                  p_in(jn)%fld(iidx(2,jc),jk+js,iblk(2,jc)), &
                  p_in(jn)%fld(iidx(3,jc),jk+js,iblk(3,jc)), &
                  p_in(jn)%fld(iidx(4,jc),jk+js,iblk(4,jc)), &
                  p_in(jn)%fld(iidx(5,jc),jk+js,iblk(5,jc)), &
                  p_in(jn)%fld(iidx(6,jc),jk+js,iblk(6,jc)), &
                  p_in(jn)%fld(iidx(7,jc),jk+js,iblk(7,jc)), &
                  p_in(jn)%fld(iidx(8,jc),jk+js,iblk(8,jc)), &
                  p_in(jn)%fld(iidx(9,jc),jk+js,iblk(9,jc)), &
                  p_in(jn)%fld(iidx(10,jc),jk+js,iblk(10,jc)))
            minval_neighb(jk,jc) =                               &
              MIN(p_in(jn)%fld(iidx(1,jc),jk+js,iblk(1,jc)), &
                  p_in(jn)%fld(iidx(2,jc),jk+js,iblk(2,jc)), &
                  p_in(jn)%fld(iidx(3,jc),jk+js,iblk(3,jc)), &
                  p_in(jn)%fld(iidx(4,jc),jk+js,iblk(4,jc)), &
                  p_in(jn)%fld(iidx(5,jc),jk+js,iblk(5,jc)), &
                  p_in(jn)%fld(iidx(6,jc),jk+js,iblk(6,jc)), &
                  p_in(jn)%fld(iidx(7,jc),jk+js,iblk(7,jc)), &
                  p_in(jn)%fld(iidx(8,jc),jk+js,iblk(8,jc)), &
                  p_in(jn)%fld(iidx(9,jc),jk+js,iblk(9,jc)), &
                  p_in(jn)%fld(iidx(10,jc),jk+js,iblk(10,jc)))
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(min_expval, max_expval) &
        !$ACC   PRIVATE(relaxed_minval, relaxed_maxval, limfac1, limfac2, limfac)
        DO jc = nshift+1, nshift+nlen
          DO jk = 1, elev

            min_expval = MIN(grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(1,1,jc) + &
                             grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(1,2,jc),  &
                             grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(2,1,jc) + &
                             grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(2,2,jc),  &
                             grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(3,1,jc) + &
                             grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(3,2,jc),  &
                             grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(4,1,jc) + &
                             grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(4,2,jc),  &
                             -epsi )

            max_expval = MAX(grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(1,1,jc) + &
                             grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(1,2,jc),  &
                             grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(2,1,jc) + &
                             grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(2,2,jc),  &
                             grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(3,1,jc) + &
                             grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(3,2,jc),  &
                             grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(4,1,jc) + &
                             grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(4,2,jc),  &
                             epsi )

            ! Allow a limited amount of over-/undershooting in the downscaled fields
            relaxed_minval = MERGE(r_ovsht_fac, ovsht_fac, &
              minval_neighb(jk,jc) > 0._wp) * minval_neighb(jk,jc)
            relaxed_maxval = MERGE(ovsht_fac, r_ovsht_fac, &
              maxval_neighb(jk,jc) > 0._wp) * maxval_neighb(jk,jc)

            limfac1 = MERGE(1._wp, ABS((relaxed_minval-val_ctr(jk,jc))/min_expval), &
              val_ctr(jk,jc) + min_expval >= relaxed_minval-epsi)
            limfac2 = MERGE(1._wp, ABS((relaxed_maxval-val_ctr(jk,jc))/max_expval), &
              val_ctr(jk,jc) + max_expval <= relaxed_maxval+epsi)
            limfac = MIN(limfac1,limfac2)

            grad_x(jk,jc) = grad_x(jk,jc)*limfac
            grad_y(jk,jc) = grad_y(jk,jc)*limfac

          ENDDO
        ENDDO
        !$ACC END PARALLEL

        IF (l_limit_nneg(jn)) THEN

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jc = nshift+1, nshift+nlen
            DO jk = 1, elev

              h_aux(jk,jc,1,jn) = MAX(0._wp, val_ctr(jk,jc) + &
                grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(1,1,jc)  + &
                grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(1,2,jc))
              h_aux(jk,jc,2,jn) = MAX(0._wp, val_ctr(jk,jc) + &
                grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(2,1,jc)  + &
                grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(2,2,jc))
              h_aux(jk,jc,3,jn) = MAX(0._wp, val_ctr(jk,jc) + &
                grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(3,1,jc)  + &
                grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(3,2,jc))
              h_aux(jk,jc,4,jn) = MAX(0._wp, val_ctr(jk,jc) + &
                grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(4,1,jc)  + &
                grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(4,2,jc))

            ENDDO
          ENDDO
          !$ACC END PARALLEL
        
        ELSE

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jc = nshift+1, nshift+nlen
            DO jk = 1, elev

              h_aux(jk,jc,1,jn) = val_ctr(jk,jc)           + &
                grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(1,1,jc) + &
                grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(1,2,jc)
              h_aux(jk,jc,2,jn) = val_ctr(jk,jc)           + &
                grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(2,1,jc) + &
                grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(2,2,jc)
              h_aux(jk,jc,3,jn) = val_ctr(jk,jc)           + &
                grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(3,1,jc) + &
                grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(3,2,jc)
              h_aux(jk,jc,4,jn) = val_ctr(jk,jc)           + &
                grad_x(jk,jc)*p_grf%dist_pc2cc_bdy(4,1,jc) + &
                grad_y(jk,jc)*p_grf%dist_pc2cc_bdy(4,2,jc)

            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

      ENDDO ! fields
    ENDDO ! blocks
    !$ACC WAIT(1)

! -------------------------------------

  DO jn = 1, nfields
    !$ACC EXIT DATA DETACH(p_in(jn)%fld) IF(i_am_accel_node)
  ENDDO


  ! Store results in p_out

  IF (l4d .AND. .NOT. PRESENT(f4din2)) THEN

    nlevtot = SIZE(f4dout1,2)*nfields
    CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_scal_grf,nfields,nlevtot, &
                           RECV4D1=f4dout1,SEND4D1=h_aux)

  ELSE IF (l4d) THEN

    nlevtot = SIZE(f4dout1,2)*n4d+SIZE(f4dout2,2)*n4d
    CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_scal_grf,nfields,nlevtot, &
                           RECV4D1=f4dout1,SEND4D1=h_aux(:,:,:,1:n4d),       &
                           RECV4D2=f4dout2,SEND4D2=h_aux(:,:,:,n4d+1:nfields)         )

  ELSE
    IF (nfields == 1) THEN
      nlevtot = SIZE(f3dout1,2)
      CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_scal_grf,nfields,nlevtot, &
                             RECV1=f3dout1,SEND1=h_aux(:,:,:,1)               )

    ELSE IF (nfields == 2) THEN
      nlevtot = SIZE(f3dout1,2)+SIZE(f3dout2,2)
      CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_scal_grf,nfields,nlevtot, &
                             RECV1=f3dout1,SEND1=h_aux(:,:,:,1),RECV2=f3dout2,&
                             SEND2=h_aux(:,:,:,2))

    ELSE IF (nfields == 3) THEN
      nlevtot = SIZE(f3dout1,2)+SIZE(f3dout2,2)+SIZE(f3dout3,2)
      CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_scal_grf,nfields,nlevtot, &
                             RECV1=f3dout1,SEND1=h_aux(:,:,:,1),RECV2=f3dout2,&
                             SEND2=h_aux(:,:,:,2),RECV3=f3dout3,SEND3=h_aux(:,:,:,3)   )

    ELSE IF (nfields == 4) THEN
      nlevtot = SIZE(f3dout1,2)+SIZE(f3dout2,2)+SIZE(f3dout3,2)+SIZE(f3dout4,2)
      CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_scal_grf,nfields,nlevtot, &
                             RECV1=f3dout1,SEND1=h_aux(:,:,:,1),RECV2=f3dout2,   &
                             SEND2=h_aux(:,:,:,2),RECV3=f3dout3,SEND3=h_aux(:,:,:,3),     &
                             RECV4=f3dout4,SEND4=h_aux(:,:,:,4)                           )

    ELSE IF (nfields == 5) THEN
      nlevtot = SIZE(f3dout1,2)+SIZE(f3dout2,2)+SIZE(f3dout3,2) + &
                SIZE(f3dout4,2)+SIZE(f3dout5,2)
      CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_scal_grf,nfields,nlevtot, &
                             RECV1=f3dout1,SEND1=h_aux(:,:,:,1),RECV2=f3dout2,&
                             SEND2=h_aux(:,:,:,2),RECV3=f3dout3,SEND3=h_aux(:,:,:,3),  &
                             RECV4=f3dout4,SEND4=h_aux(:,:,:,4),RECV5=f3dout5,         &
                             SEND5=h_aux(:,:,:,5)                                      )

    ELSE IF (nfields == 6) THEN
      nlevtot = SIZE(f3dout1,2)+SIZE(f3dout2,2)+SIZE(f3dout3,2) + &
                SIZE(f3dout4,2)+SIZE(f3dout5,2)+SIZE(f3dout6,2)
      CALL exchange_data_grf(p_pc%comm_pat_coll_interpol_scal_grf,nfields,nlevtot, &
                             RECV1=f3dout1,SEND1=h_aux(:,:,:,1),RECV2=f3dout2,&
                             SEND2=h_aux(:,:,:,2),RECV3=f3dout3,SEND3=h_aux(:,:,:,3),  &
                             RECV4=f3dout4,SEND4=h_aux(:,:,:,4),RECV5=f3dout5,         &
                             SEND5=h_aux(:,:,:,5),RECV6=f3dout6,SEND6=h_aux(:,:,:,6)   )

    ENDIF
  ENDIF
  !$ACC END DATA
END SUBROUTINE interpol_scal_grf

!-------------------------------------------------------------------------
END MODULE mo_grf_bdyintp
