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

! For runs that check result consistency we fix the different
! contractions that the Intel compiler performs on some loops (at
! least in version 16.0) for the vectorized part and the
! non-vectorized parts
!----------------------------

MODULE mo_aes_diagnostics
  USE mo_kind                ,ONLY: wp
  USE mo_exception           ,ONLY: message
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: message, message_text
  USE mo_sync,                ONLY: global_max, global_min
  USE mo_run_config,          ONLY: iqv, iqc, iqi, iqr, iqs, iqg, iqh
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_statistics          ,ONLY: levels_horizontal_mean
  USE mo_name_list_output_init, ONLY: isRegistered
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, prm_field
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_impl_constants      ,ONLY: min_rlcell_int
  USE mo_impl_constants_grf  ,ONLY: grf_bdywidth_c

  PUBLIC aes_global_diagnostics, aes_diag_output_minmax_micro

CONTAINS
  SUBROUTINE aes_global_diagnostics(patch)
    TYPE(t_patch)  ,TARGET ,INTENT(in) :: patch
    REAL(wp)                           :: scr(nproma,patch%alloc_cell_blocks)

    REAL(wp) :: tas_gmean, rsdt_gmean, rsut_gmean, rlut_gmean, prec_gmean, evap_gmean, radtop_gmean, fwfoce_gmean
    TYPE(t_aes_phy_field), POINTER    :: field
    INTEGER  :: jb, jbs, jbe, jc, jcs, jce, rls, rle

    ! global mean t2m, tas_gmean, if requested for output
    tas_gmean = 0.0_wp
    IF ( isRegistered("tas_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%tas(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & tas_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%tas_gmean = tas_gmean

    ! global mean toa incident shortwave radiation, rsdt
    rsdt_gmean = 0.0_wp
    IF ( isRegistered("rsdt_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%rsdt(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & rsdt_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%rsdt_gmean = rsdt_gmean

    ! global mean toa outgoing shortwave radiation, rsut
    rsut_gmean = 0.0_wp
    IF ( isRegistered("rsut_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%rsut(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & rsut_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%rsut_gmean = rsut_gmean

    ! global mean toa outgoing longwave radiation, rlut
    rlut_gmean = 0.0_wp
    IF ( isRegistered("rlut_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%rlut(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & rlut_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%rlut_gmean = rlut_gmean

    ! global mean precipitation flux, prec
    prec_gmean = 0.0_wp
    IF ( isRegistered("prec_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%pr(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & prec_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%prec_gmean = prec_gmean

    ! global mean evaporation flux, evap
    evap_gmean = 0.0_wp
    IF ( isRegistered("evap_gmean") ) THEN
      call levels_horizontal_mean( prm_field(patch%id)%evap(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & evap_gmean, lopenacc=.TRUE.)
    END IF
    prm_field(patch%id)%evap_gmean = evap_gmean

    ! global mean toa total radiation, radtop, derived variable
    radtop_gmean = 0.0_wp
    IF ( isRegistered("radtop_gmean") ) THEN

      field => prm_field(patch%id)

      ! Compute row and block bounds for derived variables
      rls = grf_bdywidth_c + 1
      rle = min_rlcell_int
      jbs = patch%cells%start_blk(rls, 1)
      jbe = patch%cells%end_blk(rle, MAX(1, patch%n_childdom))

      !$ACC DATA PRESENT(field%rsdt, field%rsut, field%rlut) &
      !$ACC   CREATE(scr)

      DO jb = jbs, jbe
        CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = jcs, jce
          scr(jc,jb) = 0.0_wp
          scr(jc,jb) = field%rsdt(jc,jb) - field%rsut(jc,jb) - field%rlut(jc,jb)
        END DO
        !$ACC END PARALLEL
      END DO
      call levels_horizontal_mean( scr(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & radtop_gmean, lopenacc=.TRUE.)

      !$ACC WAIT(1)
      !$ACC END DATA

      NULLIFY(field)
    END IF
    prm_field(patch%id)%radtop_gmean = radtop_gmean

    ! global mean freshwater flux over ocean area, fwfoce, derived variable
    fwfoce_gmean = 0.0_wp
    IF ( isRegistered("fwfoce_gmean") ) THEN

      field => prm_field(patch%id)

      !$ACC DATA PRESENT(field%pr, field%evap, field%sftof) &
      !$ACC   CREATE(scr)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jb = 1, patch%alloc_cell_blocks
        DO jc = 1, nproma
          scr(jc,jb) = 0.0_wp
          scr(jc,jb) = (field%pr(jc,jb) + field%evap(jc,jb))*field%sftof(jc,jb)
        END DO
      END DO
      !$ACC END PARALLEL

      call levels_horizontal_mean( scr(:,:), &
          & patch%cells%area(:,:), &
          & patch%cells%owned, &
          & fwfoce_gmean, lopenacc=.TRUE.)

      !$ACC WAIT(1)
      !$ACC END DATA

      NULLIFY(field)
    END IF
    prm_field(patch%id)%fwfoce_gmean = fwfoce_gmean

    ! global mean ice cover fraction, icefrc - not set in atmosphere
 !  icefrc_gmean = 0.0_wp
 !  IF ( isRegistered("icefrc_gmean") ) THEN
 !    call levels_horizontal_mean( prm_field(patch%id)%icefrc(:,:), &
 !        & patch%cells%area(:,:), &
 !        & patch%cells%owned, &
 !        & icefrc_gmean)
 !  END IF
 !  prm_field(patch%id)%icefrc_gmean = icefrc_gmean

  END SUBROUTINE aes_global_diagnostics

  SUBROUTINE aes_diag_output_minmax_micro (patch,lpos)

    TYPE(t_patch)  ,TARGET ,INTENT(in) :: patch
    LOGICAL                ,INTENT(in) :: lpos

    ! Local variables
    REAL(wp), DIMENSION(patch%nblks_c) :: &
         & qvmax, qcmax, qrmax, qimax, qsmax, qhmax, qgmax, tmax, wmax, &
         & qvmin, qcmin, qrmin, qimin, qsmin, qhmin, qgmin, tmin, wmin
    REAL(wp) :: &
         & qvmaxi, qcmaxi, qrmaxi, qimaxi, qsmaxi, qhmaxi, qgmaxi, tmaxi, wmaxi, &
         & qvmini, qcmini, qrmini, qimini, qsmini, qhmini, qgmini, tmini, wmini

    ! loop indices
    INTEGER :: jc,jk,jb

    INTEGER :: nlev                    !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    IF(lpos) THEN
       ! called before microphysics
       CALL message('mo_aes_diagnostics:','input min/max values of microphysics')
    ELSE
       ! called after microphysics
       CALL message('mo_aes_diagnostics:','output min/max values of microphysics')
    ENDIF

    nlev = patch%nlev

    ! Exclude the nest boundary zone 

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_nchdom  = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)


    ! Find local min/max
    wmax  = 0.0_wp
    wmin  = 0.0_wp
    tmax  = -9999.0_wp
    tmin  =  9999.0_wp
    qvmax = 0.0_wp
    qvmin = 0.0_wp
    qcmax = 0.0_wp
    qcmin = 0.0_wp
    qrmax = 0.0_wp
    qrmin = 0.0_wp
    qimax = 0.0_wp
    qimin = 0.0_wp
    qsmax = 0.0_wp
    qsmin = 0.0_wp
    qgmax = 0.0_wp
    qgmin = 0.0_wp
    qhmax = 0.0_wp
    qhmin = 0.0_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
         DO jc = i_startidx, i_endidx
            wmax(jb)  = MAX(wmax(jb), prm_field(patch%id)%wa(jc,jk,jb))
            wmin(jb)  = MIN(wmin(jb), prm_field(patch%id)%wa(jc,jk,jb))
            tmax(jb)  = MAX(tmax(jb), prm_field(patch%id)%ta(jc,jk,jb))
            tmin(jb)  = MIN(tmin(jb), prm_field(patch%id)%ta(jc,jk,jb))
            qvmax(jb) = MAX(qvmax(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqv))
            qvmin(jb) = MIN(qvmin(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqv))
            qcmax(jb) = MAX(qcmax(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqc))
            qcmin(jb) = MIN(qcmin(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqc))
            qrmax(jb) = MAX(qrmax(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqr))
            qrmin(jb) = MIN(qrmin(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqr))
            qimax(jb) = MAX(qimax(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqi))
            qimin(jb) = MIN(qimin(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqi))
            qsmax(jb) = MAX(qsmax(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqs))
            qsmin(jb) = MIN(qsmin(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqs))
            qgmax(jb) = MAX(qgmax(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqg))
            qgmin(jb) = MIN(qgmin(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqg))
            qhmax(jb) = MAX(qhmax(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqh))
            qhmin(jb) = MIN(qhmin(jb),prm_field(patch%id)%qtrc_phy(jc,jk,jb,iqh))
         ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Take maximum/minimum over blocks
    wmaxi = MAXVAL(wmax(i_startblk:i_endblk))
    wmini = MINVAL(wmin(i_startblk:i_endblk))
    tmaxi  = MAXVAL(tmax(i_startblk:i_endblk))
    tmini  = MINVAL(tmin(i_startblk:i_endblk))
    qvmaxi = MAXVAL(qvmax(i_startblk:i_endblk))
    qvmini = MINVAL(qvmin(i_startblk:i_endblk))
    qcmaxi = MAXVAL(qcmax(i_startblk:i_endblk))
    qcmini = MINVAL(qcmin(i_startblk:i_endblk))
    qrmaxi = MAXVAL(qrmax(i_startblk:i_endblk))
    qrmini = MINVAL(qrmin(i_startblk:i_endblk))
    qimaxi = MAXVAL(qimax(i_startblk:i_endblk))
    qimini = MINVAL(qimin(i_startblk:i_endblk))
    qsmaxi = MAXVAL(qsmax(i_startblk:i_endblk))
    qsmini = MINVAL(qsmin(i_startblk:i_endblk))
    qgmaxi = MAXVAL(qgmax(i_startblk:i_endblk))
    qgmini = MINVAL(qgmin(i_startblk:i_endblk))
    qhmaxi = MAXVAL(qhmax(i_startblk:i_endblk))
    qhmini = MINVAL(qhmin(i_startblk:i_endblk))

    ! Take maximum/minimum over all PEs
    wmaxi  = global_max(wmaxi)
    wmini  = global_min(wmini)
    tmaxi  = global_max(tmaxi)
    tmini  = global_min(tmini)
    qvmaxi = global_max(qvmaxi)
    qvmini = global_min(qvmini)
    qcmaxi = global_max(qcmaxi)
    qcmini = global_min(qcmini)
    qrmaxi = global_max(qrmaxi)
    qrmini = global_min(qrmini)
    qimaxi = global_max(qimaxi)
    qimini = global_min(qimini)
    qsmaxi = global_max(qsmaxi)
    qsmini = global_min(qsmini)
    qgmaxi = global_max(qgmaxi)
    qgmini = global_min(qgmini)
    qhmaxi = global_max(qhmaxi)
    qhmini = global_min(qhmini)

    ! Standard output

     WRITE(message_text,'(A10,9A11)')   '  var: ', 'w','qv','qc','qr','qi','qs','qg','qh','temp'
       CALL message("",TRIM(message_text))
     WRITE(message_text,'(A10,9E11.3)') '  max: ', wmaxi,qvmaxi,qcmaxi,qrmaxi,qimaxi,qsmaxi,qgmaxi,qhmaxi,tmaxi
       CALL message("",TRIM(message_text))
     WRITE(message_text,'(A10,9E11.3)') '  min: ', wmini,qvmini,qcmini,qrmini,qimini,qsmini,qgmini,qhmini,tmini
       CALL message("",TRIM(message_text))

  END SUBROUTINE aes_diag_output_minmax_micro

END MODULE mo_aes_diagnostics
