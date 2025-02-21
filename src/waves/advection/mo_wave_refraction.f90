! Refraction of spectral surface wave energy
!
! Literature: The WAM model - A third Generation Ocean Wave Prediction Model, 1988
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
!
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
MODULE mo_wave_refraction

  USE mo_kind,                ONLY: wp, vp
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell
  USE mo_model_domain,        ONLY: t_patch
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_grid_config,         ONLY: grid_sphere_radius
  USE mo_parallel_config,     ONLY: nproma
  USE mo_math_constants,      ONLY: pi, pi2
  USE mo_loopindices,         ONLY: get_indices_c

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_refraction'

  PUBLIC :: wave_refraction

CONTAINS

  !>
  !! Calculate grid, depth and current refraction
  !!
  !! Based on WAM shallow water with depth and current refraction
  !! P_SPHER_SHALLOW_CURR
  !!
  SUBROUTINE wave_refraction(p_patch, wave_config, dtime, wave_num_c, gv_c, depth, &
    &                        depth_grad, tracer_now, tracer_new)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         routine = modname//':wave_refraction'

    TYPE(t_patch),               INTENT(IN)   :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)   :: wave_config
    REAL(wp),                    INTENT(IN)   :: dtime               ! integration time step [s]
    REAL(wp),                    INTENT(IN)   :: wave_num_c(:,:,:)
    REAL(wp),                    INTENT(IN)   :: gv_c(:,:,:)         ! group velocity at cell centers
    REAL(wp),                    INTENT(IN)   :: depth(:,:)
    REAL(vp),                    INTENT(IN)   :: depth_grad(:,:,:,:) ! bathymetry gradient (2,jc,jk,jb)
    REAL(wp),                    INTENT(IN)   :: tracer_now(:,:,:,:) ! energy before transport
    REAL(wp),                    INTENT(INOUT):: tracer_new(:,:,:,:)


    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd,jk
    INTEGER :: jt,jtm1,jtp1                       !< tracer index

    REAL(wp) :: DELTHR, DELTH, DELTR, DELTH0, sm, sp, ak, akd, DTP, DTM, dDTC, temp, tsihkd
    REAL(wp) :: thdd(nproma,wave_config%ndirs)
    REAL(wp) :: delta_ref(nproma,wave_config%nfreqs*wave_config%ndirs)

    wc => wave_config

    DELTH = 2.0_wp*pi/REAL(wc%ndirs,wp)
    DELTR = DELTH * grid_sphere_radius
    DELTH0 = 0.5_wp * dtime / DELTR
    DELTHR = 0.5_wp * dtime / DELTH

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jc,jt,i_startidx,i_endidx,temp,ak,akd,tsihkd,thdd, &
!$OMP            sm,sp,jtm1,jtp1,dtp,dtm,dDTC,delta_ref) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx

            temp = (SIN(wc%dirs(jd)) + SIN(wc%dirs(wc%dir_neig_ind(2,jd)))) * depth_grad(2,jc,jk,jb) &
                 - (COS(wc%dirs(jd)) + COS(wc%dirs(wc%dir_neig_ind(2,jd)))) * depth_grad(1,jc,jk,jb)

            ak = wave_num_c(jc,jb,jf)
            akd = ak * depth(jc,jb)

            IF (akd <= 10.0_wp) THEN
              tsihkd = (pi2 * wc%freqs(jf))/SINH(2.0_wp*akd)
            ELSE
              tsihkd = 0.0_wp
            END IF

            thdd(jc,jd) = temp * tsihkd

          END DO !jc
        END DO !jd


        DO jd = 1,wc%ndirs

          sm = DELTH0 * (SIN(wc%dirs(jd)) + SIN(wc%dirs(wc%dir_neig_ind(1,jd)))) !index of direction - 1
          sp = DELTH0 * (SIN(wc%dirs(jd)) + SIN(wc%dirs(wc%dir_neig_ind(2,jd)))) !index of direction + 1

          jt   = wc%tracer_ind(jd,jf)
          jtm1 = wc%tracer_ind(wc%dir_neig_ind(1,jd),jf)
          jtp1 = wc%tracer_ind(wc%dir_neig_ind(2,jd),jf)

          DO jc = i_startidx, i_endidx

            DTP = SIN(p_patch%cells%center(jc,jb)%lat) / COS(p_patch%cells%center(jc,jb)%lat) * gv_c(jc,jb,jf)

            DTM = DTP * SM + thdd(jc,wc%dir_neig_ind(1,jd)) * DELTHR
            DTP = DTP * SP + thdd(jc,jd) * DELTHR

            dDTC = -MAX(0._wp , DTP) + MIN(0._wp , DTM)
            DTP  = -MIN(0._wp , DTP)
            DTM  =  MAX(0._wp , DTM)

            delta_ref(jc,jt) = dDTC * tracer_now(jc,jk,jb,jt) &
                 + DTM * tracer_now(jc,jk,jb,jtm1)  &
                 + DTP * tracer_now(jc,jk,jb,jtp1)
          END DO !jc
        END DO !jd
      END DO !jf

      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          jt = wc%tracer_ind(jd,jf)
          DO jc = i_startidx, i_endidx
            tracer_new(jc,jk,jb,jt) = tracer_new(jc,jk,jb,jt) + delta_ref(jc,jt)
          END DO !jc
        END DO !jd
      END DO !jf

    END DO !jb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE wave_refraction


END MODULE mo_wave_refraction
