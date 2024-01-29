!
! Subroutine interface_aes_car calls Cariolle's linearized ozone scheme.
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

MODULE mo_interface_aes_car

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: t_datetime => datetime

  USE mo_aes_phy_dims        ,ONLY: aes_phy_dims
  USE mo_aes_phy_config      ,ONLY: aes_phy_config, aes_phy_tc
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, prm_field, &
    &                               t_aes_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_car

  USE mo_run_config          ,ONLY: io3
  USE mo_physical_constants  ,ONLY: amd, amo3
  USE mo_bcs_time_interpolation ,ONLY: t_time_interpolation_weights, &
       &                               calculate_time_interpolation_weights
  USE mo_lcariolle     ,ONLY: t_avi, t_time_interpolation, lcariolle_do3dt 
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_aes_car

CONTAINS

  SUBROUTINE interface_aes_car(jg, jb, jcs, jce)

    ! Arguments
    !
    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    ! Pointers
    !
    TYPE(t_aes_phy_field), POINTER :: field
    TYPE(t_aes_phy_tend),  POINTER :: tend

    ! Local variables
    !
    INTEGER  :: nlev
    INTEGER  :: nproma
    !
    TYPE(t_datetime), POINTER :: datetime
    REAL(wp) :: pdtime
    LOGICAL  :: is_in_sd_ed_interval
    LOGICAL  :: is_active
    !
    INTEGER  :: fc_car
    !
    REAL(wp) :: tend_o3_car(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    !
    INTEGER                             :: jc, jk
    !
    TYPE(t_time_interpolation)          :: time_interpolation
    TYPE(t_time_interpolation_weights)  :: current_time_interpolation_weights
    TYPE(t_avi)                         :: avi

    IF (ltimer) call timer_start(timer_car)

    nlev    = aes_phy_dims(jg)%nlev
    nproma  = aes_phy_dims(jg)%nproma

    datetime             => aes_phy_tc(jg)%datetime
    pdtime               =  aes_phy_tc(jg)%dt_phy_sec
    is_in_sd_ed_interval =  aes_phy_tc(jg)%is_in_sd_ed_interval_car
    is_active            =  aes_phy_tc(jg)%is_active_car

    fc_car    =  aes_phy_config(jg)%fc_car

    ! associate pointers
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       !$ACC DATA PRESENT(field, field%qtrc_phy, tend) &
       !$ACC   CREATE(tend_o3_car)
       IF ( is_active ) THEN
          !
          current_time_interpolation_weights = calculate_time_interpolation_weights(datetime)
          time_interpolation% imonth1 = current_time_interpolation_weights% month1_index
          time_interpolation% imonth2 = current_time_interpolation_weights% month2_index
          time_interpolation% weight1 = current_time_interpolation_weights% weight1
          time_interpolation% weight2 = current_time_interpolation_weights% weight2
          !
          avi%ldown = .TRUE.
          avi%tmprt => field% ta  (:,:,jb)
          avi%pres  => field% pfull(:,:,jb)
          !
          ALLOCATE(avi%o3_vmr(nproma,nlev), avi%vmr2molm2(nproma,nlev), avi%cell_center_lat(nproma), avi%lday(nproma))
          !
          ! Note: ICON has no sources and sinks in te equation for air density. This implies
          !       that the total air mass is conserved. The parameterized turbulent mass flux
          !       at the surface and precipitation have no effect on the atmospheric mass.
          !       Therefore let us use here the total air mass as dry air mass.
          !
          !$ACC DATA CREATE(avi, avi%o3_vmr, avi%vmr2molm2, avi%cell_center_lat, avi%lday)
          !
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1, nlev
            DO jc = jcs, jce
              avi%o3_vmr(jc,jk)    = field% qtrc_phy(jc,jk,jb,io3)*amd/amo3
              avi%vmr2molm2(jc,jk) = field% mair(jc,jk,jb) / amd * 1.e3_wp
            END DO
          END DO
          !$ACC END PARALLEL LOOP

          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
          DO jc = jcs, jce
            avi%cell_center_lat(jc) =  field% clat(jc,jb)
            avi%lday(jc)            =  field% cosmu0(jc,jb) > 1.e-3_wp
          END DO
          !$ACC END PARALLEL LOOP
          !
          CALL lcariolle_do3dt(                              &
               & jcs, jce, nproma, nlev, time_interpolation, &
               & avi   = avi,                                &
               & do3dt = tend_o3_car(:,:)                    )
          !
          !$ACC WAIT
          !$ACC END DATA
          !
          DEALLOCATE(avi%o3_vmr, avi%vmr2molm2, avi%cell_center_lat, avi%lday)
          !
          ! store in memory for output or recycling
          IF (ASSOCIATED(tend% o3_car)) THEN
             !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
             DO jk = 1, nlev
               DO jc = jcs, jce
                  tend% o3_car(jc,jk,jb) = tend_o3_car(jc,jk)
               END DO
             END DO
             !$ACC END PARALLEL LOOP
          END IF
          !
       ELSE
          !
          ! retrieve from memory for recycling
          !
          IF (ASSOCIATED(tend% o3_car)) THEN
             !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
             DO jk = 1, nlev
               DO jc = jcs, jce
                 tend_o3_car(jc,jk) = tend% o3_car(jc,jk,jb)
               END DO
             END DO
             !$ACC END PARALLEL LOOP
          END IF
          !
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_car)
       CASE(1)
          !
          ! accumulate tendencies for later updating the model state
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1, nlev
            DO jc = jcs, jce
               tend% qtrc_phy(jc,jk,jb,io3) = tend% qtrc_phy(jc,jk,jb,io3) + tend_o3_car(jc,jk)*amo3/amd
            END DO
          END DO
          !$ACC END PARALLEL LOOP
          !
          ! update physics state for input to the next physics process
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1, nlev
             DO jc = jcs, jce
                field% qtrc_phy(jc,jk,jb,io3) = field% qtrc_phy(jc,jk,jb,io3) + tend_o3_car(jc,jk)*amo3/amd*pdtime
             END DO
          END DO
          !$ACC END PARALLEL LOOP
       END SELECT
       !$ACC WAIT
       !$ACC END DATA
       !
    ELSE
       !
       IF (ASSOCIATED(tend% o3_car)) THEN
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1, nlev
            DO jc = jcs, jce
              tend% o3_car(jc,jk,jb) = 0.0_wp
            END DO
          END DO
          !$ACC END PARALLEL LOOP
       END IF
       !
    END IF
       
    ! disassociate pointers
    NULLIFY(field)
    NULLIFY(tend )
    
    IF (ltimer) call timer_stop(timer_car)

  END SUBROUTINE interface_aes_car

END MODULE mo_interface_aes_car
