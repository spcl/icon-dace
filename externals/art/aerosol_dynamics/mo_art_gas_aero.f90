!
! mo_art_gas_aero
! Interaction aerosol gas phase
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

MODULE mo_art_gas_aero
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_nonhydro_types,                ONLY: t_nh_prog, t_nh_diag
  USE mo_loopindices,                   ONLY: get_indices_c
  USE mo_run_config,                    ONLY: iqv
  USE mo_util_phys,                     ONLY: rel_hum
  USE mo_impl_constants,                ONLY: SUCCESS
  
! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom
  USE art_isorropia,                    ONLY: isoropia
  
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_gas_aero'

  PUBLIC :: art_prepare_isorropia
  PUBLIC :: art_finalize_isorropia
  PUBLIC :: art_finalizeNaCl_isorropia
  PUBLIC :: art_organize_isorropia

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
  SUBROUTINE art_organize_isorropia(jg, tracer, iTRH2SO4, iTRNH3, iTRHNO3, iTRHCL, &
       &                            i_startblk, i_endblk, i_rlstart, i_rlend,      &
       &                            kstart, kend, p_diag, p_prog, p_patch)
    !<
    ! SUBROUTINE art_organize_isorropia_new
    ! Organization of isorropia functionality
    !
    ! Part of Module: mo_art_gas_aero
    ! Author: Heike Vogel, KIT
    !         Sven Werchner, KIT
    ! Initial Release: 2017-XX-XX
    ! Modifications:
    ! YYYY-MM-DD: <name>, <institution>
    ! 2021-06-02: Ulrike Niemeier, MPI-M
    ! - the original routine is extremeley slow. This implementation
    !   provides a much faster alternative
    !>
    
    INTEGER, INTENT(in)               :: &
         &  jg,                          & !< domain index
         &  iTRH2SO4, iTRNH3, iTRHNO3,   & !< Indices of necessary gases
         &  iTRHCL,                      &
         &  i_startblk, i_endblk,        &
         &  i_rlstart, i_rlend,          &         
         &  kstart, kend
    TYPE(t_nh_diag), INTENT(in)       :: &
         &  p_diag                         !< list of diagnostic fields
    TYPE(t_nh_prog), INTENT(in)       :: &
         &  p_prog
    REAL(wp), INTENT(inout)           :: &
         &  tracer(:,:,:,:)                !< tracer-field
    TYPE(t_patch), INTENT(in)         :: & 
         &  p_patch                        !< patch on which computation is performed
    
    !Local variables
    REAL(wp)                          :: &
         &  wi(8), cntrl(2),             & !< Isorropia input vector
         &  aerliq(15),wt(8),gas(3),     &
         &  aersld(19),other(9), rh,     &
         &  tempi,                       &
         &  ah2o, ano3, anh4, ana, acl,  & !< Isorropia output temporaries
         &  sumso4, sumna, sumcl
    INTEGER                           :: &
         &  istart, iend,                & !< bounds for loop over cells in block
         &  jc, jk, jb                     !< Loop indices
    CHARACTER(len=15)                 :: &
         &  scasi
    !!THIS IS RELIMINARY!!!
    REAL(wp), PARAMETER               :: &
         &  mwh2so4  = 0.09807354_wp,    & !< Molecular weight of H2SO4 (kg mol-1)
         &  mwnh3 = 0.0170_wp,           & !< Molecular weight of nh3 (kg mol-1)
         &  mwhno3  = 0.0630049_wp,      & !< Molecular weight of hno3 (kg mol-1)
         &  mwso4 = 0.0960576_wp,        & !< Molecular weight of so4 (kg mol-1)
         &  mwnh4  = 0.01803858_wp,      & !< Molecular weight of nh4 (kg mol-1)
         &  mwno3 = 0.0620049_wp,        & !< Molecular weight of no3 (kg mol-1)
         &  mwcl  = 0.035453_wp,         & !< Molecular weight of cl (kg mol-1)
         &  mwna = 0.0229898_wp,         & !< Molecular weight of na (kg mol-1)
         &  mwhcl = 0.03646_wp,          & !< Molecular weight of hcl (kg mol-1)
         &  mwh2o = 0.018_wp               !< Molecular weight of water (kg mol-1)
    
    
    wi(:)                = 0._wp
    aerliq(:)            = 0._wp
    aersld(:)            = 0._wp
    wt(:)                = 0._wp
    other(:)             = 0._wp
    gas(:)               = 0._wp
    tracer(:,:,:,iTRHCL) = MAX(tracer(:,:,:,iTRHCL),1.e-25_wp)
    
!$omp parallel do default (shared) private(jb, istart, iend)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, istart, iend, i_rlstart, i_rlend)
      DO jk = kstart, kend
!NEC$ ivdep
        DO jc = istart, iend
          ! prepare
          CALL art_prepare_isorropia(jg, tracer(jc,jk,jb,:) , "so4", wi(2))
          sumso4 = wi(2)
          wi(2) = wi(2) * p_prog%rho(jc,jk,jb) / mwso4 * 1.e-09_wp
          CALL art_prepare_isorropia(jg, tracer(jc,jk,jb,:) , "nh4", wi(3))
          wi(3) = wi(3) * p_prog%rho(jc,jk,jb) * 1.e-09_wp / mwnh4                                   &
               &   + tracer(jc,jk,jb,iTRNH3)  * p_prog%rho(jc,jk,jb) / mwnh3
          CALL art_prepare_isorropia(jg, tracer(jc,jk,jb,:) , "no3", wi(4))
          wi(4) = wi(4) * p_prog%rho(jc,jk,jb) * 1.e-09_wp / mwno3                                   &
               &   + tracer(jc,jk,jb,iTRHNO3) * p_prog%rho(jc,jk,jb) / mwhno3
          CALL art_prepare_isorropia(jg, tracer(jc,jk,jb,:) , "na" , wi(1)) ! check the index
          sumna = wi(1)
          wi(1) = wi(1) * p_prog%rho(jc,jk,jb) / mwna * 1.e-09_wp
          CALL art_prepare_isorropia(jg, tracer(jc,jk,jb,:) , "cl" , wi(5))
          sumcl = wi(5)                                   ! sum of cl component in aerosol phase only
          wi(5) = wi(5) * p_prog%rho(jc,jk,jb) * 1.e-09_wp / mwcl                                    &
               &     + tracer(jc,jk,jb,iTRHCL) * p_prog%rho(jc,jk,jb) / mwhcl
          
          rh = rel_hum(p_diag%temp(jc,jk,jb), tracer(jc,jk,jb,iqv), &
               &          p_prog%exner(jc,jk,jb))
          rh = rh * 0.01_wp
          tempi = p_diag%temp(jc,jk,jb)
          cntrl(1) = 0._wp
          cntrl(2) = 1._wp
          IF (wi(2)> 1.e-25_wp) THEN

            CALL isoropia (wi, rh, tempi, cntrl, wt, gas, aerliq, aersld, scasi, other)
            
            ! finalize
            !! update gases
            tracer(jc,jk,jb,iTRNH3)  = gas(1) * mwnh3  / p_prog%rho(jc,jk,jb)
            tracer(jc,jk,jb,iTRHNO3) = gas(2) * mwhno3 / p_prog%rho(jc,jk,jb)
            tracer(jc,jk,jb,iTRHCL)  = gas(3) * mwhcl  / p_prog%rho(jc,jk,jb)
            
            ah2o = aerliq(8) / p_prog%rho(jc,jk,jb) * mwh2o * 1.e+09_wp
            ano3 = aerliq(7) / p_prog%rho(jc,jk,jb) * mwno3 + aerliq(11) / p_prog%rho(jc,jk,jb) * mwhno3
            ano3 = ano3 * 1.e+09_wp
            anh4 = aerliq(3) / p_prog%rho(jc,jk,jb) * mwnh4 + aerliq(9)  / p_prog%rho(jc,jk,jb) * mwnh3
            anh4 = anh4 * 1.e+09_wp
            ana  = aerliq(2) / p_prog%rho(jc,jk,jb) * mwna * 1.e+09_wp
            acl  = aerliq(4) / p_prog%rho(jc,jk,jb) * mwcl * 1.e+09_wp
            
            CALL art_finalize_isorropia    (jg, tracer(jc,jk,jb,:), "h2o", sumso4, ah2o)
            CALL art_finalize_isorropia    (jg, tracer(jc,jk,jb,:), "nh4", sumso4, anh4)
            CALL art_finalize_isorropia    (jg, tracer(jc,jk,jb,:), "no3", sumso4, ano3)
            ! Different redistribution scheme
            CALL art_finalizeNaCl_isorropia(jg, tracer(jc,jk,jb,:), "na" , sumna , ana )
            CALL art_finalizeNaCl_isorropia(jg, tracer(jc,jk,jb,:), "cl" , sumcl , acl )
          ENDIF ! wi(2)>1.e-25_wp
        ENDDO ! jc
      ENDDO ! jk
    ENDDO !jb
!$omp end parallel do
  
END SUBROUTINE art_organize_isorropia
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_prepare_isorropia(jg, tracer, species, wi_entry)
!<
! SUBROUTINE art_prepare_isorropia
! Calculates wi-vector entries for isorropia
!
! Part of Module: mo_art_gas_aero
! Author: Heike Vogel, KIT
!         Sven Werchner, KIT
! Initial Release: 2017-XX-XX
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>

  INTEGER, INTENT(in)            :: &
    &  jg                             !< domain index
  REAL(wp), INTENT(inout)        :: &
    &  tracer(:)                      !< tracer-field
  CHARACTER(LEN=*),INTENT(in)    :: &
    &  species
  REAL(wp), INTENT(out)          :: &
    &  wi_entry

  !Local variables
  TYPE(t_mode), POINTER   :: this_mode
  CHARACTER(LEN=20)       :: sname
  INTEGER                 :: index, ierror
  
  this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
  wi_entry = 0.0_wp
  DO WHILE(ASSOCIATED(this_mode))
    ! Select type of mode
    SELECT TYPE (fields=>this_mode%fields)
      TYPE IS (t_fields_2mom)
        sname = TRIM(species)//'_'//fields%name
        CALL p_art_data(jg)%dict_tracer%get(sname,index,ierror)
        IF (ierror /= SUCCESS) THEN
          this_mode => this_mode%next_mode
          CYCLE
        ENDIF
        wi_entry = wi_entry + tracer(index)  
    END SELECT
    this_mode => this_mode%next_mode
  END DO

END SUBROUTINE art_prepare_isorropia
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_finalize_isorropia(jg, tracer, species, sumso4, amass)
!<
! SUBROUTINE art_finalize_isorropia
! redistributes aerosol compounds
!
! Part of Module: mo_art_gas_aero
! Author: Heike Vogel, KIT
!         Sven Werchner, KIT
! Initial Release: 2017-XX-XX
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>

  INTEGER, INTENT(in)            :: &
    &  jg                             !< domain index
  REAL(wp), INTENT(inout)        :: &
    &  tracer(:)                      !< tracer-field
  CHARACTER(LEN=*),INTENT(in)    :: &
    &  species
  REAL(wp), INTENT(in)          :: &
    &  sumso4, amass

  !Local variables
  TYPE(t_mode), POINTER   :: this_mode
  CHARACTER(LEN=20)       :: sname, so4name
  INTEGER                 :: sindex, so4index, ierror
  
  this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

  DO WHILE(ASSOCIATED(this_mode))
    ! Select type of mode
    SELECT TYPE (fields=>this_mode%fields)
      TYPE IS (t_fields_2mom)
        sname   = TRIM(species)//'_'//fields%name
        so4name = 'so4_'//fields%name
        CALL p_art_data(jg)%dict_tracer%get(sname  ,sindex  ,ierror)
        IF (ierror /= SUCCESS) THEN
          this_mode => this_mode%next_mode
          CYCLE
        ENDIF
        CALL p_art_data(jg)%dict_tracer%get(so4name,so4index,ierror)
        IF (ierror /= SUCCESS) THEN
          this_mode => this_mode%next_mode
          CYCLE
        ENDIF
        tracer(sindex) = tracer(so4index) / sumso4 * amass
    END SELECT
    this_mode => this_mode%next_mode
  END DO
  
    

END SUBROUTINE art_finalize_isorropia
!!
!!-------------------------------------------------------------------------
!!

SUBROUTINE art_finalizeNaCl_isorropia(jg, tracer, species, sum, amass)
!<
! SUBROUTINE art_finalizeNaCl_isorropia
! redistributes Na and Cl
!
! Part of Module: mo_art_gas_aero
! Author: Heike Vogel, KIT
!         Sven Werchner, KIT
! Initial Release: 2017-XX-XX
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>

  INTEGER, INTENT(in)            :: &
    &  jg                             !< domain index
  REAL(wp), INTENT(inout)        :: &
    &  tracer(:)                      !< tracer-field
  CHARACTER(LEN=*),INTENT(in)    :: &
    &  species
  REAL(wp), INTENT(in)          :: &
    &  sum, amass

  !Local variables
  TYPE(t_mode), POINTER   :: this_mode
  CHARACTER(LEN=20)       :: sname
  INTEGER                 :: sindex, ierror
  
  this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

  DO WHILE(ASSOCIATED(this_mode))
    ! Select type of mode
    SELECT TYPE (fields=>this_mode%fields)
      TYPE IS (t_fields_2mom)
        sname   = TRIM(species)//'_'//fields%name
        CALL p_art_data(jg)%dict_tracer%get(sname, sindex, ierror)
        IF (ierror /= SUCCESS) THEN
          this_mode => this_mode%next_mode
          CYCLE
        ENDIF

        tracer(sindex) = tracer(sindex) / sum * amass
    END SELECT
    this_mode => this_mode%next_mode
  END DO
  
    

END SUBROUTINE art_finalizeNaCl_isorropia
!!
!!-------------------------------------------------------------------------
!!

END MODULE mo_art_gas_aero
