!
! Provides calls to ART-routines dealing with aerosol dynamics
!
! This module provides an interface to the ART-routines.
! The interface is written in such a way, that ICON will compile and run
! properly, even if the ART-routines are not available at compile time.
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

MODULE mo_art_aerodyn_interface

  USE mo_kind,                          ONLY: wp
  USE mo_fortran_tools,                 ONLY: init, set_acc_host_or_device
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_impl_constants,                ONLY: min_rlcell_int, SUCCESS
  USE mo_impl_constants_grf,            ONLY: grf_bdywidth_c
  USE mo_loopindices,                   ONLY: get_indices_c
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_nonhydro_types,                ONLY: t_nh_diag
  USE mo_run_config,                    ONLY: lart, iqv
  USE mo_parallel_config,               ONLY: nproma
  USE mo_nonhydro_types,                ONLY: t_nh_metrics, t_nh_prog, t_nh_diag
  USE mo_nonhydrostatic_config,         ONLY: kstart_tracer
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag
  USE mo_tracer_metadata_types,         ONLY: t_aero_meta
  USE mtime,                            ONLY: datetime
  USE mo_util_phys,                     ONLY: rel_hum
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_aeroInt
  USE mo_aes_phy_memory,                ONLY: t_aes_phy_tend
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_radio, t_fields_2mom
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_gas_aero,                  ONLY: art_organize_isorropia
  USE mo_art_nucl_aero,                 ONLY: art_nucl_kw
  USE mo_art_cond_aero,                 ONLY: art_prepare_cond_so4, art_finalize_cond_so4
  USE mo_art_aerosol_utilities,         ONLY: art_calc_number_from_mass
  USE mo_art_seas_watercont,            ONLY: art_watercont

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_aerodyn_interface'

  PUBLIC  :: art_aerodyn_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_aerodyn_interface( p_patch,p_dtime,p_prog_list,p_prog, &
  &                               p_metrics,p_diag,tracer, prm_diag, tend, lacc )

  !>
  !! Interface for ART-routines treating aerosol dynamics of any kind
  !! (nucleation, coagulaton, condensation etc.)
  !!
  !! This interface calls the ART-routines, if ICON has been
  !! built including the ART-package. Otherwise, this is simply a dummy
  !! routine.
  !!
  !! @par Revision History
  !! Initial revision by Heike Vogel, KIT (2017-10-19)
  !!
  TYPE(t_patch), TARGET, INTENT(IN) :: &
    &  p_patch                 !< patch on which computation is performed
  REAL(wp), INTENT(IN)    :: &
    &  p_dtime                 !< time step
  TYPE(t_var_list_ptr), INTENT(IN) :: &
    &  p_prog_list             !< current prognostic state list
  TYPE(t_nh_prog), INTENT(IN) :: &
    &  p_prog
  TYPE(t_nh_metrics), INTENT(IN) :: &
    &  p_metrics               !< NH metrics state
  TYPE(t_nh_diag), INTENT(IN) :: &
    &  p_diag                  !< list of diagnostic fields
  REAL(wp), INTENT(INOUT) :: &
    &  tracer(:,:,:,:)         !< tracer mixing ratios (specific concentrations)
  TYPE(t_nwp_phy_diag),OPTIONAL, INTENT(IN)  :: &
    &  prm_diag                !< NH metrics state
  TYPE(t_aes_phy_tend) , OPTIONAL,  POINTER  :: tend
  LOGICAL, OPTIONAL, INTENT(IN) :: lacc
! Local variables
  REAL(wp), ALLOCATABLE   :: &
    &  rh(:,:),              & !< Relative humidty (%)
    &  totmasscond(:,:,:),   & !< Mass change of mode due to condensation
!    &  ch2so4(:,:),          & !< TEMP! concentration of H2SO4
    &  nucmass(:,:),         & !< Nucleated sulfate mass (mug kg-1)
    &  vseas(:,:,:),         & !< Nacl concentration of mode i
    &  vso4seas(:,:,:),      & !< so4 concentration of mode i
    &  vh2oseas(:,:,:)         !< h2o water content of mode i
  REAL(wp)                :: &
    &  nucnmb                  !< Nucleated sulfate number (kg-1)
  INTEGER                 :: &
    &  jb, jc, jk,           & !< loop index
    &  jg,                   & !< domain index
    &  i_startblk, i_endblk, & !< Start and end of block loop
    &  istart, iend,         & !< Start and end of nproma loop
    &  i_rlstart, i_rlend,   & !< Relaxation start and end
    &  i_nchdom,             & !< Number of child domains
    &  kstart,               & !< top level for all aerosol processes (transport, dynamics, ...)
    &  nlev,                 & !< Number of levels (equals index of lowest full level)
    &  tindex,               & !< Index of tracer
    &  ierror                  !< Error return value
  INTEGER                 :: &
    &  iTRH2SO4,             & !< index H2SO4
    &  iTRNH3,               & !< index NH3
    &  iTRHNO3,              & !< index HNO3
    &  iTRHCL
  REAL(wp), PARAMETER     :: &
    &  mwh2so4  = 0.09807354_wp      !< Molecular weight of H2SO4 (kg mol-1)

  CHARACTER(LEN=20)       :: sname

  TYPE(t_mode), POINTER   :: this_mode

  LOGICAL                 :: lzacc             ! OpenACC flag

  !-----------------------------------------------------------------------

  CALL set_acc_host_or_device(lzacc, lacc)

  IF(lart) THEN

    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_aeroInt)

    ! --- Get the loop indizes
    i_nchdom   = MAX(1,p_patch%n_childdom)
    jg         = p_patch%id
    nlev       = p_patch%nlev
    i_rlstart  = grf_bdywidth_c+1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

    IF (art_config(jg)%lart_chem) THEN
      CALL p_art_data(jg)%dict_tracer%get('TRH2SO4',iTRH2SO4,ierror)
      IF (ierror /= SUCCESS) iTRH2SO4    = 0
      CALL p_art_data(jg)%dict_tracer%get('TRNH3',iTRNH3,ierror)
      IF (ierror /= SUCCESS) iTRNH3    = 0
      CALL p_art_data(jg)%dict_tracer%get('TRHNO3',iTRHNO3,ierror)
      IF (ierror /= SUCCESS) iTRHNO3    = 0
      CALL p_art_data(jg)%dict_tracer%get('TRHCL',iTRHCL,ierror)
      IF (ierror /= SUCCESS) iTRHCL    = 0
    ENDIF

    IF (art_config(jg)%lart_aerosol) THEN
! Aerosol dynamics
      ALLOCATE(rh(nproma,nlev))
      ALLOCATE(totmasscond(nproma,nlev,p_patch%nblks_c))
      ALLOCATE(nucmass(nproma,nlev))
!      ALLOCATE(ch2so4(nproma,nlev))  ! TEMP
      ALLOCATE(vseas(nproma,nlev,p_patch%nblks_c))            
      ALLOCATE(vso4seas(nproma,nlev,p_patch%nblks_c))             
      ALLOCATE(vh2oseas(nproma,nlev,p_patch%nblks_c))                   

      !$ACC ENTER DATA CREATE(rh, totmasscond) IF(lzacc)
      CALL init(rh, lacc=lzacc, opt_acc_async=.TRUE.)
      CALL init(totmasscond, lacc=lzacc, opt_acc_async=.TRUE.)
      kstart = MINVAL(kstart_tracer(jg,p_art_data(jg)%aero%itr_start:p_art_data(jg)%aero%itr_end))

      ! ISORROPIA
      !emisselem = p_art_data(jg)%emiss%get_element_emiss(iTRHCL,emiss_elem,ierror )
      !  mol_weight = emisselem%mol_weight
      IF (art_config(jg)%iart_isorropia > 0) THEN
        CALL art_organize_isorropia(jg, tracer, iTRH2SO4, iTRNH3, iTRHNO3, iTRHCL, &
             &                      i_startblk, i_endblk, i_rlstart, i_rlend,      &
             &                      kstart, nlev, p_diag, p_prog, p_patch)
      ENDIF

      IF (art_config(jg)%iart_seas_water > 0) THEN
! water content of seasalt aerosol
        NULLIFY(this_mode)
        this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

        DO WHILE(ASSOCIATED(this_mode))
          ! Select type of mode
          SELECT TYPE (fields=>this_mode%fields)
            TYPE IS (t_fields_2mom)
                  
              IF(fields%info%l_watercont) THEN
                vseas = 1.e-19_wp
                vso4seas = 1.e-19_wp
                vh2oseas = 1.e-19_wp

                sname = TRIM("nacl"//'_'//fields%name)
                CALL p_art_data(jg)%dict_tracer%get(sname,tindex,ierror)
                IF (ierror == SUCCESS) vseas = tracer(:,:,:,tindex)
                sname = TRIM("so4"//'_'//fields%name)
                CALL p_art_data(jg)%dict_tracer%get(sname,tindex,ierror)
                IF (ierror == SUCCESS) vso4seas = tracer(:,:,:,tindex)

                DO jb = i_startblk, i_endblk
                  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                    &                istart, iend, i_rlstart, i_rlend)
                  rh(:,:) = rel_hum(p_diag%temp(:,:,jb), tracer(:,:,jb,iqv), &
                      &                 p_prog%exner(:,:,jb))

                  CALL art_watercont(rh,p_diag%temp(:,:,jb),p_prog%rho(:,:,jb),vseas(:,:,jb),    &
                    &                   vso4seas(:,:,jb),vh2oseas(:,:,jb),                       &
                    &                   istart,iend,kstart,nlev)

                ENDDO !jb

                sname = TRIM("h2o"//'_'//fields%name)
                CALL p_art_data(jg)%dict_tracer%get(sname,tindex,ierror)
                IF (ierror == SUCCESS) tracer(:,:,:,tindex) = vh2oseas(:,:,:)
            END IF
          END SELECT
          this_mode => this_mode%next_mode
        ENDDO !associated(this_mode)

      ENDIF


! Preparation routines
      NULLIFY(this_mode)
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)
          TYPE IS (t_fields_2mom)
! Condensation
            SELECT CASE(fields%info%icondensation)
              CASE(1)
                DO jb = i_startblk, i_endblk
                  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                    &                istart, iend, i_rlstart, i_rlend)

                  ! Update diameter
                  CALL fields%modal_param(p_art_data(jg)%air_prop%art_free_path(:,:,jb),          &
                    &                     istart, iend, kstart, nlev, jb, tracer(:,:,jb,:))

                  ! TEMP:
!                  DO jc = istart, iend
!                    DO jk = kstart, nlev
                    ! Kerminen and Wexler result for 273.15 K and 90 % relative humidity
!                    ch2so4(jc,jk) = 0.004665_wp / p_prog%rho(jc,jk,jb) * 1000._wp
!                    ENDDO
!                  ENDDO
                  !END TEMP
                  CALL art_prepare_cond_so4(p_diag%temp(:,:,jb), p_diag%pres(:,:,jb),    &
                    &                       p_prog%rho(:,:,jb),                          &
                    &                       tracer(:,:,jb,iTRH2SO4),                     &
!                    &                       ch2so4(:,:),                                 & !TEMP
                    &                       tracer(:,:,jb,fields%itr0),                  &
                    &                       fields%diameter(:,:,jb), fields%info%sg_ini, &
                    &                       istart, iend, kstart, nlev,                  &
                    &                       fields%dmdt_condso4(:,:,jb))
                  DO jk = kstart, nlev ! calc total condensated h2so4
!NEC$ ivdep
                    DO jc = istart, iend
                      totmasscond(jc,jk,jb) = totmasscond(jc,jk,jb)                       &
                        &                   + fields%dmdt_condso4(jc,jk,jb) * p_dtime
                    ENDDO !jc
                  ENDDO !jk

                ENDDO !jb
            END SELECT
        END SELECT
        this_mode => this_mode%next_mode
      ENDDO !associated(this_mode)
      NULLIFY(this_mode)

! Aerosol dynamics
! Updating routines
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

! Condensation
      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)
          TYPE IS (t_fields_2mom)
            SELECT CASE(fields%info%icondensation)
              CASE(1)
                DO jb = i_startblk, i_endblk
                  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                    &                istart, iend, i_rlstart, i_rlend)

                  ! TEMP:
!                  DO jc = istart, iend
!                    DO jk = kstart, nlev
!                    ! Kerminen and Wexler result for 273.15 K and 90 % relative humidity
!                      ch2so4(jc,jk) = 0.004665_wp / p_prog%rho(jc,jk,jb) * 1000._wp
!                    ENDDO
!                  ENDDO
                  !END TEMP

                  CALL art_finalize_cond_so4(fields%dmdt_condso4(:,:,jb),                         &
                    &                        tracer(:,:,jb,iTRH2SO4),                             &
!                    &                       ch2so4(:,:),                                          & !TEMP
                    &                        totmasscond(:,:,jb),                                 &
                    &                        tracer(:,:,jb,fields%info%itrcond),                  &
                    &                        p_dtime, istart, iend, kstart, nlev)
                ENDDO !jb
            END SELECT
        END SELECT
        this_mode => this_mode%next_mode
      ENDDO !associated(this_mode)

! Nucleation
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)
          TYPE IS (t_fields_2mom)
            IF (fields%info%inucleation == 1) THEN

              DO jb = i_startblk, i_endblk
                CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                  &                istart, iend, i_rlstart, i_rlend)

                DO jk = kstart, nlev
!NEC$ ivdep
                  DO jc = istart, iend
                    rh(jc,jk) = rel_hum(p_diag%temp(jc,jk,jb), tracer(jc,jk,jb,iqv), &
                      &                 p_prog%exner(jc,jk,jb))
                    ! TEMP:
                    ! Kerminen and Wexler result for 273.15 K and 90 % relative humidity
!                    ch2so4(jc,jk) = 0.004665_wp / p_prog%rho(jc,jk,jb)
                    ! END TEMP
                  ENDDO
                ENDDO

                CALL art_nucl_kw(rh(:,:),p_diag%temp(:,:,jb),p_prog%rho(:,:,jb), &
                  &              tracer(:,:,jb,iTRH2SO4),istart,iend,kstart,nlev,     &
                  &              nucmass(:,:),totmasscond(:,:,jb))

                DO jk = kstart, nlev
!NEC$ ivdep
                  DO jc = istart, iend
                    tracer(jc,jk,jb,fields%info%itrnucl) = tracer(jc,jk,jb,fields%info%itrnucl) &
                      &                                  + nucmass(jc,jk)
                    CALL art_calc_number_from_mass(fields%info%diameter_ini_mass,               &
                      &                            fields%info%sg_ini, p_prog_list,             &
                      &                            fields%info%itrnucl, nucmass(jc,jk), nucnmb)
                    tracer(jc,jk,jb,fields%itr0) = tracer(jc,jk,jb,fields%itr0) &
                      &                          + nucnmb
                  ENDDO
                ENDDO
              ENDDO
              EXIT
            ENDIF
        END SELECT
        this_mode => this_mode%next_mode
      ENDDO !associated(this_mode)

! Mode shift
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)
          TYPE IS (t_fields_2mom)

            IF (art_config(jg)%iart_modeshift==1) THEN
              DO jb = i_startblk, i_endblk
                CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                  &                istart, iend, i_rlstart, i_rlend)
                ! Update diameter
                CALL fields%modal_param(p_art_data(jg)%air_prop%art_free_path(:,:,jb),            &
                  &                     istart, iend, kstart, nlev, jb, tracer(:,:,jb,:))
                CALL fields%mode_shift(istart,iend,kstart,nlev,jb,tracer(:,:,jb,:))
              ENDDO !jb
            ENDIF ! iart_modeshift
        END SELECT
        this_mode => this_mode%next_mode
      ENDDO !associated(this_mode)

      NULLIFY(this_mode)

      !$ACC WAIT
      !$ACC EXIT DATA DELETE(rh, totmasscond) IF(lzacc)
      DEALLOCATE(rh)
      DEALLOCATE(totmasscond)
      DEALLOCATE(nucmass)
      DEALLOCATE(vseas)
      DEALLOCATE(vso4seas)
      DEALLOCATE(vh2oseas)

    ENDIF !lart_aerosol

    IF (timers_level > 3) CALL timer_stop(timer_art_aeroInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  ENDIF !lart


END SUBROUTINE art_aerodyn_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_aerodyn_interface
