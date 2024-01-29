! This module contains routines for the vertical interpolation of
! surface/soil fields provided by external analyses to the ICON grid
!
! Soil moisture is read from IFS2ICON as soil moisture index and then
! converted back to soil moisture mass [m] using TERRA soil types.
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

#include "consistent_fma.inc"
MODULE mo_nwp_sfc_interp

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma 
  USE mo_initicon_types,      ONLY: t_init_state
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, ibot_w_so, zml_soil, dzsoil_icon => dzsoil
  USE mo_run_config,          ONLY: msg_level
  USE mo_physical_constants,  ONLY: grav, dtdz_standardatm
  USE sfc_terra_data,         ONLY: cporv, cadp, cfcap, cpwp
  USE mo_ext_data_state,      ONLY: ext_data
  USE mo_exception,           ONLY: finish, message_text, message

  IMPLICIT NONE
  PRIVATE


  ! SUBROUTINE
  PUBLIC :: process_sfcfields
  PUBLIC :: smi_to_wsoil
  PUBLIC :: wsoil_to_smi
  PUBLIC :: wsoil2smi

CONTAINS


  !-------------
  !>
  !! SUBROUTINE process_sfcfields
  !! Routine to convert surface fields interpolated horizontally by IFS2ICON
  !! to the ICON prognostic variables. Important ingredients are 
  !! - height adjustment of temperatures (partly done)
  !! - vertical interpolation of soil temperature and moisture 
  !! - conversion of soil moisture information
  !! - height adjustment of snow cover information (not yet done)
  !!
  !! Other open items - just to document them somewhere
  !! - IFS2ICON needs to take into account land-sea-mask information for horizontal
  !!   interpolation of surface fields (already available)
  !! - And, the most complicated problem, lakes/islands not present at all in the
  !!   source data, is not yet addressed at all here!

   SUBROUTINE process_sfcfields(p_patch, initicon)


    TYPE(t_patch),          INTENT(IN)       :: p_patch
    CLASS(t_init_state),    INTENT(INOUT)    :: initicon

    ! LOCAL VARIABLES
    CHARACTER(LEN=*), PARAMETER       :: routine = 'process_sfcfields'

    INTEGER  :: jg, jb, jk, jc, jk1, idx0(nlev_soil-1)
    INTEGER  :: nlen, nlev, nlev_in, nlevsoil_in

    REAL(wp) :: tcorr1(nproma),tcorr2(nproma),wfac,wfac_vintp(nlev_soil-1),wfac_snow,snowdep

    ! Soil layer center depths in IFS
    REAL(wp), PARAMETER :: zsoil_ifs(4)=(/ 0.035_wp,0.175_wp,0.64_wp,1.945_wp/)

!-------------------------------------------------------------------------

    jg   = p_patch%id
    nlev = p_patch%nlev

    nlev_in     = initicon%atm_in%nlev
    nlevsoil_in = initicon%sfc_in%nlevsoil 

    IF (nlev_in == 0) THEN
      CALL finish(routine, "Number of input levels <nlev_in> not yet initialized.")
    END IF
    IF (nlevsoil_in /= SIZE(zsoil_ifs,1)) THEN
      CALL finish(routine, "Number of soil levels <nlevsoil_in> does not match soil level heights <zsoil_ifs>.")
    END IF

    ! Vertical interpolation indices and weights
!PREVENT_INCONSISTENT_IFORT_FMA
    DO jk = 1, nlev_soil-1
      IF (zml_soil(jk) < zsoil_ifs(1)) THEN
        idx0(jk)       = 0
        wfac_vintp(jk) = 1._wp - zml_soil(jk)/zsoil_ifs(1)
      ELSE IF (zml_soil(jk) > zsoil_ifs(nlevsoil_in)) THEN
        idx0(jk)       = nlevsoil_in
        wfac_vintp(jk) = 1._wp - (zml_soil(jk)-zsoil_ifs(nlevsoil_in))/&
                                 (zml_soil(8) -zsoil_ifs(nlevsoil_in))
      ELSE
        DO jk1 = 1, nlevsoil_in-1
          IF (zml_soil(jk) > zsoil_ifs(jk1) .AND. zml_soil(jk) <= zsoil_ifs(jk1+1)) THEN
            idx0(jk)       = jk1
            wfac_vintp(jk) = (zsoil_ifs(jk1+1)-zml_soil(jk)) / &
                             (zsoil_ifs(jk1+1)-zsoil_ifs(jk1))
          ENDIF
        ENDDO
      ENDIF
    ENDDO

!PREVENT_INCONSISTENT_IFORT_FMA
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen,wfac,tcorr1,tcorr2,snowdep,wfac_snow)
    DO jb = 1, p_patch%nblks_c
      nlen = MERGE(nproma, p_patch%npromz_c, jb /= p_patch%nblks_c)

      ! 2D fields that do not require height adjustment
      ! these fields are simply copied
      ! Note: land-sea-mask is currently unused
      DO jc = 1, nlen
        initicon%sfc%skinres(jc,jb) = initicon%sfc_in%skinres(jc,jb)
        initicon%sfc%ls_mask(jc,jb) = initicon%sfc_in%ls_mask(jc,jb)
        initicon%sfc%seaice(jc,jb) &
          & = MERGE(initicon%sfc_in%seaice(jc,jb), -999.9_wp, &
          &         initicon%sfc_in%seaice(jc,jb) >= 0._wp)
        ! Remark (GZ): the second condition is a workaround until a proper treatment of missing values
        !              becomes available in prep_icon
        initicon%sfc%sst(jc,jb) &
          & = MERGE(initicon%sfc_in%sst(jc,jb), -999.9_wp, &
          &               initicon%sfc_in%sst(jc,jb) > 10._wp &
          &         .AND. initicon%sfc_in%sst(jc,jb) < 305._wp)
      ENDDO

      ! 2D fields that require height adjustment
      ! unfortunately, skin temperature is the only variable for which it is
      ! intuitively clear what to do ...
      DO jc = 1, nlen
        ! Adjust skin temperature with the difference between the atmospheric
        ! temperatures at the lowest model level
        initicon%sfc%tskin(jc,jb)      = initicon%sfc_in%tskin(jc,jb) +       &
          (initicon%atm%temp(jc,nlev,jb) - initicon%atm_in%temp(jc,nlev_in,jb))

        ! Height adjustment for snow variables is not yet implemented
        initicon%sfc%tsnow(jc,jb)    = initicon%sfc_in%tsnow(jc,jb) 
        initicon%sfc%snowweq(jc,jb)  = initicon%sfc_in%snowweq(jc,jb) 
        initicon%sfc%snowdens(jc,jb) = initicon%sfc_in%snowdens(jc,jb) 
        initicon%sfc%snowalb(jc,jb)  = initicon%sfc_in%snowalb(jc,jb) 
      ENDDO

      ! Height adjustment of soil temperatures
      DO jc = 1, nlen
        ! correction for ground level
        tcorr1(jc) = initicon%atm%temp(jc,nlev,jb) - initicon%atm_in%temp(jc,nlev_in,jb)
        ! climatological correction for deep soil levels
        tcorr2(jc) = dtdz_standardatm &
          &        * (initicon%const%topography_c(jc,jb)-initicon%sfc_in%phi(jc,jb)/grav)
      ENDDO

      DO jk = 1, nlevsoil_in
        wfac = REAL(jk,wp)/REAL(nlevsoil_in,wp) ! weight for climatological correction
        DO jc = 1, nlen
          initicon%sfc_in%tsoil(jc,jb,jk) = initicon%sfc_in%tsoil(jc,jb,jk) + &
            wfac*tcorr2(jc) + (1._wp-wfac)*tcorr1(jc)
        ENDDO
      ENDDO

      ! Fill extra levels of incoming data to simplify vertical interpolation
      DO jc = 1, nlen
        ! factor depending on snow depth to weight initialization of soil top temperature
        ! between skin temperature and soil level 1 temperature
        ! For a snow depth of more than 25 cm, it is assumed that the soil top temperature
        ! is the same as the soil level 1 temperature
        snowdep = 1000._wp*initicon%sfc_in%snowweq(jc,jb) / & ! snow depth in m
                  MAX(25._wp,initicon%sfc_in%snowdens(jc,jb))
        wfac_snow = SQRT(MIN(1._wp,4._wp*snowdep))
        initicon%sfc_in%tsoil(jc,jb,0) = (1._wp-wfac_snow)*initicon%sfc%tskin(jc,jb) + &
                                          wfac_snow*initicon%sfc_in%tsoil(jc,jb,1) ! already height-adjusted

        ! Copy climatological deep-soil temperature to extra soil level nlevsoil_in+1
        ! These are limited to -60 deg C because less is definitely nonsense
        initicon%sfc_in%tsoil(jc,jb,nlevsoil_in+1) = MAX(213.15_wp,ext_data(jg)%atm%t_cl(jc,jb))


        initicon%sfc_in%wsoil(jc,jb,0) = initicon%sfc_in%wsoil(jc,jb,1) ! no-gradient condition for moisture

        ! assume no-gradient condition for soil moisture
        initicon%sfc_in%wsoil(jc,jb,nlevsoil_in+1) = initicon%sfc_in%wsoil(jc,jb,nlevsoil_in)
      ENDDO

      ! Vertical interpolation of multi-layer soil fields from IFS levels to TERRA levels
      DO jk = 1, nlev_soil-1
        DO jc = 1, nlen
          initicon%sfc%tsoil(jc,jk,jb) = wfac_vintp(jk) *initicon%sfc_in%tsoil(jc,jb,idx0(jk))+ &
                                  (1._wp-wfac_vintp(jk))*initicon%sfc_in%tsoil(jc,jb,idx0(jk)+1)
          initicon%sfc%wsoil(jc,jk,jb) = wfac_vintp(jk) *initicon%sfc_in%wsoil(jc,jb,idx0(jk))+ &
                                  (1._wp-wfac_vintp(jk))*initicon%sfc_in%wsoil(jc,jb,idx0(jk)+1)
        ENDDO
      ENDDO
      !
      ! Fill top and bottom TERRA levels of tsoil
      DO jc = 1, nlen
        ! copy soil top temperature from incoming data to soil level 0
        initicon%sfc%tsoil(jc,0,jb) = initicon%sfc_in%tsoil(jc,jb,0)
        !
        ! Copy climatological deep-soil temperature to soil level nlev_soil
        ! These are limited to -60 deg C because less is definitely nonsense
        initicon%sfc%tsoil(jc,nlev_soil,jb) = MAX(213.15_wp,ext_data(jg)%atm%t_cl(jc,jb))
      ENDDO

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! convert SMI (from IFS analysis) to W_SO
    !
    CALL smi_to_wsoil(p_patch, initicon%sfc%wsoil)

  END SUBROUTINE process_sfcfields



  !>
  !! SUBROUTINE smi_to_wsoil
  !!
  !! Conversion of soil moisture index into TERRA soil moisture [m]
  !!   soil moisture index = (soil moisture - wilting point) / (field capacity - wilting point)
  !!   safety: min=air dryness point, max=pore volume
  !!
  !! Required input: soil moisture index
  !! Output: soil water content [m H2O]. Input field is overwritten
  !!
  SUBROUTINE smi_to_wsoil(p_patch, wsoil)

    TYPE(t_patch), INTENT(IN)    :: p_patch
    REAL(wp),      INTENT(INOUT) :: wsoil(:,:,:)!soil moisture index in
                                                !soil moisture mass out [m H2O]
                                                ! (nproma,nlev_soil,nblks)
    ! LOCAL VARIABLES
    !
    REAL(wp) :: zwsoil(nproma)        ! local soil moisture field
    INTEGER  :: i_count, ic           ! counter
    INTEGER  :: jg, jb, jk, jc, nlen
    LOGICAL  :: lerr                  ! error flag
    INTEGER  :: ist

    CHARACTER(LEN=*), PARAMETER       :: routine = 'mo_nwp_sfc_interp:smi_to_wsoil'
!-------------

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a)') 'convert SMI to W_SO'
      CALL message('smi_to_wsoil:', TRIM(message_text))
    ENDIF

    jg = p_patch%id

    ! initialize error flag
    lerr = .FALSE.

!$OMP PARALLEL REDUCTION(.or.: lerr)
!$OMP DO PRIVATE(jb,jk,jc,nlen,ic,i_count,zwsoil,ist)
    DO jb = 1, p_patch%nblks_c
      nlen = MERGE(nproma, p_patch%npromz_c, jb /= p_patch%nblks_c)

      ! loop over target (ICON) land points only
      i_count = ext_data(jg)%atm%list_land%ncount(jb)


      ! Conversion of soil moisture index SMI into TERRA soil moisture [m]
      !   soil moisture index = (soil moisture - wilting point) / (field capacity - wilting point)
      !   safety: min=air dryness point, max=pore volume
      ! conversion is only done for hydrological active layers. Remaining layers are filled 
      ! based on a zero gradient assumption.
      DO jk = 1, ibot_w_so

        zwsoil(:) = 0._wp

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count

          jc  = ext_data(jg)%atm%list_land%idx(ic,jb)
          ist = ext_data(jg)%atm%soiltyp(jc,jb)

          ! Catch problematic coast cases: ICON-land but Ocean for source dataset
          ! 
          ! can we do better than this??
          IF ( wsoil(jc,jk,jb) <= -999._wp )  THEN   ! check for missing value
            ! set dummy value: 0.5*(fcap+pwp)
            zwsoil(jc) = 0.5_wp * (cfcap(ist)+cpwp(ist)) * dzsoil_icon(jk)

          ELSE
            SELECT CASE(ist)
            CASE (1,2)  ! ice,rock
              ! set wsoil to 0 for ice and rock
              zwsoil(jc) = 0._wp

            CASE (3,4,5,6,7,8)  ! soil types with non-zero water content
              zwsoil(jc) = dzsoil_icon(jk) * MIN(cporv(ist), &
                & MAX((wsoil(jc,jk,jb)*(cfcap(ist) - cpwp(ist)) + cpwp(ist)),cadp(ist)))

            CASE (9,10) ! sea water, sea ice
              ! ERROR landpoint has soiltype sea water or sea ice
              lerr = .TRUE.

            END SELECT
          END IF
        ENDDO  ! ic
        ! overwrite wsoil
        wsoil(1:nlen,jk,jb) = zwsoil(1:nlen)
      ENDDO  ! jk

      ! assume no-gradient condition for hydraulical non-active layers
      DO jk = ibot_w_so+1, nlev_soil
        DO jc = 1, nlen
          wsoil(jc,jk,jb) = wsoil(jc,jk-1,jb) *  dzsoil_icon(jk)/dzsoil_icon(jk-1)
        ENDDO
      ENDDO  ! jk


    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    IF (lerr) CALL finish(routine, &
         "Landpoint has invalid soiltype (sea water or sea ice)")

  END SUBROUTINE smi_to_wsoil


  !-------------
  !! SUBROUTINE wsoil_to_smi
  !!
  !! Conversion of TERRA soil moisture into soil moisture index
  !!   soil moisture index = (soil moisture - wilting point) / (field capacity - wilting point)
  SUBROUTINE wsoil_to_smi(p_patch, wsoil)


    TYPE(t_patch), INTENT(IN)    :: p_patch
    REAL(wp),      INTENT(INOUT) :: wsoil(:,:,:) ! input: soil moisture mass [m H2O]
                                                 ! output: soil moisture index 

    ! LOCAL VARIABLES
    !
    REAL(wp) :: smi(nproma)           ! local storage for smi
    INTEGER  :: i_count, ic           ! counter
    INTEGER  :: jg, jb, jk, jc, nlen, slt
    LOGICAL  :: lerr                  ! error flag

!-------------

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a)') 'convert W_SO to SMI'
      CALL message('wsoil_to_smi:', TRIM(message_text))
    ENDIF

    jg = p_patch%id
    lerr = .FALSE.
!$OMP PARALLEL REDUCTION(.or.: lerr)
!$OMP DO PRIVATE(jb,jk,jc,nlen,ic,i_count,smi,slt)
    DO jb = 1, p_patch%nblks_c
      nlen = MERGE(nproma, p_patch%npromz_c, jb /= p_patch%nblks_c)

      ! loop over target (ICON) land points only
      i_count = ext_data(jg)%atm%list_land%ncount(jb)

      ! Conversion of TERRA soil moisture [m] into soil moisture index SMI
      !   soil moisture index = (soil moisture - wilting point) / (field capacity - wilting point)
      !   safety: min=air dryness point, max=pore volume
      DO jk = 1, nlev_soil-1

        smi(:) = 0._wp

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data(jg)%atm%list_land%idx(ic,jb)
          slt = ext_data(jg)%atm%soiltyp(jc,jb)
          SELECT CASE(slt)
          CASE (1,2)  !ice,rock
            ! set wsoil to 0 for ice and rock
            smi(jc) = 0._wp
          CASE(3:8)
            ! 3: sand
            ! 4: sandyloam
            ! 5: loam
            ! 6: clayloam
            ! 7: clay
            ! 8: peat
            smi(jc) = (wsoil(jc,jk,jb)/dzsoil_icon(jk) - cpwp(slt)) &
                 / (cfcap(slt) - cpwp(slt))
          CASE(9,10)!sea water, sea ice
            ! ERROR landpoint has soiltype sea water or sea ice
            lerr = .TRUE.
          END SELECT

        ENDDO  ! ic
        ! overwrite wsoil with smi
        wsoil(1:nlen,jk,jb) = smi(1:nlen)
      ENDDO  ! jk


      ! assume no-gradient condition for soil moisture reservoir layer
      DO jc = 1, nlen
        wsoil(jc,nlev_soil,jb) = wsoil(jc,nlev_soil-1,jb)
      ENDDO

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    IF (lerr) CALL finish("wsoil_to_smi", &
         "Landpoint has invalid soiltype (sea water or sea ice)")

  END SUBROUTINE wsoil_to_smi


  SUBROUTINE wsoil2smi(wsoil, dzsoil, soiltyp, smi, ierr)
    !$ACC ROUTINE SEQ
    !
    REAL(wp), INTENT(IN) :: wsoil    !< soil moisture mass [m H2O]
    REAL(wp), INTENT(IN) :: dzsoil   !< soil layer thickness [m]
    INTEGER , INTENT(IN) :: soiltyp  !< soiltype
    REAL(wp), INTENT(OUT):: smi      !< soil moisture index
    INTEGER , INTENT(OUT):: ierr     !< error code

    ierr = 0

    SELECT CASE(soiltyp)
      CASE (1,2)  !ice,rock
      ! set wsoil to 0 for ice and rock
      smi = 0._wp

      CASE(3,4,5,6,7,8)  !sand,sandyloam,loam,clayloam,clay,peat
      smi = (wsoil/dzsoil - cpwp(soiltyp))/(cfcap(soiltyp) - cpwp(soiltyp))

      CASE(9,10)!sea water, sea ice
      ! ERROR landpoint has soiltype sea water or sea ice
      ierr = -1

    END SELECT

  END SUBROUTINE wsoil2smi

END MODULE mo_nwp_sfc_interp
