!
! Subroutine aes_phy_diag contains small diagnostic routines,
!  which are executed on a single block of cells.
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

MODULE mo_aes_phy_diag

  USE mo_kind,                ONLY: wp

  USE mo_aes_phy_dims,        ONLY: aes_phy_dims
  USE mo_aes_phy_memory,      ONLY: t_aes_phy_field, prm_field,     &
    &                               t_aes_phy_tend,  prm_tend,      &
    &                               cdimissval

  USE mo_physical_constants,  ONLY: cpd, cpv, cvd, cvv, Tf, tmelt
  USE mo_run_config,          ONLY: iqv
  USE mo_aes_cop_config,      ONLY: aes_cop_config
  USE mo_aes_vdf_config,      ONLY: aes_vdf_config
  USE mo_aes_sfc_indices,     ONLY: nsfc_type, iwtr, iice, ilnd

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: surface_fractions, &
       &    droplet_number,    &
       &    cpair_cvair_qconv, &
       &    initialize,        &
       &    finalize

CONTAINS

  !---------------------------------------------------------------------
  SUBROUTINE surface_fractions(jg, jb, jcs, jce)

    ! Arguments
    !
    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    ! Local variables
    !
    TYPE(t_aes_phy_field), POINTER :: field
    
    REAL(wp) :: zfrw (aes_phy_dims(jg)%nproma) !< cell area fraction of open water
    REAL(wp) :: zfri (aes_phy_dims(jg)%nproma) !< cell area fraction of ice covered water
    REAL(wp) :: zfrl (aes_phy_dims(jg)%nproma) !< cell area fraction of land
    INTEGER  :: jc

    field => prm_field(jg)
 
    ! 3.3 Weighting factors for fractional surface coverage
    !     Accumulate ice portion for diagnostics

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) CREATE(zfrw, zfri, zfrl)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc=jcs,jce

      ! fraction of solid land in the grid box, i.e. land without lakes if lakes are used
      ! see mo_aes_phy_init or input data set for details.
      !
      IF (ilnd.LE.nsfc_type) THEN
         zfrl(jc) = MAX(MIN(            &
              &     field%lsmask(jc,jb) &
              &     ,1._wp),0._wp)
      ELSE
         zfrl(jc) = 0._wp
      END IF

      ! fraction of open water in the grid box, for sea and lakes
      !
      IF (iwtr.LE.nsfc_type) THEN
         zfrw(jc) = MAX(MIN(                                                                          &
              &      (1._wp-field%lsmask(jc,jb)-field%alake(jc,jb))*(1._wp-field%seaice(jc,jb))       & ! ocean
              &     +field%alake(jc,jb)                            *(1._wp-field%lake_ice_frc(jc,jb)) & ! lakes
              &     ,1._wp),0._wp)
         !
         ! security for water temperature with changing ice mask
         ! (over lakes; over ocean this is not an issue since, over ocean,
         ! ts_tile(iwtr) is overwritten again with SST from ocean in 
         ! coupling interface after update_surface)
         IF (zfrw(jc) > 0._wp .AND. field%ts_tile(jc,jb,iwtr) == cdimissval) THEN
           ! lake was completely frozen in previous time step but only partially 
           ! frozen in current time step
           field%ts_tile(jc,jb,iwtr) = tmelt
         END IF
      ELSE
         zfrw(jc) = 0._wp
      END IF

      ! fraction of ice covered water in the grid box, for sea and lakes
      !
      IF (iice.LE.nsfc_type) THEN
         zfri(jc) = 1._wp-zfrl(jc)-zfrw(jc)
         !
         ! security for ice temperature with changing ice mask
         IF(zfri(jc) > 0._wp .AND. field%ts_tile(jc,jb,iice) == cdimissval ) THEN
            field%ts_tile(jc,jb,iice)  = tmelt + Tf    ! = 271.35 K
         END IF
      ELSE
         zfri(jc) = 0._wp
      END IF

    END DO

    ! 3.4 Merge three pieces of information into one array for vdiff
    IF (ilnd.LE.nsfc_type) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc=jcs,jce
        field%frac_tile(jc,jb,ilnd) = zfrl(jc)
      END DO
    END IF
    IF (iwtr.LE.nsfc_type) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc=jcs,jce
        field%frac_tile(jc,jb,iwtr) = zfrw(jc)
      END DO
    END IF
    IF (iice.LE.nsfc_type) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc=jcs,jce
        field%frac_tile(jc,jb,iice) = zfri(jc)
      END DO
    END IF
    !$ACC END PARALLEL

    !$ACC WAIT(1)

    NULLIFY(field)

  END SUBROUTINE surface_fractions
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE droplet_number(jg, jb, jcs, jce)

    ! Arguments
    !
    INTEGER, INTENT(in)  :: jg, jb, jcs, jce

    ! Local variables
    !
    INTEGER :: nlev
    !
    TYPE(t_aes_phy_field), POINTER :: field
    !
    INTEGER  :: jc, jk
    LOGICAL  :: lland(aes_phy_dims(jg)%nproma)
    LOGICAL  :: lglac(aes_phy_dims(jg)%nproma)
    REAL(wp) :: zprat, zn1, zn2, zcdnc

    ! Shortcuts to components of aes_cop_config
    !
    REAL(wp) :: cn1lnd, cn2lnd, cn1sea, cn2sea

    nlev   = aes_phy_dims(jg)%nlev

    cn1lnd = aes_cop_config(jg)% cn1lnd
    cn2lnd = aes_cop_config(jg)% cn2lnd
    cn1sea = aes_cop_config(jg)% cn1sea
    cn2sea = aes_cop_config(jg)% cn2sea

    field  => prm_field(jg)

    !$ACC DATA PRESENT(field) &
    !$ACC   CREATE(lland, lglac)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc=jcs,jce
      lland(jc) = field%sftlf (jc,jb) > 0._wp
      lglac(jc) = field%sftgif(jc,jb) > 0._wp
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(zprat, zn1, zn2, zcdnc)
    DO jk = 1,nlev
      DO jc = jcs,jce
        !
        zprat=(MIN(8._wp,80000._wp/field%pfull(jc,jk,jb)))**2

        IF (lland(jc).AND.(.NOT.lglac(jc))) THEN
          zn1= cn1lnd
          zn2= cn2lnd
        ELSE
          zn1= cn1sea
          zn2= cn2sea
        END IF
        IF (field%pfull(jc,jk,jb).LT.80000._wp) THEN
          zcdnc=1.e6_wp*(zn1+(zn2-zn1)*(EXP(1._wp-zprat)))
        ELSE
          zcdnc=zn2*1.e6_wp
        END IF
        field% acdnc(jc,jk,jb) = zcdnc
        !
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

    NULLIFY(field)

  END SUBROUTINE droplet_number
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE cpair_cvair_qconv(jg, jb, jcs, jce)

    ! Arguments
    !
    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    ! Local variables
    !
    INTEGER :: nlev
    !
    TYPE(t_aes_phy_field), POINTER :: field
    !
    INTEGER :: jc, jk

    nlev   = aes_phy_dims(jg)%nlev

    field  => prm_field(jg)
    
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,nlev
      DO jc=jcs,jce
        field%cpair(jc,jk,jb) = cpd+(cpv-cpd)*field%qtrc_phy(jc,jk,jb,iqv)
        field%cvair(jc,jk,jb) = cvd+(cvv-cvd)*field%qtrc_phy(jc,jk,jb,iqv)
        !
        field%qconv(jc,jk,jb) = 1._wp/(field%mair(jc,jk,jb)*field%cpair(jc,jk,jb))
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT(1)
    
    NULLIFY(field)

  END SUBROUTINE cpair_cvair_qconv
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE initialize   (jg, jb, jcs, jce)

    ! Arguments
    !
    INTEGER, INTENT(in) :: jg, jb, jcs, jce

    ! Local variables
    !
    INTEGER :: nlev
    !
    TYPE(t_aes_phy_field), POINTER      :: field
    INTEGER                             :: jc, jk

    nlev   = aes_phy_dims(jg)%nlev

    field => prm_field(jg)
    
    IF (ASSOCIATED(field% q_phy)) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG
      DO jk = 1,nlev
        !$ACC LOOP VECTOR
        DO jc=jcs,jce
          field% q_phy(jc,jk,jb) = 0._wp
        END DO
      END DO
      !$ACC END PARALLEL
    END IF
    IF (ASSOCIATED(field% q_phy_vi)) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc=jcs,jce
        field% q_phy_vi(jc, jb) = 0._wp
      END DO
      !$ACC END PARALLEL
    END IF

    !$ACC WAIT(1)

    NULLIFY(field)

  END SUBROUTINE initialize
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE finalize     (jg, jb, jcs, jce)

    ! Arguments
    !
    INTEGER, INTENT(in) :: jg, jb, jcs, jce

    ! Local variables
    !
    INTEGER :: nlev
    !
    TYPE(t_aes_phy_field), POINTER :: field
    TYPE(t_aes_phy_tend) , POINTER :: tend
    !
    INTEGER :: jc, jk

    nlev  = aes_phy_dims(jg)%nlev

    field => prm_field(jg)
    tend  => prm_tend (jg)

    !$ACC DATA PRESENT(tend%ta_phy, field%cpair, field%cvair)

    IF (.NOT. aes_vdf_config(jg)%use_tmx) THEN
      ! convert the temperature tendency from physics, as computed for constant pressure conditions,
      ! to constant volume conditions, as needed for the coupling to the dynamics
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1,nlev
        DO jc=jcs,jce
          tend% ta_phy(jc,jk,jb) = tend% ta_phy(jc,jk,jb) * field% cpair(jc,jk,jb) / field% cvair(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

    !$ACC END DATA
    !$ACC WAIT(1)

    NULLIFY(field)
    NULLIFY(tend )

  END SUBROUTINE finalize
  !---------------------------------------------------------------------

END MODULE mo_aes_phy_diag
