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

MODULE mo_icon_fluxes_sw

  USE mo_kind,          ONLY: wp
  USE mo_optical_props, ONLY: ty_optical_props
  USE mo_fluxes,        ONLY: ty_fluxes_broadband
  IMPLICIT NONE

  TYPE, EXTENDS(ty_fluxes_broadband), PUBLIC :: ty_icon_fluxes_sw
    REAL(wp), DIMENSION(:), POINTER :: &
      vis_dn_dir_sfc => NULL(), &
      par_dn_dir_sfc => NULL(), &
      nir_dn_dir_sfc => NULL(), &
      vis_dn_dff_sfc => NULL(), &
      par_dn_dff_sfc => NULL(), &
      nir_dn_dff_sfc => NULL(), &
      vis_up_sfc => NULL(), &
      par_up_sfc => NULL(), &
      nir_up_sfc => NULL()
    REAL(wp), ALLOCATABLE :: &
      band_weight(:), & ! adjustment for current Earth/Sun distance
      frc_par(:), &
      frc_vis(:)
    
  CONTAINS
    PROCEDURE, PUBLIC :: reduce      => reduce_icon
    PROCEDURE, PUBLIC :: are_desired => are_desired_icon
    FINAL             :: del
  END TYPE ty_icon_fluxes_sw

  PUBLIC :: set_fractions

CONTAINS

  SUBROUTINE set_fractions(this, optical_props, solar_constant, ssi_fraction)
    TYPE(ty_icon_fluxes_sw), INTENT(INOUT) :: this
    CLASS(ty_optical_props), INTENT(IN) :: optical_props
    REAL(wp), INTENT(IN) :: solar_constant, ssi_fraction(:)

    INTEGER :: i, nbndsw
    REAL(wp), PARAMETER :: nir_vis_boundary   = 14500._wp
    REAL(wp), ALLOCATABLE :: wavenum(:,:), delwave(:)

    nbndsw = optical_props%get_nband()
    
    IF (.not. ALLOCATED(this%frc_par)) THEN
      ALLOCATE(this%frc_par(nbndsw))
      ALLOCATE(this%frc_vis(nbndsw))
      ALLOCATE(this%band_weight(nbndsw))
      ALLOCATE(wavenum(2,nbndsw))
      ALLOCATE(delwave(nbndsw))
      
      wavenum = optical_props%get_band_lims_wavenumber()
      delwave = wavenum(2,:) - wavenum(1,:)

      this%frc_par(1:nbndsw) = 0.0
      this%frc_par(9) = 0.533725_wp
      this%frc_par(10) = 1.0_wp
      this%frc_par(11) = 0.550164_wp

      DO i = 1, nbndsw 
        this%frc_vis(i) = MAX(0.0_wp, MIN(1.0_wp, &
          (wavenum(2,i) - nir_vis_boundary) / delwave(i) ))
      ENDDO

      ! DA TODO: move to the GPU
      !$ACC ENTER DATA COPYIN(this)
      !$ACC ENTER DATA COPYIN(this%frc_par, this%frc_vis, this%band_weight)
    ENDIF

    ! --- weight radiation within a band for the solar cycle ---
    ! solar_constant contains TSI (the "solar constant") scaled with the
    ! Sun-Earth distance. ssi_fraction contains the relative contribution
    ! of each band to TSI. ssi_default is the (originally only
    ! implicitly defined) solar flux in the 14 bands.

    ! This routine is called from within RRTMGP, so it shouldn't be async
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) FIRSTPRIVATE(nbndsw, solar_constant) &
    !$ACC   COPYIN(ssi_fraction)
    !$ACC LOOP GANG VECTOR
    DO i = 1, nbndsw
      this%band_weight(i) = solar_constant*ssi_fraction( MOD(i, nbndsw)+1 ) ! / ssi_default(:)
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT(1)
  END SUBROUTINE set_fractions

  FUNCTION reduce_icon(this, gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir) result(error_msg)
    CLASS(ty_icon_fluxes_sw),          INTENT(INOUT) :: this
    REAL(kind=wp), DIMENSION(:,:,:),   INTENT(IN   ) :: gpt_flux_up ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    REAL(kind=wp), DIMENSION(:,:,:),   INTENT(IN   ) :: gpt_flux_dn ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    CLASS(ty_optical_props),           INTENT(IN   ) :: spectral_disc  !< derived type with spectral information
    LOGICAL,                           INTENT(IN   ) :: top_at_1
    REAL(kind=wp), DIMENSION(:,:,:), OPTIONAL, &
                                       INTENT(IN   ) :: gpt_flux_dn_dir! Direct flux down
    CHARACTER(len=128)                               :: error_msg
    ! ------
    INTEGER :: nlev, ncol, ngpt, nbndsw, isfc, band, gpt, limits(2), jl
    INTEGER :: band2gpt(2, spectral_disc%get_nband())

    error_msg = this%ty_fluxes_broadband%reduce(gpt_flux_up, gpt_flux_dn, &
      spectral_disc, top_at_1, gpt_flux_dn_dir)
    IF (TRIM(error_msg) /= '') RETURN

    ncol = SIZE(gpt_flux_up,1)
    nlev = SIZE(gpt_flux_up,2)
    ngpt = SIZE(gpt_flux_up,3)
    nbndsw = spectral_disc%get_nband()
    IF (top_at_1) THEN
      isfc = nlev
    ELSE
      isfc = 1
    ENDIF

    !DA TODO: this has to run on GPU
    band2gpt(:,:) = spectral_disc%get_band_lims_gpoint()

    ! This routine is called from within RRTMGP, so it shouldn't be async
    !$ACC WAIT

    !ACCWA Cray CCE (<=16.0.1) needs explicit present clauses
    !$ACC PARALLEL ASYNC(1) DEFAULT(NONE) PRESENT(this) &
    !$ACC   PRESENT(this%vis_dn_dir_sfc, this%par_dn_dir_sfc, this%nir_dn_dir_sfc) &
    !$ACC   PRESENT(this%vis_dn_dff_sfc, this%par_dn_dff_sfc, this%nir_dn_dff_sfc) &
    !$ACC   PRESENT(this%vis_up_sfc, this%par_up_sfc, this%nir_up_sfc) &
    !$ACC   PRESENT(gpt_flux_dn_dir, gpt_flux_up, gpt_flux_dn) &
    !$ACC   FIRSTPRIVATE(ncol, nbndsw, isfc) &
    !$ACC   COPYIN(band2gpt) VECTOR_LENGTH(64)
    !$ACC LOOP GANG VECTOR PRIVATE(limits)
    DO jl = 1, ncol
      this%vis_dn_dir_sfc(jl) = 0.0_wp
      this%par_dn_dir_sfc(jl) = 0.0_wp
      this%nir_dn_dir_sfc(jl) = 0.0_wp
      this%vis_dn_dff_sfc(jl) = 0.0_wp
      this%par_dn_dff_sfc(jl) = 0.0_wp
      this%nir_dn_dff_sfc(jl) = 0.0_wp
      this%vis_up_sfc(jl) = 0.0_wp
      this%par_up_sfc(jl) = 0.0_wp
      this%nir_up_sfc(jl) = 0.0_wp

      !$ACC LOOP SEQ
      DO band = 1, nbndsw
        limits(:) = band2gpt(:, band)
        !$ACC LOOP SEQ
        DO gpt = limits(1), limits(2)
          this%vis_dn_dir_sfc(jl) = this%vis_dn_dir_sfc(jl) + &
            this%frc_vis(band) * gpt_flux_dn_dir(jl,isfc,gpt)
          this%par_dn_dir_sfc(jl) = this%par_dn_dir_sfc(jl) + &
            this%frc_par(band) * gpt_flux_dn_dir(jl,isfc,gpt)
          this%nir_dn_dir_sfc(jl) = this%nir_dn_dir_sfc(jl) + &
            (1.0_wp - this%frc_vis(band)) * gpt_flux_dn_dir(jl,isfc,gpt)

          this%vis_dn_dff_sfc(jl) = this%vis_dn_dff_sfc(jl) + &
            this%frc_vis(band) * ( &
            gpt_flux_dn(jl,isfc,gpt) - gpt_flux_dn_dir(jl,isfc,gpt))
          this%par_dn_dff_sfc(jl) = this%par_dn_dff_sfc(jl) + &
            this%frc_par(band) * ( &
            gpt_flux_dn(jl,isfc,gpt) - gpt_flux_dn_dir(jl,isfc,gpt))
          this%nir_dn_dff_sfc(jl) = this%nir_dn_dff_sfc(jl) + &
            (1.0_wp - this%frc_vis(band)) * ( &
            gpt_flux_dn(jl,isfc,gpt) - gpt_flux_dn_dir(jl,isfc,gpt))

          this%vis_up_sfc(jl) = this%vis_up_sfc(jl) + &
            this%frc_vis(band) * gpt_flux_up(jl,isfc,gpt)
          this%par_up_sfc(jl) = this%par_up_sfc(jl) + &
            this%frc_par(band) * gpt_flux_up(jl,isfc,gpt)
          this%nir_up_sfc(jl) = this%nir_up_sfc(jl) + &
            (1.0_wp - this%frc_vis(band)) * gpt_flux_up(jl,isfc,gpt)
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT(1)
  
  END FUNCTION reduce_icon

  FUNCTION are_desired_icon(this)
    CLASS(ty_icon_fluxes_sw), INTENT(IN) :: this
    LOGICAL                              :: are_desired_icon

    are_desired_icon = this%ty_fluxes_broadband%are_desired()
  END FUNCTION are_desired_icon

  SUBROUTINE del(this)
    TYPE(ty_icon_fluxes_sw), INTENT(INOUT) :: this
    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(this%frc_par, this%frc_vis, this%band_weight)
    !$ACC EXIT DATA DELETE(this)
  END SUBROUTINE del
  ! --------------------------------------------------------------------------------------
END MODULE mo_icon_fluxes_sw
