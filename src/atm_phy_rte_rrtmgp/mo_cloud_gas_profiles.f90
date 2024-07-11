! Module determining gas profiles used in radiative transfer
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

MODULE mo_cloud_gas_profiles

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: max_dom
  USE mo_run_config,           ONLY: ntracer
  USE mo_parallel_config,      ONLY: nproma
  USE mo_aes_phy_config,       ONLY: aes_phy_config
  USE mo_aes_cov_config,       ONLY: aes_cov_config
  USE mo_aes_rad_config,       ONLY: aes_rad_config
  USE mo_coupling_config,      ONLY: is_coupled_to_o3
  USE mo_run_config,           ONLY: iqv, iqc, iqi, iqs, ico2, io3
  USE mo_bc_greenhouse_gases,  ONLY: ghg_co2vmr, ghg_ch4vmr, ghg_n2ovmr, ghg_cfcvmr
  USE mo_bc_ozone,             ONLY: ext_ozone
  USE mo_o3_util,              ONLY: o3_pl2ml, o3_timeint
  USE mo_physical_constants,   ONLY: amd, amw, amco2, amch4, amn2o, amo2, amo3, amc11, amc12
  USE mtime,                   ONLY: datetime
  USE mo_exception,            ONLY: finish

  
  IMPLICIT NONE

  PRIVATE

  PUBLIC              :: gas_profiles, cloud_profiles, snow_profiles, &
                       & init_gas_profiles
  
  INTEGER, PARAMETER  :: ngases=8
  TYPE t_gas
    CHARACTER(LEN=5)     :: name        !< name of gas (chemical composition)
    INTEGER              :: irad        !< integer code for gas profile choice
    ! irad=0: gas VMR is set to zero
    ! irad=1: gas VMR taken from interactively transported tracer of icon
    ! irad=2: gas VMR taken from namelist as stored in
    !         aes_rad_config(:)%vmr_<gas>
    ! irad=3: gas VMR taken from greenhouse gas scenario as given by variables
    !         ghg_<gas>mmr (but is volume mixing ratio!) of
    !         mo_bc_greenhouse_gases
    ! irad=4: time independent profile read from netcdf file will be used
    ! irad=5: time dependent transient profile read from netcdf file and
    !         interpolated with respect to time will be used
    ! irad=6: time dependent but climatological profile read from netcdf file
    !         and interpolated with respect to time will be used
    !         includes an annual cycle but no interannual variability
    ! irad>=12: hyperbolic tangent form of profiles:
    ! irad=12: as irad=2, but modified by a hyperbolic tangent function
    ! irad=13: as irad=3, but modified by a hyperbolic tangent function 
    INTEGER              :: itrac=-999  !< index of tracer in icon 
    REAL(wp)             :: vmr         !< globally constant VMR given in namelist
    REAL(wp)             :: vmr_scenario!< globally constant VMR given by greenouse gas scenario
    REAL(wp)             :: mmr2vmr
    REAL(wp)             :: vpp(3)      !< vertical profile parameters for tanh-profile
    REAL(wp)             :: frad        !< multiplication factor of gas concentration 
  END TYPE t_gas

  TYPE (t_gas), TARGET   :: gas(ngases,max_dom)
  REAL(wp), PARAMETER    :: missing_value=-999999._wp
#ifdef _OPENACC
  LOGICAL, PARAMETER     :: use_acc = .TRUE.
#else  
  LOGICAL, PARAMETER     :: use_acc = .FALSE.
#endif

CONTAINS

  SUBROUTINE init_gas_profiles ()
    ! vertical profile parameters (vpp) of CH4 and N2O
    REAL(wp), PARAMETER :: vpp_ch4(3) = (/1.25e-01_wp,  683.0_wp, -1.43_wp/)
    REAL(wp), PARAMETER :: vpp_n2o(3) = (/1.20e-02_wp, 1395.0_wp, -1.43_wp/)
    INTEGER             :: jg

! general settings (do not depend on time, but on the domain)
    DO jg=1,max_dom
!   H2O
      gas(1,jg)%name = 'h2o  '
      gas(1,jg)%irad = aes_rad_config(jg)% irad_h2o
      gas(1,jg)%vmr  = missing_value
      gas(1,jg)%itrac= iqv
      gas(1,jg)%frad = aes_rad_config(jg)% frad_h2o
      gas(1,jg)%mmr2vmr = amd/amw
!   CO2    
      gas(2,jg)%name = 'co2  '
      gas(2,jg)%irad = aes_rad_config(jg)% irad_co2
      gas(2,jg)%vmr  = aes_rad_config(jg)% vmr_co2
      gas(2,jg)%itrac= MIN(ico2,ntracer)
      gas(2,jg)%frad = aes_rad_config(jg)% frad_co2
      gas(2,jg)%mmr2vmr = amd/amco2
!   CH4
      gas(3,jg)%name = 'ch4  '
      gas(3,jg)%irad = aes_rad_config(jg)% irad_ch4
      gas(3,jg)%vmr  = aes_rad_config(jg)% vmr_ch4
      gas(3,jg)%vpp(:)  = vpp_ch4(:)
      gas(3,jg)%frad = aes_rad_config(jg)% frad_ch4
      gas(3,jg)%mmr2vmr = amd/amch4
!   O2
      gas(4,jg)%name = 'o2   '
      gas(4,jg)%irad = aes_rad_config(jg)% irad_o2
      gas(4,jg)%vmr  = aes_rad_config(jg)% vmr_o2
      gas(4,jg)%frad = aes_rad_config(jg)% frad_o2
      gas(4,jg)%mmr2vmr = amd/amo2
!   O3
      gas(5,jg)%name = 'o3   '
      gas(5,jg)%irad = aes_rad_config(jg)% irad_o3
      gas(5,jg)%vmr  = missing_value
      gas(5,jg)%itrac= MIN(io3,ntracer)
      gas(5,jg)%frad = aes_rad_config(jg)% frad_o3
      gas(5,jg)%mmr2vmr = amd/amo3
      IF (is_coupled_to_o3() .AND. &
        (gas(5,jg)%irad == 5 .OR. gas(5,jg)%irad == 6)) &
         gas(5,jg)%irad = 4 ! one time step from coupler
!   N2O
      gas(6,jg)%name = 'n2o  '
      gas(6,jg)%irad = aes_rad_config(jg)% irad_n2o
      gas(6,jg)%vmr  = aes_rad_config(jg)% vmr_n2o
      gas(6,jg)%vpp(:)  = vpp_n2o(:)
      gas(6,jg)%frad = aes_rad_config(jg)% frad_n2o
      gas(6,jg)%mmr2vmr = amd/amn2o
!   CFC11
      gas(7,jg)%name = 'cfc11'
      gas(7,jg)%irad = aes_rad_config(jg)% irad_cfc11
      gas(7,jg)%vmr  = aes_rad_config(jg)% vmr_cfc11
      gas(7,jg)%frad = aes_rad_config(jg)% frad_cfc11
      gas(7,jg)%mmr2vmr = amd/amc11
!   CFC12
      gas(8,jg)%name = 'cfc12'
      gas(8,jg)%irad = aes_rad_config(jg)% irad_cfc12
      gas(8,jg)%vmr  = aes_rad_config(jg)% vmr_cfc12
      gas(8,jg)%frad = aes_rad_config(jg)% frad_cfc12
      gas(8,jg)%mmr2vmr = amd/amc12
    END DO

    ! vpp(:) values are not used in the OpenACC code, so no need to copy them 
    !$ACC ENTER DATA COPYIN(gas)
  END SUBROUTINE init_gas_profiles

  
  SUBROUTINE gas_profiles ( jg,               jb,               jcs,        &
                          & jce,              kbdim,            klev,       &
                          & ntracer,          this_datetime,    pp_hl,      &
                          & pp_fl,            xq_trc,                       &
                          & xvmr_vap,         xvmr_co2,         xvmr_o3,    &
                          & xvmr_o2,          xvmr_ch4,         xvmr_n2o,   &
                          & xvmr_cfc                                        )

    INTEGER, INTENT(IN) :: jg,      &!> domain index
                         & jb,      &!> block index
                         & jcs,     &!> start index of column (in a block), >= 1
                         & jce,     &!> end index of column (in a block), >=jcs
                         & kbdim,   &!> maximum block length
                         & klev,    &!> number of model levels
                         & ntracer   !> number of tracers
    REAL(wp), INTENT(IN)    :: &
    & pp_hl(kbdim,klev+1),       & !> pressure at half levels [Pa]
    & pp_fl(kbdim,klev),         & !> Pressure at full levels [Pa]
    & xq_trc(kbdim,klev,ntracer)   !> tracer mass fraction [kg/kg], see remark below:
    !
    ! Note: ICON has no sources and sinks in te equation for air density. This implies
    !       that the total air mass is conserved. The parameterized turbulent mass flux
    !       at the surface and precipitation have no effect on the atmospheric mass.
    !       Therefore let let us use here the total air mass as dry air mass.
    !
    TYPE(datetime),POINTER, INTENT(IN)  :: this_datetime !< actual time step
    REAL(wp), INTENT(INOUT) ::      &
    & xvmr_vap(kbdim,klev),         & !< water vapor mass in layer [kg/m2]
    & xvmr_co2(kbdim,klev),         & !< CO2 volume mixing ratio
    & xvmr_o3(kbdim,klev),          & !< O3  volume mixing ratio
    & xvmr_o2(kbdim,klev),          & !< O2  volume mixing ratio
    & xvmr_ch4(kbdim,klev),         & !< CH4 volume mixing ratio
    & xvmr_n2o(kbdim,klev),         & !< N2O volume mixing ratio
    & xvmr_cfc(kbdim,klev,2)          !< CFC volume mixing ratio

    INTEGER             :: igas, jl, jk
    REAL(wp)            :: gas_profile(kbdim,klev,ngases)
    REAL(wp),ALLOCATABLE:: zo3_timint(:,:) !< intermediate value of ozon

    TYPE (t_gas),POINTER:: dom_gas(:) !< gas(:,jg)

    REAL(wp)            :: ghg_cfcvmr1, ghg_cfcvmr2

!   Remark: the order of gases is relevant for this subroutine only.

! general settings (do not depend on time, but on the domain)
    dom_gas => gas(:,jg)

    !$ACC DATA PRESENT(xvmr_vap, xvmr_co2, xvmr_ch4, xvmr_o2, xvmr_o3, xvmr_n2o, xvmr_cfc) &
    !$ACC   PRESENT(xq_trc, dom_gas) &
    !$ACC   CREATE(gas_profile)

! settings depending on time

    ! Pass ghg_cfcvmr as kernel parameters, without explicit copying
    ghg_cfcvmr1 = ghg_cfcvmr(1)
    ghg_cfcvmr2 = ghg_cfcvmr(2)
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
!   CO2
    dom_gas(2)%vmr_scenario = ghg_co2vmr
!   CH4
    dom_gas(3)%vmr_scenario = ghg_ch4vmr
!   N2O
    dom_gas(6)%vmr_scenario = ghg_n2ovmr
!   CFC11
    dom_gas(7)%vmr_scenario = ghg_cfcvmr1
!   CFC12
    dom_gas(8)%vmr_scenario = ghg_cfcvmr2
    !$ACC END KERNELS

    DO igas=1,ngases
      SELECT CASE (dom_gas(igas)%irad)
      CASE (0) ! gas concentration is 0
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
        DO jk=1,klev
          gas_profile(jcs:jce,jk,igas) = 0.0_wp
        END DO
        !$ACC END PARALLEL LOOP
      CASE (1) ! gas is taken from interactive (so transported) tracer
        !note that trasported species are in MMR ...
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
        DO jk=1,klev
          gas_profile(jcs:jce,jk,igas) = xq_trc(jcs:jce,jk,dom_gas(igas)%itrac) * dom_gas(igas)%mmr2vmr
        END DO
        !$ACC END PARALLEL LOOP

      CASE (2,12) ! gas concentration is set from namelist value
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
        DO jk=1,klev
            gas_profile(jcs:jce,jk,igas) = dom_gas(igas)%vmr
        END DO
        !$ACC END PARALLEL LOOP
      CASE (3,13) ! gas concentration is taken from greenhouse dom_gas scenario
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
        DO jk=1,klev
          gas_profile(jcs:jce,jk,igas) = dom_gas(igas)%vmr_scenario
        END DO
        !$ACC END PARALLEL LOOP
      
      !  O3
      CASE (4) ! ozone is constant in time in climatology, first time is used
        IF (igas /= 5) &
          & CALL finish('mo_cloud_gas_profiles: gas_profiles', 'Only ozone is supported for irad = 4')

        CALL o3_pl2ml ( jcs=jcs, jce=jce, kbdim = kbdim,      &
            &          nlev_pres = ext_ozone(jg)%nplev_o3,    &
            &          klev = klev,                           &
            &          pfoz = ext_ozone(jg)%plev_full_o3,     &
            &          phoz = ext_ozone(jg)%plev_half_o3,     &
            &          ppf  = pp_fl(:,:),                     &! in  app1
            &          pph  = pp_hl(:,:),                     &! in  aphp1
            &          o3_time_int = ext_ozone(jg)%o3_plev(:,:,jb,1),&! in
            &          o3_clim     = gas_profile(:,:,igas),   &
            &          opt_use_acc = use_acc                  )

      CASE (5,6)
        IF (igas /= 5) &
          & CALL finish('mo_cloud_gas_profiles: gas_profiles', 'Only ozone is supported for irad = 5 or 6')

        ALLOCATE(zo3_timint(kbdim,ext_ozone(jg)%nplev_o3))
        !$ACC DATA CREATE(zo3_timint)
        CALL o3_timeint(jcs=jcs, jce=jce, kbdim=kbdim,          &
              &          nlev_pres=ext_ozone(jg)%nplev_o3,      &
              &          ext_o3=ext_ozone(jg)%o3_plev(:,:,jb,:),&
              &          current_date=this_datetime,            &
              &          o3_time_int=zo3_timint,                &
              &          opt_use_acc = use_acc                  )
        CALL o3_pl2ml ( jcs=jcs, jce=jce, kbdim = kbdim,        &
              &          nlev_pres = ext_ozone(jg)%nplev_o3,    &
              &          klev = klev,                           &
              &          pfoz = ext_ozone(jg)%plev_full_o3,     &
              &          phoz = ext_ozone(jg)%plev_half_o3,     &
              &          ppf  = pp_fl(:,:),                     &
              &          pph  = pp_hl(:,:),                     &
              &          o3_time_int = zo3_timint,              &
              &          o3_clim     = gas_profile(:,:,igas),   &
              &          opt_use_acc = use_acc                  )
        !$ACC WAIT
        !$ACC END DATA
        DEALLOCATE(zo3_timint)
      END SELECT
    END DO

! For dom_gas(igas)%irad >= 12, multiply with tanh-profile
    DO igas=1,ngases
       IF (dom_gas(igas)%irad .ge. 12) THEN
!          WRITE(0,*) 'igas=',igas,'dom_gas(igas)%name=',dom_gas(igas)%name
!          WRITE(0,*) 'dom_gas(igas)%vpp=',dom_gas(igas)%vpp
         CALL tanh_profile(jcs,       jce,                                &
              &            pp_fl,     dom_gas(igas)%vpp, gas_profile(:,:,igas))
      END IF
    END DO
    ! For all gases, multiply by the frad scaling factor.
    ! This is meant for experiments asking for "4xCO2" e.g.

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
    DO igas=1,ngases
      DO jk=1,klev
        DO jl=jcs,jce
          gas_profile(jl,jk,igas)=MAX(gas_profile(jl,jk,igas)*dom_gas(igas)%frad,0._wp)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP

! Set output fields as asked by icon
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO jk=1,klev
      DO jl=jcs,jce
!       H2O
        xvmr_vap(jl,jk)   = gas_profile(jl,jk,1)
!       CO2    
        xvmr_co2(jl,jk)   = gas_profile(jl,jk,2)
!       CH4
        xvmr_ch4(jl,jk)   = gas_profile(jl,jk,3)
!       O2
        xvmr_o2(jl,jk)    = gas_profile(jl,jk,4)
!       O3
        xvmr_o3(jl,jk)    = gas_profile(jl,jk,5)
!       N2O
        xvmr_n2o(jl,jk)   = gas_profile(jl,jk,6)
!       CFC11
        xvmr_cfc(jl,jk,1) = gas_profile(jl,jk,7)
!       CFC12
        xvmr_cfc(jl,jk,2) = gas_profile(jl,jk,8)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    
    !$ACC WAIT
    !$ACC END DATA
  END SUBROUTINE gas_profiles

  SUBROUTINE cloud_profiles(jg,            jcs,           jce,         &
       &                    klev,          xq_trc,        xm_air,      &
       &                    xm_liq,        xm_ice,        xc_frc       )
    INTEGER, INTENT(IN)    :: jg             ! domain index
    INTEGER, INTENT(IN)    :: jcs            ! start index in block of col.
    INTEGER, INTENT(IN)    :: jce            ! end index in block of col.
    INTEGER, INTENT(IN)    :: klev           ! number of vertical levels
    REAL(wp), INTENT(IN)   :: xq_trc(:,:,:)  ! tracer mass mixing ratios
    REAL(wp), INTENT(IN)   :: xm_air(:,:)    ! air mass in layer
    REAL(wp), INTENT(OUT)  :: xm_liq(:,:)    ! cloud water
    REAL(wp), INTENT(OUT)  :: xm_ice(:,:)    ! cloud ice
    REAL(wp), INTENT(OUT)  :: xc_frc(:,:)    ! cloud fraction in layer

    INTEGER                    :: jl, jk, jks
    REAL(wp)                   :: frad, xrad
    REAL(wp)                   :: cqx

    !$ACC DATA PRESENT(xm_liq, xm_ice, xq_trc, xm_air, xc_frc)
    !
    SELECT CASE (aes_rad_config(jg)%irad_h2o)
    !
    CASE (0) ! no clouds in radiation
      !
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk=1,klev
        DO jl = jcs,jce
          xc_frc(jl,jk) = 0.0_wp
          xm_liq(jl,jk) = 0.0_wp
          xm_ice(jl,jk) = 0.0_wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !
    CASE (1) ! clouds in radiation in layers jks_cloudy to klev, where the mass fraction of
      !        cloud water + cloud ice in air exceeds cqx, with the cloud mass scaled by frad_h2o
      !  
      jks  = aes_phy_config(jg)%jks_cloudy
      cqx  = aes_cov_config(jg)%cqx
      frad = aes_rad_config(jg)%frad_h2o
      !
      IF (frad == 0.0_wp) THEN
         xrad = 0.0_wp
      ELSE
         xrad = 1.0_wp
      END IF
      
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk=1,jks-1 ! no clouds in layers 1 to jks-1
        DO jl = jcs,jce
          xc_frc(jl,jk) = 0.0_wp
          xm_liq(jl,jk) = 0.0_wp
          xm_ice(jl,jk) = 0.0_wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk=jks,klev ! clouds in layers jks to klev
        DO jl = jcs,jce
          !
          ! clouds exist (=1) where the mass fraction of clouds in air exceed cqx, if frad /= 0,
          ! otherwise no clouds in radiation
          xc_frc(jl,jk) = xrad*MERGE(1.0_wp, 0.0_wp, xq_trc(jl,jk,iqc) + xq_trc(jl,jk,iqi) >= cqx)
          !
          ! apply this cloud cover as mask, scale cloud mass in layer by frad, and make >= 0
          xm_liq(jl,jk) = xc_frc(jl,jk)*MAX(xq_trc(jl,jk,iqc)*xm_air(jl,jk)*frad,0._wp) 
          xm_ice(jl,jk) = xc_frc(jl,jk)*MAX(xq_trc(jl,jk,iqi)*xm_air(jl,jk)*frad,0._wp)
          !
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !
    END SELECT
    !
    !$ACC WAIT
    !$ACC END DATA
  END SUBROUTINE cloud_profiles

  SUBROUTINE snow_profiles(jg,            jcs,           jce,         &
       &                   klev,          xq_trc,xm_air, xm_snw       )
    INTEGER, INTENT(IN)    :: jg             ! domain index
    INTEGER, INTENT(IN)    :: jcs            ! start index in block of col.
    INTEGER, INTENT(IN)    :: jce            ! end index in block of col.
    INTEGER, INTENT(IN)    :: klev           ! number of vertical levels
    REAL(wp), INTENT(IN)   :: xq_trc(:,:,:)  ! in  tracer  mass fraction [kg/kg]
    REAL(wp), INTENT(IN)   :: xm_air(:,:)    ! air mass in layer
    REAL(wp), INTENT(OUT)  :: xm_snw(:,:)    ! snow     

    INTEGER                    :: jl, jk
    REAL(wp)                   :: frad

    frad = aes_rad_config(jg)% frad_h2o

    !$ACC DATA PRESENT(xm_snw, xq_trc, xm_air)
    SELECT CASE (aes_rad_config(jg)%irad_h2o)
    CASE (0)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk=1,klev
        DO jl = jcs,jce
            xm_snw(jl,jk) = 0._wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    CASE (1)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk=1,klev
        DO jl = jcs,jce
            xm_snw(jl,jk) = MAX(xq_trc(jl,jk,iqs) * xm_air(jl,jk) * frad, 0._wp)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END SELECT
    !$ACC WAIT
    !$ACC END DATA
  END SUBROUTINE snow_profiles

  SUBROUTINE tanh_profile(jcs,      jce,                 &
       &                  pressure, vpp,     gas_profile )
    INTEGER,  INTENT (IN)    :: jcs              ! start index
    INTEGER,  INTENT (IN)    :: jce              ! end   index and horizontal dimension
    REAL(wp), INTENT (IN)    :: pressure(:,:)    ! pressure at full levels (layer midpoints)
    REAL(wp), INTENT (IN)    :: vpp(3)           ! vertical profile parameters
    REAL(wp), INTENT (INOUT) :: gas_profile(:,:) ! gas profile being modified

    REAL(wp) :: zx_m, zx_p, vpp2, vpp3

    zx_p = 1._wp+vpp(1)
    zx_m = 1._wp-vpp(1)
    vpp2 = vpp(2)
    vpp3 = vpp(3)

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    gas_profile(jcs:jce,:)=(1._wp-(zx_m/zx_p)* &
         & TANH(LOG(pressure(jcs:jce,:)/vpp2)/vpp3)) * &
         & gas_profile(jcs:jce,:) * zx_p * 0.5_wp
    !$ACC END KERNELS
  END SUBROUTINE tanh_profile
END MODULE mo_cloud_gas_profiles
