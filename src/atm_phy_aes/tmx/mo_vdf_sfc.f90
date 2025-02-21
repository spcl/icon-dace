!
! Classes and functions for the turbulent mixing package (tmx)
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

MODULE mo_vdf_sfc

  USE mo_kind,              ONLY: wp, sp, i1, i4
  USE mo_exception,         ONLY: message, finish
  USE mo_fortran_tools,     ONLY: init, copy
  USE mtime,                ONLY: t_datetime => datetime
  USE mo_master_config,     ONLY: isrestart  ! TODO: use config instead of USE
  USE mo_tmx_process_class, ONLY: t_tmx_process
  USE mo_tmx_field_class,   ONLY: t_tmx_field, t_domain, isfc_oce, isfc_ice, isfc_lnd
  USE mo_variable,          ONLY: t_variable
  USE mo_variable_list,     ONLY: t_variable_list, t_variable_set


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_vdf_sfc, t_vdf_sfc_inputs, t_vdf_sfc_config, t_vdf_sfc_diagnostics, average_tiles

  TYPE, EXTENDS(t_tmx_process) :: t_vdf_sfc
  CONTAINS
    PROCEDURE :: Compute
    PROCEDURE :: Compute_diagnostics
    PROCEDURE :: Update_diagnostics
  END TYPE t_vdf_sfc
  
  INTERFACE t_vdf_sfc
    MODULE PROCEDURE t_vdf_sfc_construct
  END INTERFACE

  TYPE, EXTENDS(t_variable_set) :: t_vdf_sfc_config
    REAL(wp), POINTER :: &
      & dtime => NULL(), &
      & cpd => NULL(), &
      & cvd => NULL(), &
      & min_sfc_wind => NULL(), &
      & min_rough    => NULL(), &
      & rough_m_oce  => NULL(), &
      & rough_m_ice  => NULL(), &
      & fsl          => NULL()    ! weight for interpolation to surface_layer mid level
    INTEGER, POINTER :: &
      & nice_thickness_classes => NULL()

  CONTAINS
    ! PROCEDURE :: Init => init_t_vdf_sfc_variable_set
    PROCEDURE :: Set_pointers => Set_pointers_config
  END TYPE t_vdf_sfc_config

  TYPE, EXTENDS(t_variable_set) :: t_vdf_sfc_inputs
    REAL(wp), POINTER :: &
      & ta(:,:) => NULL(), &
      & tv(:,:) => NULL(), &
      & ua(:,:) => NULL(), &
      & va(:,:) => NULL(), &
      & qa(:,:) => NULL(), &
      & pa(:,:) => NULL(), &
      & psfc(:,:) => NULL(), &
      & rsfl(:,:) => NULL(), &
      & ssfl(:,:) => NULL(), &
      & rlds(:,:)     => NULL(), &
      & rsds(:,:)     => NULL(), &
      & rvds_dir(:,:) => NULL(), &
      & rnds_dir(:,:) => NULL(), &
      & rpds_dir(:,:) => NULL(), &
      & rvds_dif(:,:) => NULL(), &
      & rnds_dif(:,:) => NULL(), &
      & rpds_dif(:,:) => NULL(), &
      & emissivity(:,:) => NULL(), &
      & cosmu0(:,:)   => NULL(), &
      & co2(:,:)      => NULL(), &
      & tsfc_tile(:,:,:) => NULL(), &
      ! & rough_h(:,:,:) => NULL(), &
      ! & rough_m(:,:,:) => NULL()
      & u_oce_current(:,:) => NULL(), &
      & v_oce_current(:,:) => NULL(), &
      & ice_thickness(:,:)   => NULL(),  & !< sea ice: ice thickness [m]
      & zf(:,:) => NULL(), & !< geom. height of lowest atm. full level [m]
      & zh(:,:) => NULL(), & !< geom. height of surface (interface level) [m]
      & fract_tile(:,:,:) => NULL()
    REAL(wp), POINTER :: &
      & dz(:,:) => NULL()


  CONTAINS
    ! PROCEDURE :: Init => init_t_vdf_sfc_variable_set
    PROCEDURE :: Set_pointers => Set_pointers_inputs
  END TYPE t_vdf_sfc_inputs

  TYPE, EXTENDS(t_variable_set) :: t_vdf_sfc_diagnostics
    REAL(wp), POINTER :: &
      & wind(:,:) => NULL(), &
      & theta_atm(:,:) => NULL(), &
      & thetav_atm(:,:) => NULL(), &
      & theta_tile(:,:,:) => NULL(), &
      & thetav_tile(:,:,:) => NULL(), &
      & rough_h_tile(:,:,:) => NULL(), &
      & rough_m_tile(:,:,:) => NULL(), &
      & rough_h(:,:) => NULL(), &
      & rough_m(:,:) => NULL(), &
      & qsat_tile(:,:,:) => NULL(), &
      & evapotrans_tile(:,:,:) => NULL(), &
      & lhfl_tile(:,:,:) => NULL(), &
      & shfl_tile(:,:,:) => NULL(), &
      & q_snocpymlt_lnd(:,:) => NULL(), &
      & ustress_tile(:,:,:) => NULL(), &
      & vstress_tile(:,:,:) => NULL(), &
      & evapotrans(:,:) => NULL(), &
      & lhfl(:,:) => NULL(), &
      & shfl(:,:) => NULL(), &
      & ustress(:,:) => NULL(), &
      & vstress(:,:) => NULL(), &
      & tsfc(:,:) => NULL(), &
      & tsfc_rad(:,:) => NULL(), &
      & lwfl_net_tile(:,:,:) => NULL(), &
      & swfl_net_tile(:,:,:) => NULL(), &
      & lwfl_up(:,:)         => NULL(), &
      & swfl_up(:,:)         => NULL(), &
      & rho_tile(:,:,:) => NULL(), &
      & kh_tile(:,:,:) => NULL(), &
      & km_tile(:,:,:) => NULL(), &
      & kh(:,:) => NULL(), &
      & km(:,:) => NULL(), &
      & kh_neutral_tile(:,:,:) => NULL(), &
      & km_neutral_tile(:,:,:) => NULL(), &
      & kh_neutral(:,:) => NULL(), &
      & km_neutral(:,:) => NULL(), &
      & moist_rich_tile(:,:,:) => NULL(), &
      & albvisdir_tile (:,:,:) => NULL(),  & !< surface albedo over tiles for visible range, direct
      & albvisdif_tile (:,:,:) => NULL(),  & !< surface albedo over tiles for visible range, diffuse
      & albnirdir_tile (:,:,:) => NULL(),  & !< surface albedo over tiles for near-IR range, direct
      & albnirdif_tile (:,:,:) => NULL(),  & !< surface albedo over tiles for near-IR range, diffuse
      & albedo_tile    (:,:,:) => NULL(),  & !< surface albedo over tiles
      & albvisdir      (:,:  ) => NULL(),  & !< surface albedo for visible range, direct, grid-box mean
      & albvisdif      (:,:  ) => NULL(),  & !< surface albedo for visible range, diffuse, grid-box mean
      & albnirdir      (:,:  ) => NULL(),  & !< surface albedo for near-IR range, direct, grid-box mean
      & albnirdif      (:,:  ) => NULL(),  & !< surface albedo for near-IR range, diffuse, grid-box mean
      & albedo         (:,:  ) => NULL(),  & !< surface albedo, grid-box mean
      !
      ! Sea ice
      & q_ice_top      (:,:)   => NULL(),  & !< sea ice: energy flux available for surface melting [W m-2]
      & q_ice_bot      (:,:)   => NULL(),  & !< sea ice: energy flux at ice-ocean interface        [W m-2]
      & snow_thickness (:,:)   => NULL(),  & !< sea ice: snow thickness                            [m]
      !
      & t2m(:,:)            => NULL(), &
      & t2m_tile(:,:,:)     => NULL(), &
      & wind10m(:,:)        => NULL(), &
      & u10m(:,:)           => NULL(), &
      & v10m(:,:)           => NULL(), &
      & wind10m_tile(:,:,:) => NULL(), &
      & u10m_tile(:,:,:)    => NULL(), &
      & v10m_tile(:,:,:)    => NULL()

      ! & csat(:,:) => NULL(), &
      ! & cair(:,:) => NULL()
    INTEGER, POINTER :: &
      & nvalid(:,:) => NULL(), &
      & indices(:,:,:) => NULL()

  CONTAINS
    ! PROCEDURE :: Init => init_t_vdf_sfc_variable_set
    PROCEDURE :: Set_pointers => Set_pointers_diagnostics
  END TYPE t_vdf_sfc_diagnostics

  LOGICAL, SAVE :: l_init = .TRUE.

  CHARACTER(len=*), PARAMETER :: modname = 'mo_vdf_sfc'

CONTAINS

  FUNCTION t_vdf_sfc_construct(name, dt, domain) RESULT(result)

    CHARACTER(len=*), INTENT(in) :: name
    REAL(wp),         INTENT(in) :: dt
    TYPE(t_domain),      POINTER :: domain
    TYPE(t_vdf_sfc),     POINTER :: result

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_vdf_sfc_construct'

    ! CALL message(routine, '')

    ALLOCATE(t_vdf_sfc::result)
    !$ACC ENTER DATA COPYIN(result)
    ! Call Init of abstract parent class
    CALL result%Init(dt=dt, name=name, domain=domain)

    ! Initialize variable sets
    ALLOCATE(t_vdf_sfc_config :: result%config)
    result%config%list = build_sfc_config_list(result%domain)
    !$ACC ENTER DATA COPYIN(result%config)
    ! CALL result%config%Init(build_sfc_config_list(result%domain))

    ALLOCATE(t_vdf_sfc_inputs :: result%inputs)
    result%inputs%list = build_sfc_input_list(result%domain)
    !$ACC ENTER DATA COPYIN(result%inputs)
    ! CALL result%inputs%Init(build_sfc_input_list(result%domain))

    ALLOCATE(t_vdf_sfc_diagnostics :: result%diagnostics)
    result%diagnostics%list = build_sfc_diagnostic_list(result%domain)
    !$ACC ENTER DATA COPYIN(result%diagnostics)
    ! CALL result%diagnostics%list%allocator()
    ! CALL result%diagnostics%Set_pointers()

  END FUNCTION t_vdf_sfc_construct

  SUBROUTINE Compute(this, datetime)

    USE mo_tmx_surface_interface, ONLY: &
      & update_land, update_sea_ice, compute_lw_rad_net, compute_sw_rad_net, compute_albedo, &
      & compute_sfc_fluxes, compute_sfc_sat_spec_humidity
    ! USE mo_vdf_diag_smag,  ONLY: compute_sfc_fluxes, compute_sfc_sat_spec_humidity
    USE mo_physical_constants, ONLY: albedoW ! TODO
    USE mo_sea_ice_nml, ONLY: albi           ! TODO

    CLASS(t_vdf_sfc), INTENT(inout), TARGET :: this
    TYPE(t_datetime), OPTIONAL, INTENT(in), POINTER :: datetime     !< date and time at beginning of time step

    TYPE(t_vdf_sfc_config),      POINTER :: conf
    TYPE(t_vdf_sfc_inputs),      POINTER :: ins
    TYPE(t_vdf_sfc_diagnostics), POINTER :: diags
    CLASS(t_variable_set),       POINTER :: set

    INTEGER :: jg, jtile, isfc, jc, jb
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & old_tsfc, tend_tsfc, new_tsfc, &
      & new_qsfc
    REAL(wp) :: &
      & new_tsfc_rad(this%domain%nproma,this%domain%nblks_c,this%domain%ntiles), &
      & new_tsfc_eff(this%domain%nproma,this%domain%nblks_c,this%domain%ntiles), &
      & lwfl_net    (this%domain%nproma,this%domain%nblks_c), &
      & swfl_net    (this%domain%nproma,this%domain%nblks_c)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Compute'

    jg = 1

    !$ACC DATA CREATE(new_tsfc_rad, new_tsfc_eff, lwfl_net, swfl_net)

    SELECT TYPE (set => this%config)
    TYPE IS (t_vdf_sfc_config)
      conf => set
    END SELECT
    
    SELECT TYPE (set => this%inputs)
    TYPE IS (t_vdf_sfc_inputs)
      ins => set
    END SELECT
    
    SELECT TYPE (set => this%diagnostics)
    TYPE IS (t_vdf_sfc_diagnostics)
      diags => set
    END SELECT
    
    ! set => this%config
    ! SELECT TYPE (set)
    ! TYPE IS (t_vdf_sfc_config)
    !   conf => set
    !   
    ! END SELECT
    ! set => this%inputs
    ! SELECT TYPE (set)
    ! TYPE IS (t_vdf_sfc_inputs)
    !   ins => set
    !   
    ! END SELECT
    ! set => this%diagnostics
    ! SELECT TYPE (set)
    ! TYPE IS (t_vdf_sfc_diagnostics)
    !   diags => set
    !   
    ! END SELECT

    old_tsfc  => this%states    %Get_ptr_r3d('surface temperature')
    tend_tsfc => this%tendencies%Get_ptr_r3d('surface temperature')
    new_tsfc  => this%new_states%Get_ptr_r3d('surface temperature')
    new_qsfc  => this%new_states%Get_ptr_r3d('saturation specific humidity')

    ASSOCIATE( &
      & dtime    => conf%dtime,     &
      & domain   => this%domain,    &
      & tsfc_rad => diags%tsfc_rad, &
      & lwfl_up  => diags%lwfl_up,  &
      & swfl_up  => diags%swfl_up,  &
      & rlds     => ins%rlds,       &
      & rsds     => ins%rsds        &
      & )

      ! DO isfc=1,SIZE(this%domain%sfc_types)
    DO jtile=1,this%domain%ntiles
      isfc = this%domain%sfc_types(jtile)

      SELECT CASE(isfc)
      CASE(isfc_oce)
        ! Ocean surface temperature is calculated outside of this, set tendency to zero
!$OMP PARALLEL
        CALL init(tend_tsfc(:,:,jtile))
        CALL copy(old_tsfc(:,:,jtile), new_tsfc(:,:,jtile))
        CALL copy(old_tsfc(:,:,jtile), new_tsfc_rad(:,:,jtile))
        CALL copy(old_tsfc(:,:,jtile), new_tsfc_eff(:,:,jtile))
!$OMP END PARALLEL
        CALL compute_sfc_sat_spec_humidity(.FALSE., this%domain, this%domain%sfc_types(jtile), &
          & diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
          & ins%psfc(:,:), new_tsfc(:,:,jtile), diags%qsat_tile(:,:,jtile))
        ! TODO: This should be replaced by routine mo_surface_ocean:update_albedo_ocean from ECHAM6.2
!$OMP PARALLEL
        CALL init(diags%albvisdir_tile(:,:,jtile), albedoW)
        CALL init(diags%albvisdif_tile(:,:,jtile), albedoW)
        CALL init(diags%albnirdir_tile(:,:,jtile), albedoW)
        CALL init(diags%albnirdif_tile(:,:,jtile), albedoW)
!$OMP END PARALLEL
      CASE(isfc_ice)
        IF (conf%nice_thickness_classes /= 1) CALL finish(routine, 'Only one ice thickness class (kice) implemented!')

        CALL update_sea_ice(this%domain, this%dt, conf%cpd, &
          & old_tsfc(:,:,jtile), &
          & diags%lwfl_net_tile(:,:,jtile), diags%swfl_net_tile(:,:,jtile), &
          & diags%lhfl_tile(:,:,jtile), diags%shfl_tile(:,:,jtile), &
          & ins%ssfl(:,:), ins%ice_thickness(:,:), &
          & ins%emissivity(:,:), &
          ! inout &
          & diags%snow_thickness(:,:), & 
          ! out &
          & new_tsfc(:,:,jtile), &
          & diags%q_ice_top(:,:), diags%q_ice_bot(:,:), &
          & diags%albvisdir_tile(:,:,jtile), diags%albvisdif_tile(:,:,jtile), &
          & diags%albnirdir_tile(:,:,jtile), diags%albnirdif_tile(:,:,jtile)  &
          & )

!$OMP PARALLEL DO PRIVATE(jc, jb) ICON_OMP_DEFAULT_SCHEDULE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jb = 1, this%domain%nblks_c
            DO jc = 1, this%domain%nproma  
              new_tsfc_rad(jc,jb,jtile) = new_tsfc(jc,jb,jtile)
              new_tsfc_eff(jc,jb,jtile) = new_tsfc(jc,jb,jtile)
              tend_tsfc(jc,jb,jtile) = (new_tsfc(jc,jb,jtile) - old_tsfc(jc,jb,jtile)) / dtime
            END DO
          END DO
          !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO

        CALL compute_sfc_sat_spec_humidity(.FALSE., this%domain, this%domain%sfc_types(jtile), &
          & diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
          & ins%psfc(:,:), new_tsfc(:,:,jtile), diags%qsat_tile(:,:,jtile))
      CASE(isfc_lnd)
        CALL update_land(jg, this%domain, datetime, this%dt, conf%cpd, &
          & ins%dz(:,:), ins%psfc(:,:), ins%ta(:,:), ins%qa(:,:), ins%pa(:,:), &
          & ins%rsfl(:,:), ins%ssfl(:,:), &
          & ins%rlds(:,:), &
          & ins%rvds_dir(:,:), ins%rnds_dir(:,:), ins%rpds_dir(:,:), &
          & ins%rvds_dif(:,:), ins%rnds_dif(:,:), ins%rpds_dif(:,:), &
          & ins%cosmu0(:,:), &
          & diags%wind(:,:), diags%wind10m_tile(:,:,jtile), diags%rho_tile(:,:,jtile), ins%co2(:,:), &
          ! out
          & new_tsfc(:,:,jtile), new_tsfc_rad(:,:,jtile), new_tsfc_eff(:,:,jtile), diags%q_snocpymlt_lnd(:,:), &
          & new_qsfc(:,:,jtile), &
          & diags%albvisdir_tile(:,:,jtile), diags%albvisdif_tile(:,:,jtile), &
          & diags%albnirdir_tile(:,:,jtile), diags%albnirdif_tile(:,:,jtile), &
          & diags%kh_tile(:,:,jtile), diags%km_tile(:,:,jtile), &
          & diags%kh_neutral_tile(:,:,jtile), diags%km_neutral_tile(:,:,jtile) &
          & )
!$OMP PARALLEL DO PRIVATE(jc, jb) ICON_OMP_DEFAULT_SCHEDULE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jb = 1, this%domain%nblks_c
          DO jc = 1, this%domain%nproma  
            tend_tsfc(jc,jb,jtile) = (new_tsfc(jc,jb,jtile) - old_tsfc(jc,jb,jtile)) / dtime
          END DO
        END DO
        !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO

      END SELECT

      ! Compute surface fluxes for heat, water vapor and momentum from new state
      ! Note: for land, latent and sensible heat fluxes are directly taken from land model
      ! because they are part of the surface energy balance equation that is solved there
      CALL compute_sfc_fluxes( &
        ! Input
        & this%domain, this%domain%sfc_types(jtile), diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
        & conf%cpd, conf%cvd, &
        & ins%ua(:,:), ins%va(:,:), &
        ! & diags%theta_atm(:,:), ins%qa(:,:), diags%wind(:,:), diags%rho_tile(:,:,jtile),  &
        & ins%ta(:,:), ins%qa(:,:), diags%wind(:,:), ins%u_oce_current(:,:), ins%v_oce_current(:,:), &
        & diags%rho_tile(:,:,jtile),  &
        ! & diags%qsat_tile(:,:,jtile), diags%theta_tile(:,:,jtile), &
        & diags%qsat_tile(:,:,jtile), new_tsfc(:,:,jtile), &
        & diags%kh_tile(:,:,jtile), diags%km_tile(:,:,jtile), &
        ! Output
        & diags%evapotrans_tile(:,:,jtile), diags%lhfl_tile(:,:,jtile), diags%shfl_tile(:,:,jtile), &
        & diags%ustress_tile(:,:,jtile), diags%vstress_tile(:,:,jtile))

      CALL compute_lw_rad_net(this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
        & ins%emissivity(:,:), ins%rlds(:,:), new_tsfc_eff(:,:,jtile), diags%lwfl_net_tile(:,:,jtile) )

      CALL compute_sw_rad_net(this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
        & ins%rvds_dir(:,:), ins%rvds_dif(:,:), ins%rnds_dir(:,:), ins%rnds_dif(:,:), &
        & diags%albvisdir_tile(:,:,jtile), diags%albvisdif_tile(:,:,jtile), &
        & diags%albnirdir_tile(:,:,jtile), diags%albnirdif_tile(:,:,jtile), &
        & diags%swfl_net_tile(:,:,jtile))

        CALL compute_albedo(this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
        & ins%rsds(:,:), ins%rvds_dir(:,:), ins%rvds_dif(:,:), ins%rnds_dir(:,:), ins%rnds_dif(:,:), &
        & diags%albvisdir_tile(:,:,jtile), diags%albvisdif_tile(:,:,jtile), &
        & diags%albnirdir_tile(:,:,jtile), diags%albnirdif_tile(:,:,jtile), &
        & diags%albedo_tile(:,:,jtile))

    END DO

    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, new_tsfc(:,:,:), diags%tsfc, 'tsfc')

    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, new_tsfc_rad(:,:,:)**4, tsfc_rad, 'tsfc_rad')
!$OMP PARALLEL DO PRIVATE(jc, jb) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
        tsfc_rad(jc,jb) = tsfc_rad(jc,jb)**0.25_wp
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

    ! Aggregate surface fluxes
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices,diags%evapotrans_tile,diags%evapotrans,'evapotrans')
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%lhfl_tile,       diags%lhfl, 'lhfl')
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%shfl_tile,       diags%shfl, 'shfl')
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%ustress_tile,    diags%ustress, 'ustress')
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%vstress_tile,    diags%vstress, 'vstress')
    !
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%albvisdir_tile, diags%albvisdir, 'albvisdir')
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%albvisdif_tile, diags%albvisdif, 'albvisdif')
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%albnirdir_tile, diags%albnirdir, 'albnirdir')
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%albnirdif_tile, diags%albnirdif, 'albnirdif')
    !
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%albedo_tile,    diags%albedo, 'albedo')
    !
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%lwfl_net_tile, lwfl_net, 'lwfl_net')
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%swfl_net_tile, swfl_net, 'swfl_net')
!$OMP PARALLEL DO PRIVATE(jc, jb) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
        lwfl_up(jc,jb) = rlds(jc,jb) - lwfl_net(jc,jb)
        swfl_up(jc,jb) = rsds(jc,jb) - swfl_net(jc,jb)
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

    END ASSOCIATE

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE Compute

  SUBROUTINE Compute_diagnostics(this, datetime)

    USE mo_vdf_diag_smag, ONLY: compute_wind_speed, compute_atm_potential_temperature, compute_sfc_density, &
      ! & compute_sfc_sat_spec_humidity, compute_sfc_fluxes, compute_sfc_roughness, &
      & compute_sfc_potential_temperature, compute_sfc_exchange_coefficients, compute_moist_richardson
    USE mo_tmx_surface_interface, ONLY: &
      & compute_lw_rad_net, compute_sw_rad_net, compute_sfc_fluxes, compute_sfc_sat_spec_humidity, compute_sfc_roughness, &
      & update_land, compute_10m_wind

    CLASS(t_vdf_sfc), INTENT(inout), TARGET :: this
    TYPE(t_datetime), OPTIONAL, INTENT(in), POINTER :: datetime

    TYPE(t_vdf_sfc_config),      POINTER :: conf
    TYPE(t_vdf_sfc_inputs),      POINTER :: ins
    TYPE(t_vdf_sfc_diagnostics), POINTER :: diags

    INTEGER :: jg, jtile

    REAL(wp) :: rough_min
    REAL(wp), POINTER, DIMENSION(:,:,:) :: old_tsfc, old_qsat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Compute_diagnostics'

    ! CALL message(routine, '')

    jg = 1

    SELECT TYPE (set => this%config)
    TYPE IS (t_vdf_sfc_config)
      conf => set
    END SELECT
    
    SELECT TYPE (set => this%inputs)
    TYPE IS (t_vdf_sfc_inputs)
      ins => set
    END SELECT
    
    SELECT TYPE (set => this%diagnostics)
    TYPE IS (t_vdf_sfc_diagnostics)
      diags => set
      
    END SELECT

    old_tsfc => this%states%Get_ptr_r3d('surface temperature')
    old_qsat => this%states%Get_ptr_r3d('saturation specific humidity')

    ASSOCIATE( &
      & dtime    => conf%dtime,       &
      & domain   => this%domain,      &
      & fract_tile => ins%fract_tile, &
      & nvalid   => diags%nvalid,     &
      & indices  => diags%indices     &
      & )

    CALL compute_valid_indices(domain, fract_tile, nvalid, indices)

    CALL compute_wind_speed(this%domain, conf%min_sfc_wind, ins%ua(:,:), ins%va(:,:), diags%wind(:,:))

    CALL compute_atm_potential_temperature(this%domain, ins%ta(:,:), ins%tv(:,:), ins%pa(:,:), &
      & diags%theta_atm(:,:), diags%thetav_atm(:,:))

    DO jtile=1,this%domain%ntiles

      ! Surface saturated humidity
      CALL compute_sfc_sat_spec_humidity(l_init .AND. .NOT. isrestart(), this%domain, this%domain%sfc_types(jtile), &
        & diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
        & ins%psfc(:,:), ins%tsfc_tile(:,:,jtile), diags%qsat_tile(:,:,jtile))

      ! Density at surface
      CALL compute_sfc_density(this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
        ! & diags%qsat_tile(:,:,jtile), ins%psfc(:,:), ins%ta(:,:), diags%rho_tile(:,:,jtile))
        & diags%qsat_tile(:,:,jtile), ins%psfc(:,:), ins%tsfc_tile(:,:,jtile), diags%rho_tile(:,:,jtile))

      IF (l_init .AND. .NOT. isrestart()) THEN

        IF (this%domain%sfc_types(jtile) == isfc_lnd) THEN
          ! Compute inital 10m wind for update_land
          CALL compute_sfc_potential_temperature(this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
            & ins%psfc(:,:), ins%tsfc_tile(:,:,jtile), diags%qsat_tile(:,:,jtile), &
            & diags%theta_tile(:,:,jtile), diags%thetav_tile(:,:,jtile))
          CALL compute_moist_richardson(this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
            & conf%fsl, ins%zf(:,:), diags%thetav_atm(:,:), diags%thetav_tile(:,:,jtile), diags%wind(:,:), &
            & diags%moist_rich_tile(:,:,jtile))
          CALL compute_sfc_roughness(.TRUE., this%domain, this%domain%sfc_types(jtile), &
            & diags%nvalid(:,jtile), diags%indices(:,:,jtile), conf%min_rough, conf%rough_m_oce, conf%rough_m_ice, &
            & diags%wind(:,:), diags%km_tile(:,:,jtile), diags%rough_h_tile(:,:,jtile), diags%rough_m_tile(:,:,jtile))
          CALL compute_sfc_exchange_coefficients( &
            ! Input
            & this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), ins%dz(:,:),         &
            & ins%qa(:,:),     &
            & diags%theta_atm(:,:), diags%wind(:,:), diags%rough_m_tile(:,:,jtile), &
            & diags%theta_tile(:,:,jtile), diags%qsat_tile(:,:,jtile), &
            ! Output
            & diags%km_tile(:,:,jtile), diags%kh_tile(:,:,jtile), &
            & diags%km_neutral_tile(:,:,jtile), diags%kh_neutral_tile(:,:,jtile))
          CALL compute_10m_wind( &
            & this%domain, this%domain%sfc_types(jtile), diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
            & ins%zf(:,:), ins%zh(:,:), &
            & ins%ua(:,:), ins%va(:,:), ins%u_oce_current(:,:), ins%v_oce_current(:,:), &
            & diags%moist_rich_tile(:,:,jtile), diags%km_tile(:,:,jtile), diags%km_neutral_tile(:,:,jtile), &
            & diags%u10m_tile(:,:,jtile), diags%v10m_tile(:,:,jtile), diags%wind10m_tile(:,:,jtile) &
            & )
          
          ! Call land in quasi-diagnostic mode with a time step of 1 second
          CALL update_land(jg, this%domain, datetime, 1._wp, conf%cpd, &
            & ins%dz(:,:), ins%psfc(:,:), ins%ta(:,:), ins%qa(:,:), ins%pa(:,:), &
            & ins%rsfl(:,:), ins%ssfl(:,:), &
            & ins%rlds(:,:), &
            & ins%rvds_dir(:,:), ins%rnds_dir(:,:), ins%rpds_dir(:,:), &
            & ins%rvds_dif(:,:), ins%rnds_dif(:,:), ins%rpds_dif(:,:), &
            & ins%cosmu0(:,:), &
            & diags%wind(:,:), diags%wind10m_tile(:,:,jtile), diags%rho_tile(:,:,jtile), ins%co2(:,:), &
            ! out
            & tsfc=ins%tsfc_tile(:,:,jtile), & !, new_tsfc_rad(:,:,jtile), new_tsfc_eff(:,:,jtile), diags%q_snocpymlt_lnd(:,:), &
            & qsat=diags%qsat_tile(:,:,jtile) &
            ! & diags%albvisdir_tile(:,:,jtile), diags%albvisdif_tile(:,:,jtile), &
            ! & diags%albnirdir_tile(:,:,jtile), diags%albnirdif_tile(:,:,jtile), &
            ! & diags%kh_tile(:,:,jtile), diags%km_tile(:,:,jtile), &
            ! & diags%kh_neutral_tile(:,:,jtile), diags%km_neutral_tile(:,:,jtile) &
            & )
        END IF
      END IF
  
      ! Surface potential temperature
      CALL compute_sfc_potential_temperature(this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
        & ins%psfc(:,:), ins%tsfc_tile(:,:,jtile), diags%qsat_tile(:,:,jtile), &
        & diags%theta_tile(:,:,jtile), diags%thetav_tile(:,:,jtile))

      ! Moist Richardson number
      CALL compute_moist_richardson(this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
        & conf%fsl, ins%zf(:,:), diags%thetav_atm(:,:), diags%thetav_tile(:,:,jtile), diags%wind(:,:), &
        & diags%moist_rich_tile(:,:,jtile))
    
      ! Surface roughness length
      IF (this%domain%sfc_types(jtile) == isfc_oce) THEN
        rough_min = conf%min_rough
      ELSE
        rough_min = 0._wp
      END IF
      ! Uses old value of km_tile before computation of new exchange coefficients (only in case of ocean)
      CALL compute_sfc_roughness(l_init .AND. .NOT. isrestart(), this%domain, this%domain%sfc_types(jtile), &
        & diags%nvalid(:,jtile), diags%indices(:,:,jtile), rough_min, conf%rough_m_oce, conf%rough_m_ice, &
        & diags%wind(:,:), diags%km_tile(:,:,jtile), diags%rough_h_tile(:,:,jtile), diags%rough_m_tile(:,:,jtile))

      ! Surface exchange coefficients
      ! Note: coefficients for land will be computed in the land model itself
      IF (this%domain%sfc_types(jtile) /= isfc_lnd) THEN
        CALL compute_sfc_exchange_coefficients( &
            ! Input
          & this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), ins%dz(:,:),         &
          & ins%qa(:,:),     &
          & diags%theta_atm(:,:), diags%wind(:,:), diags%rough_m_tile(:,:,jtile), &
          & diags%theta_tile(:,:,jtile), diags%qsat_tile(:,:,jtile), &
            ! Output
          & diags%km_tile(:,:,jtile), diags%kh_tile(:,:,jtile), &
          & diags%km_neutral_tile(:,:,jtile), diags%kh_neutral_tile(:,:,jtile))
      END IF

      IF (this%domain%sfc_types(jtile) == isfc_ice) THEN
        ! The ice model needs the latent and sensible heat fluxes as input (based on old state)
        CALL compute_sfc_fluxes( &
          ! Input
          & this%domain, isfc_ice, diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
          & conf%cpd, conf%cvd, &
          & ins%ua(:,:), ins%va(:,:), &
          ! & diags%theta_atm(:,:), ins%qa(:,:), diags%wind(:,:), diags%rho_tile(:,:,jtile),  &
          & ins%ta(:,:), ins%qa(:,:), diags%wind(:,:), ins%u_oce_current(:,:), ins%v_oce_current(:,:), &
          & diags%rho_tile(:,:,jtile),  &
          ! & diags%qsat_tile(:,:,jtile), diags%theta_tile(:,:,jtile), &
          & diags%qsat_tile(:,:,jtile), ins%tsfc_tile(:,:,jtile), &
          & diags%kh_tile(:,:,jtile), diags%km_tile(:,:,jtile), &
          ! Output
          & diags%evapotrans_tile(:,:,jtile), diags%lhfl_tile(:,:,jtile), diags%shfl_tile(:,:,jtile), &
          & diags%ustress_tile(:,:,jtile), diags%vstress_tile(:,:,jtile))

          ! The ice model needs the shortwave net surface flux as input (based on old state)
        CALL compute_sw_rad_net(this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
          & ins%rvds_dir(:,:), ins%rvds_dif(:,:), ins%rnds_dir(:,:), ins%rnds_dif(:,:), &
          & diags%albvisdir_tile(:,:,jtile), diags%albvisdif_tile(:,:,jtile), &
          & diags%albnirdir_tile(:,:,jtile), diags%albnirdif_tile(:,:,jtile), &
          ! Out
          & diags%swfl_net_tile(:,:,jtile))

          ! The ice model needs the longwave net surface flux as input (based on old state)
        CALL compute_lw_rad_net(this%domain, diags%nvalid(:,jtile), diags%indices(:,:,jtile), &
          & ins%emissivity(:,:), ins%rlds(:,:), ins%tsfc_tile(:,:,jtile), diags%lwfl_net_tile(:,:,jtile) )
      END IF

    END DO

    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%rough_m_tile, diags%rough_m, 'rough_m')
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%rough_h_tile, diags%rough_h, 'rough_h')

    END ASSOCIATE

  END SUBROUTINE Compute_diagnostics

  SUBROUTINE Update_diagnostics(this)

    USE mo_tmx_surface_interface, ONLY: compute_2m_temperature, compute_10m_wind

    CLASS(t_vdf_sfc), INTENT(inout), TARGET :: this

    TYPE(t_vdf_sfc_config),      POINTER :: conf
    TYPE(t_vdf_sfc_inputs),      POINTER :: ins
    TYPE(t_vdf_sfc_diagnostics), POINTER :: diags

    INTEGER :: jg, jtile
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & new_ta, new_qv, new_qc, new_qi, new_ua, new_va, &
      & new_tsf
    REAL(wp) :: new_tsfc(this%domain%nproma,this%domain%nblks_c,this%domain%ntiles)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Update_diagnostics'

    jg = 1

    SELECT TYPE (set => this%config)
    TYPE IS (t_vdf_sfc_config)
      conf => set
    END SELECT
    
    SELECT TYPE (set => this%inputs)
    TYPE IS (t_vdf_sfc_inputs)
      ins => set
    END SELECT
    
    SELECT TYPE (set => this%diagnostics)
    TYPE IS (t_vdf_sfc_diagnostics)
      diags => set
    END SELECT
    
    
    ! CALL message(routine, 'start')

    ! state => this%states%Get_ptr_r3d('surface temperature')
    ! tend => this%tendencies%Get_ptr_r3d('surface temperature')

    ! new_tsfc => this%atmo%new_states%Get_ptr_r3d('surface temperature')
    ! new_qv   => this%atmo%new_states%Get_ptr_r3d('water vapor')
    ! new_qc   => this%atmo%new_states%Get_ptr_r3d('cloud water')
    ! new_qi   => this%atmo%new_states%Get_ptr_r3d('cloud ice')
    ! new_ua   => this%atmo%new_states%Get_ptr_r3d('zonal wind')
    ! new_va   => this%atmo%new_states%Get_ptr_r3d('meridional wind')

    ! DO jtile=1,this%domain%ntiles

    ! END DO

    ! TODO: Update roughness length for heat and momentum

    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%km_tile, diags%km)
    CALL average_tiles(this%domain, ins%fract_tile, diags%nvalid, diags%indices, diags%kh_tile, diags%kh)

    IF (l_init) l_init = .FALSE.

    ! CALL message(routine, 'end')

  END SUBROUTINE Update_diagnostics

  SUBROUTINE compute_valid_indices(domain, fract_tile, nvalid, indices)

    USE mo_index_list, ONLY: generate_index_list_batched, generate_index_list

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp),       INTENT(in)          :: fract_tile(:,:,:)
    INTEGER,        INTENT(out)         :: nvalid(:,:), indices(:,:,:)

    INTEGER :: ib, ic, ics, ice, jsfc, ntiles
    INTEGER(i1) :: pfrc_test(domain%nproma,domain%nblks_c, domain%ntiles)

    ntiles = domain%ntiles

    !$ACC DATA CREATE(pfrc_test) PRESENT(indices, nvalid)

!$OMP PARALLEL

    ! CALL init(nvalid)
    ! CALL init(indices)

!$OMP DO PRIVATE(ib, ic, ics, ice, jsfc) ICON_OMP_RUNTIME_SCHEDULE
    DO ib = domain%i_startblk_c, domain%i_endblk_c
      ics = domain%i_startidx_c(ib)
      ice = domain%i_endidx_c  (ib)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jsfc = 1, ntiles
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = ics, ice
          pfrc_test(ic,ib,jsfc) = MERGE(1_i1, 0_i1, fract_tile(ic,ib,jsfc) > 0.0_wp)
          indices(ic,ib,jsfc) = 0
        END DO
        !$ACC END LOOP
        nvalid(ib,jsfc) = 0
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL
    END DO
!$OMP END DO

!$OMP END PARALLEL

    ! TODO: for some reason, I can't get generate_index_list_batched to work on GPUs

    !$ACC WAIT(1)
    !$ACC UPDATE HOST(pfrc_test, indices, nvalid)

!$OMP PARALLEL DO PRIVATE(ib, ics, ice) ICON_OMP_RUNTIME_SCHEDULE
    DO ib = domain%i_startblk_c, domain%i_endblk_c
      ics = domain%i_startidx_c(ib)
      ice = domain%i_endidx_c  (ib)

      CALL generate_index_list_batched( &
        & pfrc_test(:,ib,:), indices(:,ib,:), ics, ice, nvalid(ib,:), opt_use_acc=.FALSE.)
  
    END DO
!$OMP END PARALLEL DO

    !$ACC UPDATE DEVICE(indices, nvalid)
    !$noACC UPDATE HOST(indices, nvalid)

    !$noACC WAIT(1)

    !$ACC END DATA

  END SUBROUTINE compute_valid_indices

  SUBROUTINE average_tiles(domain, fract_tile, nvalid, indices, var_in, var_out, msg)

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), INTENT(in)  :: &
      & fract_tile(:,:,:)
    REAL(wp), TARGET, INTENT(in) :: &
      & var_in(:,:,:)
    INTEGER,  INTENT(in)  :: &
      & nvalid(:,:),         &
      & indices(:,:,:)
    REAL(wp), INTENT(out) :: &
      & var_out(:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: msg

    INTEGER :: jb, jls, js, jsfc
    REAL(wp), POINTER :: ptr3d(:,:,:)
    LOGICAL :: not_is_present

    CHARACTER(len=*), PARAMETER :: routine = modname//':average_tiles'

    ! IF (PRESENT(msg)) CALL message(routine, 'working on '//msg)

    IF (domain%ntiles < 1) CALL finish(routine, 'This should not happen - ntiles < 1')

    not_is_present = .TRUE.
    ptr3d => var_in
    !$ACC DATA COPYIN(ptr3d) IF(not_is_present)

    IF (domain%ntiles == 1) THEN

!$OMP PARALLEL
      CALL copy(ptr3d(:,:,1), var_out(:,:))
!$OMP END PARALLEL

    ELSE

!$OMP PARALLEL
      CALL init(var_out)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb,jls,js,jsfc) ICON_OMP_RUNTIME_SCHEDULE
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jb = domain%i_startblk_c,domain%i_endblk_c
        !$ACC LOOP GANG(STATIC: 1)
        DO jsfc = 1, domain%ntiles
          !$ACC LOOP VECTOR
          DO jls = 1, nvalid(jb,jsfc)
            js=indices(jls,jb,jsfc)
            var_out(js,jb) = var_out(js,jb) + fract_tile(js,jb,jsfc) * ptr3d(js,jb,jsfc)
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
!$OMP END PARALLEL DO

    END IF

    !$ACC END DATA

    ! IF (PRESENT(msg)) CALL message(routine, 'average complete for '//msg)

  END SUBROUTINE average_tiles

  FUNCTION build_sfc_config_list(domain) RESULT(configlist)

    TYPE(t_variable_list) :: configlist
    TYPE(t_domain),        INTENT(in)    :: domain

    INTEGER :: nproma, nblks_c, nlev, nlevp1
    INTEGER :: shape_0d(1)

    configlist = t_variable_list('sfc_config')

    nproma  = domain%nproma
    nblks_c = domain%nblks_c
    nlev    = domain%nlev
    nlevp1  = nlev + 1

    shape_0d = [0]
    CALL configlist%append(t_variable('time step', shape_0d, "s", type_id="real"))
    CALL configlist%append(t_variable('cpd', shape_0d, "", type_id="real"))
    CALL configlist%append(t_variable('cvd', shape_0d, "", type_id="real"))
    CALL configlist%append(t_variable('minimum surface wind speed', shape_0d, "m/s", type_id="real"))
    CALL configlist%append(t_variable('ocean roughness length', shape_0d, "m", type_id="real"))
    CALL configlist%append(t_variable('ice roughness length', shape_0d, "m", type_id="real"))
    CALL configlist%append(t_variable('minimal roughness length', shape_0d, "m", type_id="real"))
    CALL configlist%append(t_variable('weight for interpolation to surface_layer mid level', shape_0d, "", type_id="real"))
    CALL configlist%append(t_variable('number of sea ice thickness classes', shape_0d, "", type_id="int"))

  END FUNCTION build_sfc_config_list

  SUBROUTINE Set_pointers_config(this)

    CLASS(t_vdf_sfc_config), INTENT(inout) :: this

    SELECT TYPE (this)
    TYPE IS (t_vdf_sfc_config)
      this%dtime => this%list%Get_ptr_r0d('time step')
      
      this%cpd => this%list%Get_ptr_r0d('cpd')
      
      this%cvd => this%list%Get_ptr_r0d('cvd')
      
      this%min_sfc_wind => this%list%Get_ptr_r0d('minimum surface wind speed')
      
      this%min_rough   => this%list%Get_ptr_r0d('minimal roughness length')
      
      this%rough_m_oce => this%list%Get_ptr_r0d('ocean roughness length')
      
      this%rough_m_ice => this%list%Get_ptr_r0d('ice roughness length')
      
      this%fsl         => this%list%Get_ptr_r0d('weight for interpolation to surface_layer mid level')
      
      this%nice_thickness_classes => this%list%Get_ptr_i0d('number of sea ice thickness classes')
      
    END SELECT

  END SUBROUTINE Set_pointers_config

  FUNCTION build_sfc_input_list(domain) RESULT(inlist)

    TYPE(t_variable_list) :: inlist
    TYPE(t_domain),        INTENT(in)    :: domain

    INTEGER :: nproma, nblks_c, ntiles
    INTEGER :: shape_2d(2), shape_3d(3)

    inlist = t_variable_list('sfc_inputs')

    nproma  = domain%nproma
    nblks_c = domain%nblks_c
    ntiles  = domain%ntiles

    shape_2d = [nproma,nblks_c]
    CALL inlist%append(t_variable('atm temperature', shape_2d, "K", type_id="real"))
    CALL inlist%append(t_variable('atm virtual temperature', shape_2d, "K", type_id="real"))
    CALL inlist%append(t_variable('atm zonal wind', shape_2d, "m/s", type_id="real"))
    CALL inlist%append(t_variable('atm meridional wind', shape_2d, "m/s", type_id="real"))
    CALL inlist%append(t_variable('atm total water', shape_2d, "kg/kg", type_id="real"))
    CALL inlist%append(t_variable('atm full level pressure', shape_2d, "Pa", type_id="real"))
    CALL inlist%append(t_variable('surface pressure', shape_2d, "Pa", type_id="real"))
    CALL inlist%append(t_variable('atm geometric height full', shape_2d, "Pa", type_id="real"))
    CALL inlist%append(t_variable('sfc geometric height half', shape_2d, "Pa", type_id="real"))
    CALL inlist%append(t_variable('reference height in surface layer times 2', shape_2d, "m", type_id="real"))
    CALL inlist%append(t_variable('surface rain flux, large-scale', shape_2d, "kg m-2 s-1", type_id="real"))
    CALL inlist%append(t_variable('surface snow flux, large-scale', shape_2d, "kg m-2 s-1", type_id="real"))

    CALL inlist%append(t_variable('surface downward longwave radiation',                shape_2d, "W m-2", type_id="real"))
    CALL inlist%append(t_variable('surface downward shortwave radiation',               shape_2d, "W m-2", type_id="real"))
    CALL inlist%append(t_variable('all-sky surface downward direct visible radiation',  shape_2d, "W m-2", type_id="real"))
    CALL inlist%append(t_variable('all-sky surface downward direct near-IR radiation',  shape_2d, "W m-2", type_id="real"))
    CALL inlist%append(t_variable('all-sky surface downward direct PAR radiation',      shape_2d, "W m-2", type_id="real"))  
    CALL inlist%append(t_variable('all-sky surface downward diffuse visible radiation', shape_2d, "W m-2", type_id="real"))
    CALL inlist%append(t_variable('all-sky surface downward diffuse near-IR radiation', shape_2d, "W m-2", type_id="real"))
    CALL inlist%append(t_variable('all-sky surface downward diffuse PAR radiation',     shape_2d, "W m-2", type_id="real"))

    CALL inlist%append(t_variable('longwave surface emissivity', shape_2d, "", type_id="real"))
    CALL inlist%append(t_variable('cosine of zenith angle', shape_2d, "", type_id="real"))
    CALL inlist%append(t_variable('atm CO2 concentration', shape_2d, "", type_id="real"))
    CALL inlist%append(t_variable(('u-component of ocean current'), shape_2d, "m s-1", type_id="real"))
    CALL inlist%append(t_variable(('v-component of ocean current'), shape_2d, "m s-1", type_id="real"))

    CALL inlist%append(t_variable('thickness of sea ice', shape_2d, "m", type_id="real"))

    shape_3d = [nproma,nblks_c,ntiles]
    CALL inlist%append(t_variable('surface temperature, tile', shape_3d, "K", type_id="real"))
    ! CALL inlist%append(t_variable('roughness length momentum', shape_3d, "m", type_id="real"))
    ! CALL inlist%append(t_variable('roughness length heat', shape_3d, "m", type_id="real"))
    CALL inlist%append(t_variable('grid box fraction of tiles', shape_3d, "", type_id="real"))

  END FUNCTION build_sfc_input_list

  SUBROUTINE Set_pointers_inputs(this)

    CLASS(t_vdf_sfc_inputs), INTENT(inout) :: this

    SELECT TYPE (this)
    TYPE IS (t_vdf_sfc_inputs)
      this%ta => this%list%Get_ptr_r2d('atm temperature')
      
      this%tv => this%list%Get_ptr_r2d('atm virtual temperature')
      
      this%ua => this%list%Get_ptr_r2d('atm zonal wind')
      
      this%va => this%list%Get_ptr_r2d('atm meridional wind')
      
      this%qa => this%list%Get_ptr_r2d('atm total water')
      
      this%pa => this%list%Get_ptr_r2d('atm full level pressure')
      
      this%psfc => this%list%Get_ptr_r2d('surface pressure')
      
      this%rsfl => this%list%Get_ptr_r2d('surface rain flux, large-scale')
      
      this%ssfl => this%list%Get_ptr_r2d('surface snow flux, large-scale')
      

      this%rlds     => this%list%Get_ptr_r2d('surface downward longwave radiation')
      
      this%rsds     => this%list%Get_ptr_r2d('surface downward shortwave radiation')
      
      this%rvds_dir => this%list%Get_ptr_r2d('all-sky surface downward direct visible radiation')
      
      this%rnds_dir => this%list%Get_ptr_r2d('all-sky surface downward direct near-IR radiation')
      
      this%rpds_dir => this%list%Get_ptr_r2d('all-sky surface downward direct PAR radiation') 
      
      this%rvds_dif => this%list%Get_ptr_r2d('all-sky surface downward diffuse visible radiation')
      
      this%rnds_dif => this%list%Get_ptr_r2d('all-sky surface downward diffuse near-IR radiation')
      
      this%rpds_dif => this%list%Get_ptr_r2d('all-sky surface downward diffuse PAR radiation')
      

      this%emissivity => this%list%Get_ptr_r2d('longwave surface emissivity')
      
      this%cosmu0     => this%list%Get_ptr_r2d('cosine of zenith angle')
      
      this%co2        => this%list%Get_ptr_r2d('atm CO2 concentration')
      

      this%tsfc_tile => this%list%Get_ptr_r3d('surface temperature, tile')
      
      ! this%rough_m => this%list%Get_ptr_r3d('roughness length momentum')
      ! this%rough_h => this%list%Get_ptr_r3d('roughness length heat')
      this%u_oce_current => this%list%Get_ptr_r2d('u-component of ocean current')
      
      this%v_oce_current => this%list%Get_ptr_r2d('v-component of ocean current')
      

      this%ice_thickness  => this%list%Get_ptr_r2d('thickness of sea ice')
      
  
      this%zf => this%list%Get_ptr_r2d('atm geometric height full')
      
      this%zh => this%list%Get_ptr_r2d('sfc geometric height half')
      
      this%dz => this%list%Get_ptr_r2d('reference height in surface layer times 2')
      
      this%fract_tile => this%list%get_ptr_r3d('grid box fraction of tiles')
      

    END SELECT

  END SUBROUTINE Set_pointers_inputs

  FUNCTION build_sfc_diagnostic_list(domain) RESULT(diaglist)

    TYPE(t_variable_list) :: diaglist
    TYPE(t_domain),        INTENT(in)    :: domain

    INTEGER :: nproma, nblks_c, ntiles
    INTEGER :: shape_2d(2), shape_3d(3)

    diaglist = t_variable_list('sfc_diagnostics')

    nproma  = domain%nproma
    nblks_c = domain%nblks_c
    ntiles  = domain%ntiles

    shape_3d = [nproma,nblks_c,ntiles]
    CALL diaglist%append(t_variable('roughness length heat, tile', shape_3d, "m", type_id="real"))
    CALL diaglist%append(t_variable('roughness length momentum, tile', shape_3d, "m", type_id="real"))
    !
    CALL diaglist%append(t_variable('sfc saturation specific humidity, tile', shape_3d, "kg/kg", type_id="real"))
    !
    CALL diaglist%append(t_variable('exchange coefficient for scalar, tile', shape_3d, "-", type_id="real"))
    CALL diaglist%append(t_variable('exchange coefficient for momentum, tile', shape_3d, "-", type_id="real"))
    CALL diaglist%append(t_variable('neutral exchange coefficient for scalar, tile', shape_3d, "-", type_id="real"))
    CALL diaglist%append(t_variable('neutral exchange coefficient for momentum, tile', shape_3d, "-", type_id="real"))
    CALL diaglist%append(t_variable('moist richardson number, tile', shape_3d, "-", type_id="real"))
    !
    CALL diaglist%append(t_variable('sfc density, tile', shape_3d, "kg/m3", type_id="real"))
    CALL diaglist%append(t_variable('sfc potential temperature, tile', shape_3d, "K", type_id="real"))
    CALL diaglist%append(t_variable('sfc virtual potential temperature, tile', shape_3d, "K", type_id="real"))
    !
    CALL diaglist%append(t_variable('sfc evapotranspiration, tile', shape_3d, "W m-2", type_id="real"))
    CALL diaglist%append(t_variable('sfc latent heat flux, tile', shape_3d, "W m-2", type_id="real"))
    CALL diaglist%append(t_variable('sfc sensible heat flux, tile', shape_3d, "W m-2", type_id="real"))
    CALL diaglist%append(t_variable('sfc zonal wind stress, tile', shape_3d, "N m-2", type_id="real"))
    CALL diaglist%append(t_variable('sfc mer. wind stress, tile', shape_3d, "N m-2", type_id="real"))
    !
    CALL diaglist%append(t_variable('sfc longwave net flux, tile', shape_3d, "W m-2", type_id="real"))
    CALL diaglist%append(t_variable('sfc shortwave net flux, tile', shape_3d, "W m-2", type_id="real"))

    CALL diaglist%append(t_variable('albedo VIS direct, tile', shape_3d, "", type_id="real"))
    CALL diaglist%append(t_variable('albedo VIS diffuse, tile', shape_3d, "", type_id="real"))
    CALL diaglist%append(t_variable('albedo NIR direct, tile', shape_3d, "", type_id="real"))
    CALL diaglist%append(t_variable('albedo NIR diffuse, tile', shape_3d, "", type_id="real"))
    CALL diaglist%append(t_variable('albedo, tile', shape_3d, "", type_id="real"))

    shape_2d = [nproma,nblks_c]
    CALL diaglist%append(t_variable('heating used to melt snow on canopy', shape_2d, "W m-2", type_id="real"))
    CALL diaglist%append(t_variable('atm wind speed', shape_2d, "m", type_id="real"))
    CALL diaglist%append(t_variable('atm potential temperature', shape_2d, "K", type_id="real"))
    CALL diaglist%append(t_variable('atm virtual potential temperature', shape_2d, "K", type_id="real"))
    !
    CALL diaglist%append(t_variable('roughness length heat', shape_2d, "m", type_id="real"))
    CALL diaglist%append(t_variable('roughness length momentum', shape_2d, "m", type_id="real"))
    !
    CALL diaglist%append(t_variable('exchange coefficient for scalar', shape_2d, "-", type_id="real"))
    CALL diaglist%append(t_variable('exchange coefficient for momentum', shape_2d, "-", type_id="real"))
    CALL diaglist%append(t_variable('neutral exchange coefficient for scalar', shape_2d, "-", type_id="real"))
    CALL diaglist%append(t_variable('neutral exchange coefficient for momentum', shape_2d, "-", type_id="real"))
    !
    CALL diaglist%append(t_variable('sfc evapotranspiration', shape_2d, "W m-2", type_id="real"))
    CALL diaglist%append(t_variable('sfc latent heat flux', shape_2d, "W m-2", type_id="real"))
    CALL diaglist%append(t_variable('sfc sensible heat flux', shape_2d, "W m-2", type_id="real"))
    CALL diaglist%append(t_variable('sfc zonal wind stress', shape_2d, "N m-2", type_id="real"))
    CALL diaglist%append(t_variable('sfc mer. wind stress', shape_2d, "N m-2", type_id="real"))
    !
    CALL diaglist%append(t_variable('sfc temperature', shape_2d, "K", type_id="real"))
    CALL diaglist%append(t_variable('sfc radiative temperature', shape_2d, "K", type_id="real"))
    CALL diaglist%append(t_variable('sfc longwave upward flux', shape_2d, "W m-2", type_id="real"))
    CALL diaglist%append(t_variable('sfc shortwave upward flux', shape_2d, "W m-2", type_id="real"))
    !
    CALL diaglist%append(t_variable('albedo VIS direct', shape_2d, "", type_id="real"))
    CALL diaglist%append(t_variable('albedo VIS diffuse', shape_2d, "", type_id="real"))
    CALL diaglist%append(t_variable('albedo NIR direct', shape_2d, "", type_id="real"))
    CALL diaglist%append(t_variable('albedo NIR diffuse', shape_2d, "", type_id="real"))
    CALL diaglist%append(t_variable('albedo', shape_2d, "", type_id="real"))
    !
    CALL diaglist%append(t_variable('energy flux available for surface melting of sea ice', shape_2d, "W m-2", type_id="real"))
    CALL diaglist%append(t_variable('energy flux at ice-ocean interface', shape_2d, "W m-2", type_id="real"))
    CALL diaglist%append(t_variable('thickness of snow on sea ice', shape_2d, "m", type_id="real"))
    !
    CALL diaglist%append(t_variable('2m temperature', shape_2d, "K", type_id="real"))
    CALL diaglist%append(t_variable('2m temperature, tile', shape_3d, "K", type_id="real"))
    CALL diaglist%append(t_variable('10m wind speed', shape_2d, "m/s", type_id="real"))
    CALL diaglist%append(t_variable('10m zonal wind', shape_2d, "m/s", type_id="real"))
    CALL diaglist%append(t_variable('10m meridional wind', shape_2d, "m/s", type_id="real"))
    CALL diaglist%append(t_variable('10m wind speed, tile', shape_3d, "m/s", type_id="real"))
    CALL diaglist%append(t_variable('10m zonal wind, tile', shape_3d, "m/s", type_id="real"))
    CALL diaglist%append(t_variable('10m meridional wind, tile', shape_3d, "m/s", type_id="real"))
    !
    shape_2d = [nblks_c,ntiles]
    CALL diaglist%append(t_variable('number of valid points per block, tile', shape_2d, "", type_id="int"))
    shape_3d = [nproma,nblks_c,ntiles]
    CALL diaglist%append(t_variable('indices of valid points on chunk per block, tile', shape_3d, "", type_id="int"))

  END FUNCTION build_sfc_diagnostic_list

  SUBROUTINE Set_pointers_diagnostics(this)

    CLASS(t_vdf_sfc_diagnostics), INTENT(inout) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_pointers_diagnostics'

    ! CALL message(routine,'')

    SELECT TYPE (this)
    TYPE IS (t_vdf_sfc_diagnostics)
      this%wind            => this%list%Get_ptr_r2d('atm wind speed')
      
      this%theta_atm       => this%list%Get_ptr_r2d('atm potential temperature')
      
      this%thetav_atm      => this%list%Get_ptr_r2d('atm virtual potential temperature')
      
      this%theta_tile      => this%list%Get_ptr_r3d('sfc potential temperature, tile')
      
      this%thetav_tile     => this%list%Get_ptr_r3d('sfc virtual potential temperature, tile')
      
      !
      this%rough_h_tile    => this%list%Get_ptr_r3d('roughness length heat, tile')
      
      this%rough_m_tile    => this%list%Get_ptr_r3d('roughness length momentum, tile')
      
      this%rough_h         => this%list%Get_ptr_r2d('roughness length heat')
      
      this%rough_m         => this%list%Get_ptr_r2d('roughness length momentum')
      
      !
      this%qsat_tile       => this%list%Get_ptr_r3d('sfc saturation specific humidity, tile')
      
      ! this%pcpt_tile        => this%list%Get_ptr_r3d('dry static energy')
      !
      this%kh_tile         => this%list%Get_ptr_r3d('exchange coefficient for scalar, tile')
      
      this%km_tile         => this%list%Get_ptr_r3d('exchange coefficient for momentum, tile')
      
      this%kh              => this%list%Get_ptr_r2d('exchange coefficient for scalar')
      
      this%km              => this%list%Get_ptr_r2d('exchange coefficient for momentum')
      
      this%kh_neutral_tile => this%list%Get_ptr_r3d('neutral exchange coefficient for scalar, tile')
      
      this%km_neutral_tile => this%list%Get_ptr_r3d('neutral exchange coefficient for momentum, tile')
      
      this%kh_neutral      => this%list%Get_ptr_r2d('neutral exchange coefficient for scalar')
      
      this%km_neutral      => this%list%Get_ptr_r2d('neutral exchange coefficient for momentum')
      
      this%moist_rich_tile => this%list%Get_ptr_r3d('moist richardson number, tile')
      
      !
      this%rho_tile        => this%list%Get_ptr_r3d('sfc density, tile')
      
      !
      this%evapotrans_tile => this%list%Get_ptr_r3d('sfc evapotranspiration, tile')
      
      this%lhfl_tile       => this%list%Get_ptr_r3d('sfc latent heat flux, tile')
      
      this%shfl_tile       => this%list%Get_ptr_r3d('sfc sensible heat flux, tile')
      
      this%q_snocpymlt_lnd => this%list%Get_ptr_r2d('heating used to melt snow on canopy')
      
      this%ustress_tile    => this%list%Get_ptr_r3d('sfc zonal wind stress, tile')
      
      this%vstress_tile    => this%list%Get_ptr_r3d('sfc mer. wind stress, tile')
      
      this%evapotrans      => this%list%Get_ptr_r2d('sfc evapotranspiration')
      
      this%lhfl            => this%list%Get_ptr_r2d('sfc latent heat flux')
      
      this%shfl            => this%list%Get_ptr_r2d('sfc sensible heat flux')
      
      this%ustress         => this%list%Get_ptr_r2d('sfc zonal wind stress')
      
      this%vstress         => this%list%Get_ptr_r2d('sfc mer. wind stress')
      
   
      this%tsfc            => this%list%Get_ptr_r2d('sfc temperature')
      
      this%tsfc_rad        => this%list%Get_ptr_r2d('sfc radiative temperature')
      
      this%lwfl_up         => this%list%Get_ptr_r2d('sfc longwave upward flux')
      
      this%swfl_up         => this%list%Get_ptr_r2d('sfc shortwave upward flux')
      
      this%lwfl_net_tile   => this%list%Get_ptr_r3d('sfc longwave net flux, tile')
      
      this%swfl_net_tile   => this%list%Get_ptr_r3d('sfc shortwave net flux, tile')
      

      this%albvisdir_tile  => this%list%Get_ptr_r3d('albedo VIS direct, tile')
      
      this%albvisdif_tile  => this%list%Get_ptr_r3d('albedo VIS diffuse, tile')
      
      this%albnirdir_tile  => this%list%Get_ptr_r3d('albedo NIR direct, tile')
      
      this%albnirdif_tile  => this%list%Get_ptr_r3d('albedo NIR diffuse, tile')
      
      this%albedo_tile     => this%list%Get_ptr_r3d('albedo, tile')
      
      this%albvisdir       => this%list%Get_ptr_r2d('albedo VIS direct')
      
      this%albvisdif       => this%list%Get_ptr_r2d('albedo VIS diffuse')
      
      this%albnirdir       => this%list%Get_ptr_r2d('albedo NIR direct')
      
      this%albnirdif       => this%list%Get_ptr_r2d('albedo NIR diffuse')
      
      this%albedo          => this%list%Get_ptr_r2d('albedo')
      

      this%q_ice_top       => this%list%Get_ptr_r2d('energy flux available for surface melting of sea ice')
      
      this%q_ice_bot       => this%list%Get_ptr_r2d('energy flux at ice-ocean interface')
      
      this%snow_thickness  => this%list%Get_ptr_r2d('thickness of snow on sea ice')
      

      this%t2m             => this%list%Get_ptr_r2d('2m temperature')
      
      this%t2m_tile        => this%list%Get_ptr_r3d('2m temperature, tile')
      
      this%wind10m         => this%list%Get_ptr_r2d('10m wind speed')
      
      this%u10m            => this%list%Get_ptr_r2d('10m zonal wind')
      
      this%v10m            => this%list%Get_ptr_r2d('10m meridional wind')
      
      this%wind10m_tile    => this%list%Get_ptr_r3d('10m wind speed, tile')
      
      this%u10m_tile       => this%list%Get_ptr_r3d('10m zonal wind, tile')
      
      this%v10m_tile       => this%list%Get_ptr_r3d('10m meridional wind, tile')
      

      this%nvalid          => this%list%get_ptr_i2d('number of valid points per block, tile')
      
      this%indices         => this%list%get_ptr_i3d('indices of valid points on chunk per block, tile')
      

    END SELECT

  END SUBROUTINE Set_pointers_diagnostics

END MODULE mo_vdf_sfc
