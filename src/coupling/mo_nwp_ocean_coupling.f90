!
! Interface between NWP physics and the ocean, through a coupler
!
! Notes on openMP parallelisation:
!   Most jb-loops cannot be optimized easily because converting a 2-D field into a 1-D field
!   requires remembering the 1-D index "ncount", which is incremented over the loop.  This
!   technique is used in the current code and openMP is not used.
!   A solution suggested by Rene Redler would be to calculate that index (nn) would look like this:
!      !ICON_OMP_PARALLEL_DO PRIVATE(jb, ic, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
!      DO jb = i_startblk, i_endblk
!        nn = (jb-1)*nproma                               ! translation to 1-d buffer fields
!        DO ic = 1, ext_data%atm%list_sea%ncount(jb)      ! number of ocean points (open water & sea ice)
!          jc = ext_data%atm%list_sea%idx(ic,jb)
!          prm_field(jg)%ocv(n,i_blk) = buffer(nn+jc,1)
!  It might also be necessary to synch the data before each loop passing data to the ocean.
!      CALL sync_patch_array(sync_c, p_patch, prm_diag%swflxsfc_t (:,:,isub_water) )
!
! Note: The variable names and numbers need to be consistent in 3 files:
!        - XML file: (supplied to model in run script)
!            <transient id="6" transient_standard_name="sea_surface_temperature"/>
!            The name will be used find the variable in mo_atmo_coupling_frame, not the number.
!        - mo_atmo_coupling_frame:
!            Variable names are associated to a variable number.
!            field_name(6) = "sea_surface_temperature"
!        - mo_nwp_ocean_coupling:
!            CALL yac_fget ( field_id(6), ... )
!            The numbers have to be consistent in both fortran files.
!       Component names in coupling.xml must (!) match with modelname_list[*].
!
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

MODULE mo_nwp_ocean_coupling

  USE mo_bc_greenhouse_gases ,ONLY: ghg_co2vmr
  USE mo_ccycle_config       ,ONLY: ccycle_config,                          &
       & CCYCLE_MODE_NONE, CCYCLE_MODE_PRESCRIBED, CCYCLE_MODE_INTERACTIVE, &
       & CCYCLE_CO2CONC_CONST, CCYCLE_CO2CONC_FROMFILE
  USE mo_dbg_nml             ,ONLY: idbg_mxmn, idbg_val
  USE mo_exception           ,ONLY: warning, message, finish
  USE mo_ext_data_types      ,ONLY: t_external_data
  USE mo_fortran_tools       ,ONLY: assert_acc_host_only, init
  USE mo_idx_list            ,ONLY: t_idx_list_blocked
  USE mo_impl_constants      ,ONLY: start_prog_cells, end_prog_cells
  USE mo_kind                ,ONLY: wp
  USE mo_lnd_nwp_config      ,ONLY: isub_water, isub_seaice, isub_lake,     &
       &                            hice_max
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_nonhydro_types      ,ONLY: t_nh_diag
  USE mo_nwp_phy_types       ,ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_types       ,ONLY: t_wtr_prog, t_lnd_diag
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_physical_constants  ,ONLY: vmr_to_mmr_co2
  USE mo_run_config          ,ONLY: ltimer
  USE mo_timer               ,ONLY: timer_start, timer_stop,                &
       &                            timer_coupling_put, timer_coupling_get, &
       &                            timer_coupling_1stget
  USE mo_util_dbg_prnt       ,ONLY: dbg_print


  IMPLICIT NONE

  PRIVATE

  !>
  !! Package of fields passed to the ocean via YAC.
  TYPE t_nwp_ocean_fields_tx

    !> Sea-water fraction [m2(ocean)/m2(gridcell)].
    REAL(wp), CONTIGUOUS, POINTER :: frac_w(:,:) => NULL()
    !> Sea-ice fraction [m2(seaice)/m2(gridcell)].
    REAL(wp), CONTIGUOUS, POINTER :: frac_i(:,:) => NULL()

    !> Resolved zonal surface stress over sea water [N/m2].
    REAL(wp), CONTIGUOUS, POINTER :: umfl_s_w(:,:) => NULL()
    !> Resolved zonal surface stress over sea ice [N/m2].
    REAL(wp), CONTIGUOUS, POINTER :: umfl_s_i(:,:) => NULL()
    !> Resolved meridional surface stress over sea water [N/m2].
    REAL(wp), CONTIGUOUS, POINTER :: vmfl_s_w(:,:) => NULL()
    !> Resolved meridional surface stress over sea ice [N/m2].
    REAL(wp), CONTIGUOUS, POINTER :: vmfl_s_i(:,:) => NULL()

    !> Moisture flux at sea water surface [kg/m2/s].
    REAL(wp), CONTIGUOUS, POINTER :: qhfl_s_w(:,:) => NULL()
    !> Moisture flux at sea ice surface [kg/m2/s].
    REAL(wp), CONTIGUOUS, POINTER :: qhfl_s_i(:,:) => NULL()

    !> Sensible heat flux at sea water surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: shfl_s_w(:,:) => NULL()
    !> Sensible heat flux at sea ice surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: shfl_s_i(:,:) => NULL()

    !> Latent heat flux at sea water surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: lhfl_s_w(:,:) => NULL()
    !> Latent heat flux at sea ice surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: lhfl_s_i(:,:) => NULL()

    !> Conductive heat flux at water-ice interface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: chfl_i(:,:) => NULL()

    !> 10m wind speed [m/s].
    REAL(wp), CONTIGUOUS, POINTER :: sp_10m(:,:) => NULL()

    !> Surface pressure [Pa].
    REAL(wp), CONTIGUOUS, POINTER :: pres_sfc(:,:) => NULL()

    !> Net short-wave flux at sea water surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: swflxsfc_w(:,:) => NULL()
    !> Net short-wave flux at sea ice surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: swflxsfc_i(:,:) => NULL()

    !> Net long-wave flux at sea water surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: lwflxsfc_w(:,:) => NULL()
    !> Net long-wave flux at sea ice surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: lwflxsfc_i(:,:) => NULL()

    !> Total rain rate [kg/m2/s].
    REAL(wp), CONTIGUOUS, POINTER :: rain_rate(:,:) => NULL()
    !> Total snow rate [kg/m2/s].
    REAL(wp), CONTIGUOUS, POINTER :: snow_rate(:,:) => NULL()

    !> Surface CO2 concentration [kg(CO2)/kg(air)].
    !! Not contiguous because it is the lowest level of a 3D field.
    REAL(wp), POINTER :: q_co2(:,:) => NULL()

  END TYPE t_nwp_ocean_fields_tx


  !>
  !! Package of fields received from the ocean via YAC.
  TYPE t_nwp_ocean_fields_rx

    !> Sea-surface temperature [K].
    REAL(wp), CONTIGUOUS, POINTER :: t_seasfc(:,:) => NULL()

    !> Sea-ice fraction [m2(seaice)/m2(ocean)].
    REAL(wp), CONTIGUOUS, POINTER :: fr_seaice(:,:) => NULL()

    !> Sea-ice thickness [m].
    REAL(wp), CONTIGUOUS, POINTER :: h_ice(:,:) => NULL()

    !> Zonal ocean surface velocity (optional) [m/s].
    REAL(wp), CONTIGUOUS, POINTER :: ocean_u(:,:) => NULL()

    !> Meridional ocean surface velocity (optional) [m/s].
    REAL(wp), CONTIGUOUS, POINTER :: ocean_v(:,:) => NULL()

    !> CO2 surface flux [kg/m2/s].
    REAL(wp), CONTIGUOUS, POINTER :: flx_co2(:,:) => NULL()

  END TYPE t_nwp_ocean_fields_rx


  PUBLIC :: nwp_couple_ocean
  PUBLIC :: couple_ocean
  PUBLIC :: t_nwp_ocean_fields_rx
  PUBLIC :: t_nwp_ocean_fields_tx

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_nwp_ocean_coupling' ! Output of module for debug

CONTAINS

  !>
  !! SUBROUTINE nwp_couple_ocean -- the interface between
  !! NWP physics and the ocean, through a coupler
  !!
  !! This subroutine is called from nwp_nh_interface.
  SUBROUTINE nwp_couple_ocean ( &
      & p_patch, pt_diag, lnd_diag, wtr_prog_new, prm_diag, ext_data, lacc &
    )

    TYPE(t_patch),                INTENT(INOUT)  :: p_patch
    TYPE(t_nh_diag),      TARGET, INTENT(INOUT)  :: pt_diag
    TYPE(t_wtr_prog),     TARGET, INTENT(INOUT)  :: wtr_prog_new
    TYPE(t_lnd_diag),     TARGET, INTENT(INOUT)  :: lnd_diag
    TYPE(t_nwp_phy_diag), TARGET, INTENT(INOUT)  :: prm_diag
    TYPE(t_external_data),        INTENT(IN)     :: ext_data
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    TYPE(t_nwp_ocean_fields_tx) :: tx
    TYPE(t_nwp_ocean_fields_rx) :: rx

    REAL(wp), TARGET :: rain_rate(nproma, p_patch%nblks_c)
    REAL(wp), TARGET :: snow_rate(nproma, p_patch%nblks_c)

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: jb, jc

    LOGICAL :: have_ice, have_hail, have_graupel


    CALL assert_acc_host_only('nwp_couple_ocean', lacc)

    CALL finish('nwp_couple_ocean', 'unintentionally called. Check your source code and configure.')
  END SUBROUTINE nwp_couple_ocean


  !>
  !! SUBROUTINE couple_ocean -- the actual interface between
  !! NWP physics and the ocean, through a coupler
  !!
  SUBROUTINE couple_ocean( p_patch, list_sea, tx, rx, lacc )

    ! Arguments

    TYPE(t_patch), INTENT(IN) :: p_patch
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_sea
    TYPE(t_nwp_ocean_fields_tx), INTENT(IN) :: tx
    TYPE(t_nwp_ocean_fields_rx), INTENT(INOUT) :: rx
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: jg                    ! grid index
    INTEGER               :: jb                    ! block loop count
    INTEGER               :: jc                    ! nproma loop count
    INTEGER               :: ic                    ! nproma loop count
    INTEGER               :: info, ierror          ! return values from cpl_put/get calls
    INTEGER               :: i_startblk, i_endblk  ! blocks
    INTEGER               :: i_startidx, i_endidx  ! slices

    REAL(wp), TARGET      :: buf(nproma, p_patch%nblks_c)


    REAL(wp), PARAMETER   :: csmall = 1.0E-5_wp    ! small number (security constant)

    REAL(wp) :: co2conc

    CALL assert_acc_host_only('couple_ocean', lacc)

    CALL finish('couple_ocean: unintentionally called. Check your source code and configure.')

  END SUBROUTINE couple_ocean


END MODULE mo_nwp_ocean_coupling
