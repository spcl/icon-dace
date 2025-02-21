! Read namelists, make sanity checks specific to each namelist and make
! a cross check once all namelists of a component are available.
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

MODULE mo_read_namelists

  USE mo_mpi                 ,ONLY: my_process_is_stdio
  USE mo_namelist            ,ONLY: open_nml_output, close_nml_output
  USE mo_nml_annotate        ,ONLY: log_nml_settings

  USE mo_time_nml            ,ONLY: read_time_namelist

  USE mo_parallel_nml        ,ONLY: read_parallel_namelist
  USE mo_run_nml             ,ONLY: read_run_namelist
  USE mo_io_nml              ,ONLY: read_io_namelist
  USE mo_gribout_nml         ,ONLY: read_gribout_namelist
  USE mo_dbg_nml             ,ONLY: read_dbg_namelist


  USE mo_grid_nml            ,ONLY: read_grid_namelist
  USE mo_grid_config         ,ONLY: init_grid_configuration
  USE mo_gridref_nml         ,ONLY: read_gridref_namelist
  USE mo_dynamics_nml        ,ONLY: read_dynamics_namelist
  USE mo_interpol_nml        ,ONLY: read_interpol_namelist
  USE mo_sleve_nml           ,ONLY: read_sleve_namelist
  USE mo_nonhydrostatic_nml  ,ONLY: read_nonhydrostatic_namelist
  USE mo_diffusion_nml       ,ONLY: read_diffusion_namelist

  USE mo_advection_nml       ,ONLY: read_transport_namelist

  USE mo_aes_phy_nml         ,ONLY: process_aes_phy_nml
  USE mo_aes_cov_nml         ,ONLY: process_aes_cov_nml
  USE mo_aes_cop_nml         ,ONLY: process_aes_cop_nml
  USE mo_aes_wmo_nml         ,ONLY: process_aes_wmo_nml
  USE mo_aes_rad_nml         ,ONLY: process_aes_rad_nml
  USE mo_aes_vdf_nml         ,ONLY: process_aes_vdf_nml
  USE mo_ccycle_nml          ,ONLY: process_ccycle_nml
  
  USE mo_nwp_phy_nml         ,ONLY: read_nwp_phy_namelist
  USE mo_nwp_tuning_nml      ,ONLY: read_nwp_tuning_namelist
  USE mo_ensemble_pert_nml   ,ONLY: read_ensemble_pert_namelist
  USE mo_radiation_nml       ,ONLY: read_radiation_namelist
  USE mo_synsat_nml          ,ONLY: read_synsat_namelist
  USE mo_synradar_nml        ,ONLY: read_synradar_namelist
  USE mo_turbdiff_nml        ,ONLY: read_turbdiff_namelist
  USE mo_lnd_nwp_nml         ,ONLY: read_nwp_lnd_namelist
  USE mo_art_nml             ,ONLY: read_art_namelist

  USE mo_turb_vdiff_nml      ,ONLY: vdiff_read_namelist

  USE mo_2mom_mcrph_nml      ,ONLY: read_2mom_mcrph_namelist

  USE mo_initicon_nml        ,ONLY: read_initicon_namelist
  USE mo_nh_testcases_nml    ,ONLY: read_nh_testcase_namelist, nh_test_name
  USE mo_aes_bubble_nml      ,ONLY: process_aes_bubble_nml
  USE mo_scm_nml             ,ONLY: read_scm_namelist
  USE mo_meteogram_nml       ,ONLY: read_meteogram_namelist

  USE mo_coupling_nml        ,ONLY: read_coupling_namelist
  USE mo_extpar_nml          ,ONLY: read_extpar_namelist

  USE mo_sea_ice_nml         ,ONLY: read_sea_ice_namelist

  USE mo_name_list_output_init ,ONLY: read_name_list_output_namelists

  USE mo_les_nml             ,ONLY: read_les_namelist
  USE mo_ls_forcing_nml      ,ONLY: read_ls_forcing_namelist

  USE mo_limarea_nml         ,ONLY: read_limarea_namelist

  USE mo_run_config          ,ONLY: iforcing
  USE mo_impl_constants      ,ONLY: iaes, ILDF_ECHAM, INWP
  USE mo_assimilation_nml    ,ONLY: read_assimilation_namelist
  USE mo_nudging_nml         ,ONLY: read_nudging_namelist
  USE mo_upatmo_nml          ,ONLY: read_upatmo_namelist
  USE mo_ser_nml             ,ONLY: read_ser_namelist

  USE mo_sppt_nml            ,ONLY: read_sppt_namelist
  USE mo_comin_nml           ,ONLY: read_comin_namelist

  !OEM
  USE mo_oem_nml             ,ONLY: read_oemctrl_namelist

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_atmo_namelists

CONTAINS

  !---------------------------------------------------------------------
  !>
  !! Read namelists for atmospheric models
  !!
  SUBROUTINE read_atmo_namelists(atm_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: atm_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    INTEGER :: tlen
    LOGICAL :: is_stdio

    !-----------------------------------------------------------------
    ! Create a new file in which all the namelist variables and their
    ! actual values used in the model run will be stored.
    !-----------------------------------------------------------------
    is_stdio = my_process_is_stdio()
    IF (is_stdio) CALL open_nml_output('NAMELIST_ICON_output_atm')

    !-----------------------------------------------------------------
    ! Read namelists that are shared by all components of the model.
    ! This means that the same namelists with the same values are
    ! read by all components of a coupled system.
    !-----------------------------------------------------------------

    CALL read_time_namelist           (TRIM(shr_namelist_filename))

    !-----------------------------------------------------------------
    ! Read namelist that are specific to the atm model.
    ! In case of a coupled simulation, the ocean model may also
    ! read some of these namelists, but probably from a different
    ! ASCII file containing different values.
    !-----------------------------------------------------------------

    ! General
    !

    tlen = LEN_TRIM(atm_namelist_filename)
    CALL read_parallel_namelist       (atm_namelist_filename(1:tlen))
    CALL read_run_namelist            (atm_namelist_filename(1:tlen))
    CALL read_io_namelist             (atm_namelist_filename(1:tlen))
    CALL read_meteogram_namelist      (atm_namelist_filename(1:tlen))
    CALL read_name_list_output_namelists (atm_namelist_filename(1:tlen))
    CALL read_dbg_namelist            (atm_namelist_filename(1:tlen))
    CALL read_synsat_namelist         (atm_namelist_filename(1:tlen))
    CALL read_synradar_namelist       (atm_namelist_filename(1:tlen))

    ! Grid
    !
    CALL read_grid_namelist           (atm_namelist_filename(1:tlen))
    CALL read_gridref_namelist        (atm_namelist_filename(1:tlen))
    CALL read_interpol_namelist       (atm_namelist_filename(1:tlen))
    CALL read_sleve_namelist          (atm_namelist_filename(1:tlen))
    !
    CALL init_grid_configuration()    ! so that the number of grids is known
    !                                 ! and arrays can be allocated

    ! Dynamics
    !
    CALL read_dynamics_namelist       (atm_namelist_filename(1:tlen))
    CALL read_nonhydrostatic_namelist (atm_namelist_filename(1:tlen))
    CALL read_diffusion_namelist      (atm_namelist_filename(1:tlen))

    ! Transport
    !
    CALL read_transport_namelist      (atm_namelist_filename(1:tlen))

    ! Physics
    !
    SELECT CASE (iforcing)
    CASE (iaes, ILDF_ECHAM)
       !
       ! AES physics ...
       CALL process_aes_phy_nml          (atm_namelist_filename(1:tlen))
       !
       ! ... and the employed parameterizations
       CALL process_aes_cov_nml          (atm_namelist_filename(1:tlen))
       CALL process_aes_cop_nml          (atm_namelist_filename(1:tlen))
       CALL process_aes_wmo_nml          (atm_namelist_filename(1:tlen))
       CALL process_aes_rad_nml          (atm_namelist_filename(1:tlen))
       CALL process_aes_vdf_nml          (atm_namelist_filename(1:tlen))
       !
       ! Carbon cycle for AES physics/JSBACH/HAMOCC
       CALL process_ccycle_nml           (atm_namelist_filename(1:tlen))
       !
       CALL read_sea_ice_namelist        (atm_namelist_filename(1:tlen))
       CALL read_art_namelist            (atm_namelist_filename(1:tlen))
       !
    CASE (INWP)
       !
       CALL read_nwp_phy_namelist        (atm_namelist_filename(1:tlen))
       CALL read_nwp_tuning_namelist     (atm_namelist_filename(1:tlen))
       CALL read_ensemble_pert_namelist  (atm_namelist_filename(1:tlen))
       CALL read_radiation_namelist      (atm_namelist_filename(1:tlen))
       CALL read_turbdiff_namelist       (atm_namelist_filename(1:tlen))
       CALL read_nwp_lnd_namelist        (atm_namelist_filename(1:tlen))

       CALL vdiff_read_namelist          (atm_namelist_filename(1:tlen))
       CALL process_ccycle_nml           (atm_namelist_filename(1:tlen))

       CALL read_art_namelist            (atm_namelist_filename(1:tlen))
       CALL read_2mom_mcrph_namelist     (atm_namelist_filename(1:tlen))

       CALL read_les_namelist            (atm_namelist_filename(1:tlen))
       CALL read_ls_forcing_namelist     (atm_namelist_filename(1:tlen))

       CALL read_sppt_namelist           (atm_namelist_filename(1:tlen))
       !
    END SELECT

    ! Upper atmosphere
    !
    CALL read_upatmo_namelist         (TRIM(atm_namelist_filename))

    ! Initial conditions
    !
    CALL read_initicon_namelist       (atm_namelist_filename(1:tlen))
    CALL read_nh_testcase_namelist    (atm_namelist_filename(1:tlen))
    IF (nh_test_name(1:10) == 'aes_bubble') THEN
      CALL process_aes_bubble_nml     (atm_namelist_filename(1:tlen))
    END IF  
    
    CALL read_scm_namelist            (atm_namelist_filename(1:tlen))

    ! Boundary conditions
    !
    CALL read_extpar_namelist         (atm_namelist_filename(1:tlen))
    CALL read_limarea_namelist        (atm_namelist_filename(1:tlen))
    CALL read_nudging_namelist        (atm_namelist_filename(1:tlen))

    ! GRIB output
    CALL read_gribout_namelist        (atm_namelist_filename(1:tlen))

    ! Coupling
    !
    CALL read_coupling_namelist       (atm_namelist_filename(1:tlen))

    ! Assimilation
    CALL read_assimilation_namelist   (atm_namelist_filename(1:tlen))

    !OEM
    CALL read_oemctrl_namelist        (atm_namelist_filename(1:tlen))

    !COMIN
    !
    !remark: every ICON component defines its own ComIn namelist,
    !        ie. `atmo_model` and `ocean_model`.
    CALL read_comin_namelist          (atm_namelist_filename(1:tlen))

    ! Serialization
    !$ser verbatim CALL read_ser_namelist(atm_namelist_filename(1:tlen))
    !-----------------------------------------------------------------
    ! Close the file in which all the namelist variables and their
    ! actual values were stored.
    !-----------------------------------------------------------------

    IF (is_stdio) CALL close_nml_output

    ! write an annotate table of all namelist settings to a text file
    IF (is_stdio) CALL log_nml_settings("nml.atmo.log")

  END SUBROUTINE read_atmo_namelists
  !-------------------------------------------------------------------------


END MODULE mo_read_namelists
