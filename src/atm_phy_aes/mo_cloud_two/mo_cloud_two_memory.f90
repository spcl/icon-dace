!
! Global variables for the two-moment bulk microphysics by Seifert and Beheng (2006)
!                  with prognostic cloud droplet number parameterization
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

MODULE mo_cloud_two_memory

  USE mo_exception               ,ONLY: message, finish
  USE mtime                      ,ONLY: timedelta, OPERATOR(>)

  USE mo_model_domain            ,ONLY: t_patch
  USE mo_parallel_config         ,ONLY: nproma
  USE mo_time_config             ,ONLY: get_dynamics_timestep
  USE mo_aes_phy_config          ,ONLY: aes_phy_tc, dt_zero
  USE mo_io_config               ,ONLY: lnetcdf_flt64_output
  USE mo_name_list_output_config ,ONLY: is_variable_in_output

  USE mo_impl_constants          ,ONLY: success, vintp_method_lin, vintp_method_lin_nlevp1
  USE mo_cdi_constants           ,ONLY: grid_unstructured_cell, grid_cell
  USE mo_cdi                     ,ONLY: grid_unstructured,                 &
       &                                datatype_pack16,                   &
       &                                datatype_flt32,  datatype_flt64,   &
       &                                datatype_int32,                    &
       &                                tstep_instant, tstep_constant

  USE mo_var_list                ,ONLY: add_var, t_var_list_ptr
  USE mo_var_list_register       ,ONLY: vlr_add, vlr_del
  USE mo_var_metadata            ,ONLY: create_vert_interp_metadata, vintp_types
  USE mo_cf_convention           ,ONLY: t_cf_var
  USE mo_grib2                   ,ONLY: grib2_var
  USE mo_zaxis_type              ,ONLY: za_reference, za_surface, za_reference_half

  USE mo_cloud_two_types         ,ONLY: t_cloud_two_input, t_cloud_two_output

  ! include definition for "__acc_attach(ptr)"
#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cloud_two_input, cloud_two_output
  PUBLIC :: cloud_two_list
  PUBLIC :: construct_cloud_two_memory
  PUBLIC :: destruct_cloud_two_memory

  CHARACTER(len=*), PARAMETER                    :: thismodule = 'mo_cloud_two_memory'
  TYPE(t_cloud_two_input)  , ALLOCATABLE, TARGET :: cloud_two_input(:)  !< shape: (ng)
  TYPE(t_cloud_two_output) , ALLOCATABLE, TARGET :: cloud_two_output(:) !< shape: (ng)
  TYPE(t_var_list_ptr)     , ALLOCATABLE         :: cloud_two_list(:)   !< shape: (ng)

CONTAINS

  !!--------------------------------------------------------------------------
  !!                SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
  !!--------------------------------------------------------------------------
  !>
  !! Top-level procedure for building the state
  !!
  SUBROUTINE construct_cloud_two_memory( patch_array )

    TYPE(t_patch),INTENT(IN) :: patch_array(:)

    INTEGER :: ng, jg, ist
    INTEGER :: nlev, nblks
    TYPE(timedelta) :: dt_dyn

    !---

    CALL message(thismodule,'Construction of cloud_two_list, cloud_two_input and cloud two_output started.')

    ! allocate pointer arrays for the pointer arrays cloud_two_memory and cloud_two_list

    ng = SIZE(patch_array)

    ALLOCATE( cloud_two_list(ng), STAT=ist)
    IF (ist/=success) CALL finish(thismodule, 'allocation of cloud_two_list(ng) failed')

    ALLOCATE( cloud_two_input(ng), STAT=ist)
    IF (ist/=success) CALL finish(thismodule, 'allocation of cloud_two_input(ng) failed')

    ALLOCATE( cloud_two_output(ng), STAT=ist)
    IF (ist/=success) CALL finish(thismodule, 'allocation of cloud_two_output(ng) failed')

    !$ACC ENTER DATA CREATE(cloud_two_input, cloud_two_output)

    ! build lists and allocate memory for all grids where the 2 moment cloud microphysics is used

    DO jg = 1,ng
       IF (aes_phy_tc(jg)%dt_two > dt_zero) THEN
          !
          nlev   = patch_array(jg)%nlev
          nblks  = patch_array(jg)%nblks_c
          dt_dyn = get_dynamics_timestep(patch_array(jg))
          !
          CALL construct_cloud_two_list( jg,                  &
               &                         nproma, nlev, nblks, &
               &                         dt_dyn,              &
               &                         cloud_two_list(jg),  &
               &                         cloud_two_input(jg), &
               &                         cloud_two_output(jg) )
          !
       END IF
    END DO

    CALL message(thismodule,'Construction of cloud_two_list, cloud_two_input and cloud two_output finished.')

  END SUBROUTINE construct_cloud_two_memory
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !>
  !! Release memory used by the state variable arrays and list arrays
  !!
  SUBROUTINE destruct_cloud_two_memory

    INTEGER :: ng   !< total number of grids
    INTEGER :: jg   !< grid index
    INTEGER :: ist  !< system status code

    !---
    CALL message(thismodule,'Destruction of cloud_two_list, cloud_two_input and cloud_two_output started.')

    ng = SIZE(cloud_two_input)

    DO jg = 1,ng
       IF (aes_phy_tc(jg)%dt_two > dt_zero) THEN
          !
          CALL vlr_del( cloud_two_list(jg) )
          !
       END IF
    END DO

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(cloud_two_input, cloud_two_output)

    DEALLOCATE( cloud_two_list, STAT=ist )
    IF (ist/=success) CALL finish(thismodule, 'deallocation of cloud_two_list failed')

    DEALLOCATE( cloud_two_input, STAT=ist )
    IF (ist/=success) CALL finish(thismodule, 'deallocation of cloud_two_input failed')

    DEALLOCATE( cloud_two_output, STAT=ist )
    IF (ist/=success) CALL finish(thismodule, 'deallocation of cloud_two_output failed')

    CALL message(thismodule,'Destruction of cloud_two_list, cloud_two_input and cloud two_output finished.')

  END SUBROUTINE destruct_cloud_two_memory
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !>
  SUBROUTINE construct_cloud_two_list( jg,                  &
       &                               nproma, nlev, nblks, &
       &                               dt_dyn,              &
       &                               cloud_two_list,      &
       &                               cloud_two_input,     &
       &                               cloud_two_output )

    INTEGER                 , INTENT(in)    :: jg                  !< grid index
    INTEGER                 , INTENT(in)    :: nproma, nlev, nblks !< size of dimensions
    TYPE(timedelta)         , INTENT(in)    :: dt_dyn              !< Dynamics timestep

    TYPE(t_var_list_ptr)    , INTENT(inout) :: cloud_two_list      !< pointers for list of variables
    TYPE(t_cloud_two_input) , INTENT(inout) :: cloud_two_input     !< pointers for input variables

    TYPE(t_cloud_two_output), INTENT(inout) :: cloud_two_output    !< pointers for output variables

    ! Local variables

    CHARACTER(len= 2) :: cg
    CHARACTER(len=20) :: listname

    INTEGER           :: shape2d(2), shape3d(3)
    INTEGER           :: datatype_grb
    INTEGER           :: datatype_flt, datatype_int

    WRITE(cg,'(i2.2)') jg
    CALL message('construct_cloud_two_list','create list and allocate memory for jg ='//cg)

    ! number of bits for data representation in grib2
    datatype_grb = datatype_pack16

    ! number of bits for data representation in netcdf
    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
    datatype_int = datatype_int32

    shape2d         = (/nproma,       nblks/)
    shape3d         = (/nproma, nlev, nblks/)

    ! define list name
    WRITE(listname,'(a,i2.2)') 'cloud_two_memory_D',jg

    ! register the cloud_two_list for grid jg
    CALL vlr_add( cloud_two_list                  ,&
         &        listname                        ,&
         &        patch_id  = jg                  ,&
         &        loutput   = .TRUE.              ,&
         &        lrestart  = .FALSE.             ,&
         &        linitial  = .FALSE. )

    ! Input parameters
    ! ----------------
    !
    ! These fields are constructed only if they are requested for output
    !
    IF ( is_variable_in_output(var_name='jcs_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'jcs_two'                                                ,&
            &        ptr         = cloud_two_input%jcs                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('column_start_index',                           &
            &                                '-',                                            &
            &                                'column start index (cloud_two input)',         &
            &                                datatype_int)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_constant                                           ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%jcs)
    END IF
    !
    IF ( is_variable_in_output(var_name='jce_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'jce_two'                                                ,&
            &        ptr         = cloud_two_input%jce                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('column_end_index',                             &
            &                                '-',                                            &
            &                                'column end index (cloud_two input)',           &
            &                                datatype_int)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_constant                                           ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%jce)
    END IF
    !
    IF ( is_variable_in_output(var_name='msg_level_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'msg_level_two'                                          ,&
            &        ptr         = cloud_two_input%msg_level                                ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('message_level',                                &
            &                                '-',                                            &
            &                                'message level (cloud_two input)',              &
            &                                datatype_int)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_constant                                           ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%msg_level)
    END IF
    !
    IF ( is_variable_in_output(var_name='pdtime_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'pdtime_two'                                             ,&
            &        ptr         = cloud_two_input%pdtime                                   ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('time_step',                                    &
            &                                's',                                            &
            &                                'time step (cloud_two input)',                  &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_constant                                           ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%pdtime)
    END IF

    !
    ! Input fields (1)
    ! -----------------
    !
    ! These fields are constructed only if they are requested for output
    !
    ! INTENT-IN variables of the parameterization
    !
    IF ( is_variable_in_output(var_name='dz_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'dz_two'                                                 ,&
            &        ptr         = cloud_two_input%dz                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('cell_thickness',                               &
            &                                'm',                                            &
            &                                'cell thickness (cloud_two input)',             &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,3,12,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%dz)
    END IF
    !
    IF ( is_variable_in_output(var_name='zh_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'zh_two'                                                 ,&
            &        ptr         = cloud_two_input%zh                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference_half                                        ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('height_of_half_levels',                        &
            &                                'm',                                            &
            &                                'height of half levels (cloud_two input)',      &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,3,12,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin_nlevp1)       ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%zh)
    END IF
    !
    IF ( is_variable_in_output(var_name='rho_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'rho_two'                                                ,&
            &        ptr         = cloud_two_input%rho                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('air_density',                                  &
            &                                'kg/m3',                                        &
            &                                'air density (cloud_two input)',                &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,3,10,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%rho)
    END IF
    !
    IF ( is_variable_in_output(var_name='pf_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'pf_two'                                                 ,&
            &        ptr         = cloud_two_input%pf                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('air_pressure',                                 &
            &                                'kg/m3',                                        &
            &                                'air pressure (cloud_two input)',               &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,3,0,                                          &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%pf)
    END IF
    !
    IF ( is_variable_in_output(var_name='cpair_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'cpair_two'                                              ,&
            &        ptr         = cloud_two_input%cpair                                    ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('air_specific_heat',                            &
            &                                'J/K/kg',                                       &
            &                                'specific heat of air at constant pressure '//  &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,0,255,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%cpair)
    END IF
    !
    ! Input fields (2)
    ! -----------------
    !
    ! These fields are constructed only if they are requested for output
    !
    ! INTENT-INOUT variables of the parameterization
    !
    IF ( is_variable_in_output(var_name='in_ta_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_ta_two'                                              ,&
            &        ptr         = cloud_two_input%ta                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('air_temperature',                              &
            &                                'K',                                            &
            &                                'air temperature (cloud_two input)',            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,0,0,                                          &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%ta)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qv_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qv_two'                                              ,&
            &        ptr         = cloud_two_input%qv                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('specific_humidity',                            &
            &                                'kg/kg',                                        &
            &                                'specific humidity (cloud_two input)',          &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,0,                                          &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qv)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qc_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qc_two'                                              ,&
            &        ptr         = cloud_two_input%qc                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_cloud_liquid_water_in_air',   &
            &                                'kg/kg',                                        &
            &                                'mass fraction of cloud liquid water in air '// &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,83,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qc)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qi_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qi_two'                                              ,&
            &        ptr         = cloud_two_input%qi                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_cloud_ice_in_air',            &
            &                                'kg/kg',                                        &
            &                                'mass fraction of cloud ice in air '//          &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,84,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qi)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qr_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qr_two'                                              ,&
            &        ptr         = cloud_two_input%qr                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_rain_in_air',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of rain in air '//               &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,85,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qr)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qs_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qs_two'                                              ,&
            &        ptr         = cloud_two_input%qs                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_snow_in_air',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of snow in air '//               &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,86,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qs)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qg_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qg_two'                                              ,&
            &        ptr         = cloud_two_input%qg                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_graupel_in_air',              &
            &                                'kg/kg',                                        &
            &                                'mass fraction of graupel in air '//            &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,32,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qg)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qh_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qh_two'                                              ,&
            &        ptr         = cloud_two_input%qh                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_hail_in_air',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of hail in air '//               &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,31,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qg)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qnc_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qnc_two'                                             ,&
            &        ptr         = cloud_two_input%qnc                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_cloud_droplet_number',        &
            &                                'kg/kg',                                        &
            &                                'mass fraction of cloud doplet number '//       &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,83,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qnc)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qni_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qni_two'                                             ,&
            &        ptr         = cloud_two_input%qni                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_cloud_ice_number',            &
            &                                'kg/kg',                                        &
            &                                'mass fraction of cloud ice number '//          &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,84,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qni)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qnr_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qnr_two'                                             ,&
            &        ptr         = cloud_two_input%qnr                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_rain_droplet_number',         &
            &                                'kg/kg',                                        &
            &                                'mass fraction of rain droplet number '//       &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,85,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qnr)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qns_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qns_two'                                             ,&
            &        ptr         = cloud_two_input%qns                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_snow_number',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of snow number '//               &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,86,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qns)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qng_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qng_two'                                             ,&
            &        ptr         = cloud_two_input%qng                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_graupel_number',              &
            &                                'kg/kg',                                        &
            &                                'mass fraction of graupel number '//            &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,32,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qng)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_qnh_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_qnh_two'                                             ,&
            &        ptr         = cloud_two_input%qnh                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_hail_number',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of hail number '//               &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,31,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%qng)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_ninact_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_ninact_two'                                          ,&
            &        ptr         = cloud_two_input%ninact                                   ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('number_of_activated_ice_nuclei',               &
            &                                'kg/kg',                                        &
            &                                'number of activated ice nuclei '//             &
            &                                '(cloud_two input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,31,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%ninact)
    END IF
    !
    IF ( is_variable_in_output(var_name='in_w_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'in_w_two'                                               ,&
            &        ptr         = cloud_two_input%w                                        ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference_half                                        ,&
            &        ldims       = (/nproma,nlev+1,nblks/)                                  ,&
            &        cf          = t_cf_var ('vertical_velocity',                            &
            &                                'Pa s-1',                                       &
            &                                'vertical velocity (cloud_two input)',          &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,2,8,                                          &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_input%w)
    END IF


    ! Output fields (1)
    ! -----------------
    !
    ! These fields are constructed only if they are requested for output
    ! or needed for recycling in the time stepping.
    !
    ! INTENT-OUT variables of the parameterization:
    ! tendencies in the atmosphere
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_ta_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_ta_two'                                            ,&
            &        ptr         = cloud_two_output%tend_ta_two                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_air_temperature_'//                &
            &                                'due_to_stratiform_cloud_and_precipitation',    &
            &                                'K/s',                                          &
            &                                'tendency of air temperature '//                &
            &                                'due to stratiform cloud and precipitation '//  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,0,193,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_ta_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qv_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qv_two'                                            ,&
            &        ptr         = cloud_two_output%tend_qv_two                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_specific_humidity_'//              &
            &                                'due_to_stratiform_cloud_and_precipitation',    &
            &                                '1/s',                                          &
            &                                'tendency of specific humidity '//              &
            &                                'due to stratiform cloud and precipitation '//  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,200,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qv_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qc_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qc_two'                                            ,&
            &        ptr         = cloud_two_output%tend_qc_two                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_mass_fraction_of_'//               &
            &                                'stratiform_cloud_liquid_water_in_air_'//       &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of mass fraction of '//               &
            &                                'stratiform cloud liquid water in air '//       &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,201,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qc_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qi_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qi_two'                                            ,&
            &        ptr         = cloud_two_output%tend_qi_two                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_mass_fraction_of_'//               &
            &                                'stratiform_cloud_ice_in_air_'//                &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of mass fraction of '//               &
            &                                'stratiform cloud ice in air '//                &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,202,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qi_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qr_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qr_two'                                            ,&
            &        ptr         = cloud_two_output%tend_qr_two                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_mass_fraction_of_'//               &
            &                                'rain_in_air_'//                                &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of mass fraction of '//               &
            &                                'rain in air '//                                &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,255,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qr_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qs_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qs_two'                                            ,&
            &        ptr         = cloud_two_output%tend_qs_two                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_mass_fraction_of_'//               &
            &                                'snow_in_air_'//                                &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of mass fraction of '//               &
            &                                'snow in air '//                                &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,255,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qs_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qg_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qg_two'                                            ,&
            &        ptr         = cloud_two_output%tend_qg_two                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_mass_fraction_of_'//               &
            &                                'graupel_in_air_'//                             &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of mass fraction of '//               &
            &                                'graupel in air '//                             &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,255,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qg_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qh_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qh_two'                                            ,&
            &        ptr         = cloud_two_output%tend_qh_two                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_mass_fraction_of_'//               &
            &                                'hail_in_air'//                                 &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of mass fraction of '//               &
            &                                'hail in air '//                                &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,255,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qh_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qnc_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qnc_two'                                           ,&
            &        ptr         = cloud_two_output%tend_qnc_two                            ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_number_of_'//                      &
            &                                'stratiform_cloud_droplets_'//                  &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of number of '//                      &
            &                                'stratiform cloud droplets '//                  &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,202,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qnc_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qni_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qni_two'                                           ,&
            &        ptr         = cloud_two_output%tend_qni_two                            ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_number_of_'//                      &
            &                                'stratiform_cloud_ice_droplets_'//              &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of number of '//                      &
            &                                'stratiform cloud ice droplets '//              &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,202,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qni_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qnr_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qnr_two'                                           ,&
            &        ptr         = cloud_two_output%tend_qnr_two                            ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_number_of_'//                      &
            &                                'rain_droplets_'//                              &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of number of '//                      &
            &                                'rain droplets '//                              &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,255,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qnr_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qns_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qns_two'                                           ,&
            &        ptr         = cloud_two_output%tend_qns_two                            ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_number_of_'//                      &
            &                                'snow_droplets'//                               &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of number of '//                      &
            &                                'snow droplets '//                              &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,255,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qns_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qng_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qng_two'                                           ,&
            &        ptr         = cloud_two_output%tend_qng_two                            ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_number_of_'//                      &
            &                                'graupel_droplets_'//                           &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of number of '//                      &
            &                                'graupel droplets '//                           &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,255,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qng_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_qnh_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_qnh_two'                                           ,&
            &        ptr         = cloud_two_output%tend_qnh_two                            ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_number_of_'//                      &
            &                                'hail_droplets_'//                              &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of number of '//                      &
            &                                'hail droplets '//                              &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,255,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_qnh_two)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_dyn .OR.                                                 &
         & is_variable_in_output(var_name='tend_ninact_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'tend_ninact_two'                                        ,&
            &        ptr         = cloud_two_output%tend_ninact_two                         ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tendency_of_number_of_'//                      &
            &                                'activated_ice_nuclei'//                        &
            &                                'due_to_cloud_microphysics',                    &
            &                                '1/s',                                          &
            &                                'tendency of number of '//                      &
            &                                'activated ice nuclei '//                       &
            &                                'due to cloud microphysics '//                  &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,255,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%tend_ninact_two)
    END IF
    ! INTENT-INOUT variables of the parameterization
    !
    IF ( is_variable_in_output(var_name='out_ta_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_ta_two'                                             ,&
            &        ptr         = cloud_two_output%ta                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('air_temperature',                              &
            &                                'K',                                            &
            &                                'air temperature (cloud_two output)',           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,0,0,                                          &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%ta)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qv_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qv_two'                                             ,&
            &        ptr         = cloud_two_output%qv                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('specific_humidity',                            &
            &                                'kg/kg',                                        &
            &                                'specific humidity (cloud_two output)',         &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,0,                                          &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qv)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qc_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qc_two'                                             ,&
            &        ptr         = cloud_two_output%qc                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_cloud_liquid_water_in_air',   &
            &                                'kg/kg',                                        &
            &                                'mass fraction of cloud liquid water in air '// &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,83,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qc)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qi_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qi_two'                                             ,&
            &        ptr         = cloud_two_output%qi                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_cloud_ice_in_air',            &
            &                                'kg/kg',                                        &
            &                                'mass fraction of cloud ice in air '//          &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,84,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qi)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qr_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qr_two'                                             ,&
            &        ptr         = cloud_two_output%qr                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_rain_in_air',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of rain in air '//               &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,85,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qr)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qs_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qs_two'                                             ,&
            &        ptr         = cloud_two_output%qs                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_snow_in_air',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of snow in air '//               &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,86,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qs)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qg_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qg_two'                                             ,&
            &        ptr         = cloud_two_output%qg                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_graupel_in_air',              &
            &                                'kg/kg',                                        &
            &                                'mass fraction of graupel in air '//            &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,32,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qg)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qh_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qh_two'                                             ,&
            &        ptr         = cloud_two_output%qh                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_hail_in_air',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of hail in air '//               &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,31,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qg)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qnc_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qnc_two'                                            ,&
            &        ptr         = cloud_two_output%qnc                                     ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_cloud_droplet_number',        &
            &                                'kg/kg',                                        &
            &                                'mass fraction of cloud doplet number '//       &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,83,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qnc)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qni_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qni_two'                                            ,&
            &        ptr         = cloud_two_output%qni                                     ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_cloud_ice_number',            &
            &                                'kg/kg',                                        &
            &                                'mass fraction of cloud ice number '//          &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,84,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qni)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qnr_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qnr_two'                                            ,&
            &        ptr         = cloud_two_output%qnr                                     ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_rain_droplet_number',         &
            &                                'kg/kg',                                        &
            &                                'mass fraction of rain droplet number '//       &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,85,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qnr)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qns_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qns_two'                                            ,&
            &        ptr         = cloud_two_output%qns                                     ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_snow_number',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of snow number '//               &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,86,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qns)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qng_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qng_two'                                            ,&
            &        ptr         = cloud_two_output%qng                                     ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_graupel_number',              &
            &                                'kg/kg',                                        &
            &                                'mass fraction of graupel number '//            &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,32,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qng)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_qnh_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_qnh_two'                                            ,&
            &        ptr         = cloud_two_output%qnh                                     ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_hail_number',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of hail number '//               &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,31,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%qng)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_ninact_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_ninact_two'                                         ,&
            &        ptr         = cloud_two_output%ninact                                  ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('number_of_activated_ice_nuclei',               &
            &                                'kg/kg',                                        &
            &                                'number of activated ice nuclei '//             &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,31,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%ninact)
    END IF
    !
    IF ( is_variable_in_output(var_name='out_w_two') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'out_w_two'                                              ,&
            &        ptr         = cloud_two_output%w                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference_half                                        ,&
            &        ldims       = (/nproma,nlev+1,nblks/)                                  ,&
            &        cf          = t_cf_var ('vertical_velocity',                            &
            &                                'Pa s-1',                                       &
            &                                'vertical velocity (cloud_two output)',         &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,2,8,                                          &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%w)
    END IF
    !
    ! fluxes at the surface
    !
    IF ( is_variable_in_output(var_name='pr_rain') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'pr_rain'                                                ,&
            &        ptr         = cloud_two_output%pr_rain                                 ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('stratiform_rainfall_flux',                     &
            &                                'kg/m2/s',                                      &
            &                                'stratiform rainfall flux '//                   &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,77,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%pr_rain)
    END IF
    !
    IF ( is_variable_in_output(var_name='pr_ice') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'pr_ice'                                                 ,&
            &        ptr         = cloud_two_output%pr_ice                                  ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('stratiform_icefall_flux',                      &
            &                                'kg/m2/s',                                      &
            &                                'stratiform icefall flux '//                    &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,78,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%pr_ice)
    END IF
    !
    IF ( is_variable_in_output(var_name='pr_snow') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'pr_snow'                                                ,&
            &        ptr         = cloud_two_output%pr_snow                                 ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('stratiform_snowfall_flux',                     &
            &                                'kg/m2/s',                                      &
            &                                'stratiform snowfall flux '//                   &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,56,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%pr_snow)
    END IF
    !
    IF ( is_variable_in_output(var_name='pr_grpl') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'pr_grpl'                                                ,&
            &        ptr         = cloud_two_output%pr_grpl                                 ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('stratiform_graupel_flux',                      &
            &                                'kg/m2/s',                                      &
            &                                'stratiform graupel flux '//                    &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,75,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%pr_grpl)
   END IF
    !
    IF ( is_variable_in_output(var_name='pr_hail') ) THEN
       CALL add_var( this_list   = cloud_two_list                                           ,&
            &        varname     = 'pr_hail'                                                ,&
            &        ptr         = cloud_two_output%pr_hail                                 ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('stratiform_hail_flux',                         &
            &                                'kg/m2/s',                                      &
            &                                'stratiform hail flux '//                       &
            &                                '(cloud_two output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,75,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_two_output%pr_hail)
   END IF
    !

  END SUBROUTINE construct_cloud_two_list
  !-------------

END MODULE mo_cloud_two_memory
