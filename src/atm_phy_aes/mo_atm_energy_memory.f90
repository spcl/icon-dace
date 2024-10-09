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

!> @brief
!! Module for variables used for the atmosphere energy diagnostics
!!
!! Author: Marco Giorgetta, MPI-M, 2024

MODULE mo_atm_energy_memory

  USE mo_exception               ,ONLY: message, finish

  USE mo_model_domain            ,ONLY: t_patch

  USE mo_parallel_config         ,ONLY: nproma
  USE mo_io_config               ,ONLY: lnetcdf_flt64_output
  USE mo_gribout_config          ,ONLY: gribout_config
  USE mo_name_list_output_config ,ONLY: is_variable_in_output

  USE mo_master_control          ,ONLY: get_my_process_name
  USE mo_var_list                ,ONLY: add_var, t_var_list_ptr
  USE mo_var_list_register       ,ONLY: vlr_add, vlr_del

  USE mo_impl_constants          ,ONLY: success
  USE mo_cdi_constants           ,ONLY: grid_unstructured_cell, grid_cell
  USE mo_cdi                     ,ONLY: grid_unstructured, grid_lonlat,    &
       &                                datatype_pack16, datatype_pack24,  &
       &                                datatype_flt32,  datatype_flt64,   &
       &                                datatype_int32,                    &
       &                                tstep_instant, tstep_constant
  USE mo_cf_convention           ,ONLY: t_cf_var
  USE mo_grib2                   ,ONLY: t_grib2_var, grib2_var
  USE mo_zaxis_type              ,ONLY: za_reference, za_surface

  USE mo_atm_energy_types        ,ONLY: t_atm_energy, t_atm_energy_config

  ! include definition for "__acc_attach(ptr)"
#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: atm_energy_config, atm_energy_list, atm_energy
  PUBLIC :: construct_atm_energy, destruct_atm_energy

  CHARACTER(len=*)         , PARAMETER           :: thismodule = 'mo_atm_energy_memory'

  TYPE(t_atm_energy_config), ALLOCATABLE         :: atm_energy_config(:)
  TYPE(t_var_list_ptr)     , ALLOCATABLE         :: atm_energy_list(:)
  TYPE(t_atm_energy)       , ALLOCATABLE, TARGET :: atm_energy(:)

CONTAINS

  !!--------------------------------------------------------------------------
  !!                SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
  !!--------------------------------------------------------------------------
  !>
  !! Top-level procedure for building the state
  !!
  SUBROUTINE construct_atm_energy(patch_array)

    TYPE(t_patch), INTENT(IN)  :: patch_array(:)

    INTEGER :: ng, jg, ist
    INTEGER :: nlev, nblks

    !---

    CALL message(thismodule,'Construction of atm_energy_config, atm_energy_list and atm_energy started.')

    !> number of grids
    !
    ng = SIZE(patch_array)

    !> allocate configuration for ng grids
    !
    ALLOCATE(atm_energy_config(ng), STAT=ist)
    IF (ist/=success) CALL finish(thismodule, 'allocation of atm_energy_config(ng) failed')

    !> allocate lists for ng grids
    !
    ALLOCATE(atm_energy_list(ng), STAT=ist)
    IF (ist/=success) CALL finish(thismodule, 'allocation of atm_energy_list(ng) failed')

    !> allocate memory for ng grids
    !
    ALLOCATE(atm_energy(ng), STAT=ist)
    IF (ist/=success) CALL finish(thismodule, 'allocation of atm_energy(ng) failed')

    !$ACC ENTER DATA CREATE(atm_energy)

    !> build lists and memory for ng grids
    !
    DO jg = 1, ng
       !
       nlev   = patch_array(jg)%nlev
       nblks  = patch_array(jg)%nblks_c
       !
       CALL construct_atm_energy_list(jg,                    &
            &                         nproma, nlev, nblks,   &
            &                         atm_energy_config(jg), &
            &                         atm_energy_list(jg),   &
            &                         atm_energy(jg))
       !
    END DO

    CALL message(thismodule,'Construction of atm_energy_config, atm_energy_list and atm_energy finished.')

  END SUBROUTINE construct_atm_energy
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !>
  !! Release memory used by the state variable arrays and list arrays
  !!
  SUBROUTINE destruct_atm_energy

    INTEGER :: ng   !< total number of grids
    INTEGER :: jg   !< grid index
    INTEGER :: ist  !< system status code

    !---
    CALL message(thismodule,'Destruction of atm_energy_config, atm_energy_list and atm_energy started.')

    ng = SIZE(atm_energy_list)

    DO jg = 1, ng
       !
       CALL vlr_del(atm_energy_list(jg))
       !
    END DO
    
    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(atm_energy)

    DEALLOCATE(atm_energy_config, STAT=ist)
    IF (ist/=success) CALL finish(thismodule, 'deallocation of atm_energy_config failed')

    DEALLOCATE(atm_energy_list, STAT=ist)
    IF (ist/=success) CALL finish(thismodule, 'deallocation of atm_energy_list failed')

    DEALLOCATE(atm_energy, STAT=ist )
    IF (ist/=success) CALL finish(thismodule, 'deallocation of atm_energy failed')

    CALL message(thismodule,'Destruction of atm_energy_config, atm_energy_list and atm_energy finished.')

  END SUBROUTINE destruct_atm_energy
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !>
  SUBROUTINE construct_atm_energy_list(jg,                  &
       &                               nproma, nlev, nblks, &
       &                               atm_energy_config,   &
       &                               atm_energy_list,     &
       &                               atm_energy)

    INTEGER             , INTENT(in)         :: jg                      !< grid index
    INTEGER             , INTENT(in)         :: nproma, nlev, nblks     !< size of dimensions

    TYPE(t_atm_energy_config), INTENT(inout) :: atm_energy_config       !< configuration switch
    TYPE(t_var_list_ptr), INTENT(inout)      :: atm_energy_list         !< list of variables with metadata
    TYPE(t_atm_energy)  , INTENT(inout)      :: atm_energy              !< contains pointers for variables

    ! Local variables

    CHARACTER(len= 2) :: cg
    CHARACTER(len=20) :: listname

    INTEGER           :: datatype_grb, datatype_flt, datatype_int

    INTEGER           :: shape0d(1), shape1d(1), shape2d(2), shape3d(3)

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    LOGICAL           :: lenergy                                        !< general switch

    LOGICAL           :: leindyn, leindynvi, leindynhi, leindynti, &
         &               leinphy, leinphyvi, leinphyhi, leinphyti, &
         &               lein,    leinvi,    leinhi,    leinti          !< logicals for internal energies

    LOGICAL           :: lekhdyn, lekhdynvi, lekhdynhi, lekhdynti, &
         &               lekhphy, lekhphyvi, lekhphyhi, lekhphyti, &
         &               lekh,    lekhvi,    lekhhi,    lekhti          !< logicals for horizontal kinetic energies

    LOGICAL           :: lekvdyn, lekvdynvi, lekvdynhi, lekvdynti, &
         &               lekvphy, lekvphyvi, lekvphyhi, lekvphyti, &
         &               lekv,    lekvvi,    lekvhi,    lekvti          !< logicals for vertical kinetic energies

    LOGICAL           :: legpdyn, legpdynvi, legpdynhi, legpdynti, &
         &               legpphy, legpphyvi, legpphyhi, legpphyti, &
         &               legp,    legpvi,    legphi,    legpti          !< logicals for geopotential energies

    LOGICAL           :: letodyn, letodynvi, letodynhi, letodynti, &
         &               letophy, letophyvi, letophyhi, letophyti, &
         &               leto,    letovi,    letohi,    letoti          !< logicals for total energies

    LOGICAL           :: ltend_eindyn, ltend_eindynvi, ltend_eindynhi, ltend_eindynti, &
         &               ltend_ekhdyn, ltend_ekhdynvi, ltend_ekhdynhi, ltend_ekhdynti, &
         &               ltend_ekvdyn, ltend_ekvdynvi, ltend_ekvdynhi, ltend_ekvdynti, &
         &               ltend_egpdyn, ltend_egpdynvi, ltend_egpdynhi, ltend_egpdynti, &
         &               ltend_etodyn, ltend_etodynvi, ltend_etodynhi, ltend_etodynti    !< logicals for dyn tendencies of energies

    LOGICAL           :: ltend_eincld, ltend_eincldvi, ltend_eincldhi, ltend_eincldti    !< logicals for cld tendencies of energies

    LOGICAL           :: ltend_einrad, ltend_einradvi, ltend_einradhi, ltend_einradti    !< logicals for rad tendencies of energies

    LOGICAL           :: ltend_eintmx, ltend_eintmxvi, ltend_eintmxhi, ltend_eintmxti, &
         &               ltend_ekhtmx, ltend_ekhtmxvi, ltend_ekhtmxhi, ltend_ekhtmxti, &
         &               ltend_ekvtmx, ltend_ekvtmxvi, ltend_ekvtmxhi, ltend_ekvtmxti    !< logicals for tmx tendencies of energies

    LOGICAL           :: ltend_einphy, ltend_einphyvi, ltend_einphyhi, ltend_einphyti, &
         &               ltend_ekhphy, ltend_ekhphyvi, ltend_ekhphyhi, ltend_ekhphyti, &
         &               ltend_ekvphy, ltend_ekvphyvi, ltend_ekvphyhi, ltend_ekvphyti, &
         &               ltend_egpphy, ltend_egpphyvi, ltend_egpphyhi, ltend_egpphyti, &
         &               ltend_etophy, ltend_etophyvi, ltend_etophyhi, ltend_etophyti    !< logicals for phy tendencies of energies


    ! Set logicals dependent on the list of output fields
    !
    lenergy   = .FALSE.
    !
    ! energies
    !
    leindyn   = is_variable_in_output(var_name='eindyn')
    leindynvi = is_variable_in_output(var_name='eindynvi')
    leindynhi = is_variable_in_output(var_name='eindynhi')
    leindynti = is_variable_in_output(var_name='eindynti')
    !
    leinphy   = is_variable_in_output(var_name='einphy')
    leinphyvi = is_variable_in_output(var_name='einphyvi')
    leinphyhi = is_variable_in_output(var_name='einphyhi')
    leinphyti = is_variable_in_output(var_name='einphyti')
    !
    lein      = is_variable_in_output(var_name='ein')
    leinvi    = is_variable_in_output(var_name='einvi')
    leinhi    = is_variable_in_output(var_name='einhi')
    leinti    = is_variable_in_output(var_name='einti')
    !
    lenergy   = lenergy .OR. leindyn .OR. leindynvi .OR. leindynhi .OR. leindynti &
         &              .OR. leinphy .OR. leinphyvi .OR. leinphyhi .OR. leinphyti &
         &              .OR. lein    .OR. leinvi    .OR. leinhi    .OR. leinti
    !
    !
    lekhdyn   = is_variable_in_output(var_name='ekhdyn')
    lekhdynvi = is_variable_in_output(var_name='ekhdynvi')
    lekhdynhi = is_variable_in_output(var_name='ekhdynhi')
    lekhdynti = is_variable_in_output(var_name='ekhdynti')
    !
    lekhphy   = is_variable_in_output(var_name='ekhphy')
    lekhphyvi = is_variable_in_output(var_name='ekhphyvi')
    lekhphyhi = is_variable_in_output(var_name='ekhphyhi')
    lekhphyti = is_variable_in_output(var_name='ekhphyti')
    !
    lekh      = is_variable_in_output(var_name='ekh')
    lekhvi    = is_variable_in_output(var_name='ekhvi')
    lekhhi    = is_variable_in_output(var_name='ekhhi')
    lekhti    = is_variable_in_output(var_name='ekhti')
    !
    lenergy   = lenergy .OR. lekhdyn .OR. lekhdynvi .OR. lekhdynhi .OR. lekhdynti &
         &              .OR. lekhphy .OR. lekhphyvi .OR. lekhphyhi .OR. lekhphyti &
         &              .OR. lekh    .OR. lekhvi    .OR. lekhhi    .OR. lekhti
    !
    !
    lekvdyn   = is_variable_in_output(var_name='ekvdyn')
    lekvdynvi = is_variable_in_output(var_name='ekvdynvi')
    lekvdynhi = is_variable_in_output(var_name='ekvdynhi')
    lekvdynti = is_variable_in_output(var_name='ekvdynti')
    !
    lekvphy   = is_variable_in_output(var_name='ekvphy')
    lekvphyvi = is_variable_in_output(var_name='ekvphyvi')
    lekvphyhi = is_variable_in_output(var_name='ekvphyhi')
    lekvphyti = is_variable_in_output(var_name='ekvphyti')
    !
    lekv      = is_variable_in_output(var_name='ekv')
    lekvvi    = is_variable_in_output(var_name='ekvvi')
    lekvhi    = is_variable_in_output(var_name='ekvhi')
    lekvti    = is_variable_in_output(var_name='ekvti')
    !
    lenergy   = lenergy .OR. lekvdyn .OR. lekvdynvi .OR. lekvdynhi .OR. lekvdynti &
         &              .OR. lekvphy .OR. lekvphyvi .OR. lekvphyhi .OR. lekvphyti &
         &              .OR. lekv    .OR. lekvvi    .OR. lekvhi    .OR. lekvti
    !
    !
    legpdyn   = is_variable_in_output(var_name='egpdyn')
    legpdynvi = is_variable_in_output(var_name='egpdynvi')
    legpdynhi = is_variable_in_output(var_name='egpdynhi')
    legpdynti = is_variable_in_output(var_name='egpdynti')
    !
    legpphy   = is_variable_in_output(var_name='egpphy')
    legpphyvi = is_variable_in_output(var_name='egpphyvi')
    legpphyhi = is_variable_in_output(var_name='egpphyhi')
    legpphyti = is_variable_in_output(var_name='egpphyti')
    !
    legp      = is_variable_in_output(var_name='egp')
    legpvi    = is_variable_in_output(var_name='egpvi')
    legphi    = is_variable_in_output(var_name='egphi')
    legpti    = is_variable_in_output(var_name='egpti')
    !
    lenergy   = lenergy .OR. legpdyn .OR. legpdynvi .OR. legpdynhi .OR. legpdynti &
         &              .OR. legpphy .OR. legpphyvi .OR. legpphyhi .OR. legpphyti &
         &              .OR. legp    .OR. legpvi    .OR. legphi    .OR. legpti
    !
    !
    letodyn   = is_variable_in_output(var_name='etodyn')
    letodynvi = is_variable_in_output(var_name='etodynvi')
    letodynhi = is_variable_in_output(var_name='etodynhi')
    letodynti = is_variable_in_output(var_name='etodynti')
    !
    letophy   = is_variable_in_output(var_name='etophy')
    letophyvi = is_variable_in_output(var_name='etophyvi')
    letophyhi = is_variable_in_output(var_name='etophyhi')
    letophyti = is_variable_in_output(var_name='etophyti')
    !
    leto      = is_variable_in_output(var_name='eto')
    letovi    = is_variable_in_output(var_name='etovi')
    letohi    = is_variable_in_output(var_name='etohi')
    letoti    = is_variable_in_output(var_name='etoti')
    !
    lenergy   = lenergy .OR. letodyn .OR. letodynvi .OR. letodynhi .OR. letodynti &
         &              .OR. letophy .OR. letophyvi .OR. letophyhi .OR. letophyti &
         &              .OR. leto    .OR. letovi    .OR. letohi    .OR. letoti
    !
    ! tendencies of energies
    !
    ltend_eindyn   = is_variable_in_output(var_name='tend_eindyn')
    ltend_eindynvi = is_variable_in_output(var_name='tend_eindynvi')
    ltend_eindynhi = is_variable_in_output(var_name='tend_eindynhi')
    ltend_eindynti = is_variable_in_output(var_name='tend_eindynti')
    !
    ltend_ekhdyn   = is_variable_in_output(var_name='tend_ekhdyn')
    ltend_ekhdynvi = is_variable_in_output(var_name='tend_ekhdynvi')
    ltend_ekhdynhi = is_variable_in_output(var_name='tend_ekhdynhi')
    ltend_ekhdynti = is_variable_in_output(var_name='tend_ekhdynti')
    !
    ltend_ekvdyn   = is_variable_in_output(var_name='tend_ekvdyn')
    ltend_ekvdynvi = is_variable_in_output(var_name='tend_ekvdynvi')
    ltend_ekvdynhi = is_variable_in_output(var_name='tend_ekvdynhi')
    ltend_ekvdynti = is_variable_in_output(var_name='tend_ekvdynti')
    !
    ltend_egpdyn   = is_variable_in_output(var_name='tend_egpdyn')
    ltend_egpdynvi = is_variable_in_output(var_name='tend_egpdynvi')
    ltend_egpdynhi = is_variable_in_output(var_name='tend_egpdynhi')
    ltend_egpdynti = is_variable_in_output(var_name='tend_egpdynti')
    !
    ltend_etodyn   = is_variable_in_output(var_name='tend_etodyn')
    ltend_etodynvi = is_variable_in_output(var_name='tend_etodynvi')
    ltend_etodynhi = is_variable_in_output(var_name='tend_etodynhi')
    ltend_etodynti = is_variable_in_output(var_name='tend_etodynti')
    !
    lenergy        = lenergy .OR. ltend_eindyn .OR. ltend_eindynvi .OR. ltend_eindynhi .OR. ltend_eindynti &
         &                   .OR. ltend_ekhdyn .OR. ltend_ekhdynvi .OR. ltend_ekhdynhi .OR. ltend_ekhdynti &
         &                   .OR. ltend_ekvdyn .OR. ltend_ekvdynvi .OR. ltend_ekvdynhi .OR. ltend_ekvdynti &
         &                   .OR. ltend_egpdyn .OR. ltend_egpdynvi .OR. ltend_egpdynhi .OR. ltend_egpdynti &
         &                   .OR. ltend_etodyn .OR. ltend_etodynvi .OR. ltend_etodynhi .OR. ltend_etodynti
    !
    !
    ltend_eincld   = is_variable_in_output(var_name='tend_eincld')
    ltend_eincldvi = is_variable_in_output(var_name='tend_eincldvi')
    ltend_eincldhi = is_variable_in_output(var_name='tend_eincldhi')
    ltend_eincldti = is_variable_in_output(var_name='tend_eincldti')
    !
    lenergy        = lenergy .OR. ltend_eincld .OR. ltend_eincldvi .OR. ltend_eincldhi .OR. ltend_eincldti
    !
    !
    ltend_einrad   = is_variable_in_output(var_name='tend_einrad')
    ltend_einradvi = is_variable_in_output(var_name='tend_einradvi')
    ltend_einradhi = is_variable_in_output(var_name='tend_einradhi')
    ltend_einradti = is_variable_in_output(var_name='tend_einradti')
    !
    lenergy        = lenergy .OR. ltend_einrad .OR. ltend_einradvi .OR. ltend_einradhi .OR. ltend_einradti
    !
    !
    ltend_eintmx   = is_variable_in_output(var_name='tend_eintmx')
    ltend_eintmxvi = is_variable_in_output(var_name='tend_eintmxvi')
    ltend_eintmxhi = is_variable_in_output(var_name='tend_eintmxhi')
    ltend_eintmxti = is_variable_in_output(var_name='tend_eintmxti')
    !
    ltend_ekhtmx   = is_variable_in_output(var_name='tend_ekhtmx')
    ltend_ekhtmxvi = is_variable_in_output(var_name='tend_ekhtmxvi')
    ltend_ekhtmxhi = is_variable_in_output(var_name='tend_ekhtmxhi')
    ltend_ekhtmxti = is_variable_in_output(var_name='tend_ekhtmxti')
    !
    ltend_ekvtmx   = is_variable_in_output(var_name='tend_ekvtmx')
    ltend_ekvtmxvi = is_variable_in_output(var_name='tend_ekvtmxvi')
    ltend_ekvtmxhi = is_variable_in_output(var_name='tend_ekvtmxhi')
    ltend_ekvtmxti = is_variable_in_output(var_name='tend_ekvtmxti')
    !
    lenergy        = lenergy .OR. ltend_eintmx .OR. ltend_eintmxvi .OR. ltend_eintmxhi .OR. ltend_eintmxti &
         &                   .OR. ltend_ekhtmx .OR. ltend_ekhtmxvi .OR. ltend_ekhtmxhi .OR. ltend_ekhtmxti &
         &                   .OR. ltend_ekvtmx .OR. ltend_ekvtmxvi .OR. ltend_ekvtmxhi .OR. ltend_ekvtmxti
    !
    !
    ltend_einphy   = is_variable_in_output(var_name='tend_einphy')
    ltend_einphyvi = is_variable_in_output(var_name='tend_einphyvi')
    ltend_einphyhi = is_variable_in_output(var_name='tend_einphyhi')
    ltend_einphyti = is_variable_in_output(var_name='tend_einphyti')
    !
    ltend_ekhphy   = is_variable_in_output(var_name='tend_ekhphy')
    ltend_ekhphyvi = is_variable_in_output(var_name='tend_ekhphyvi')
    ltend_ekhphyhi = is_variable_in_output(var_name='tend_ekhphyhi')
    ltend_ekhphyti = is_variable_in_output(var_name='tend_ekhphyti')
    !
    ltend_ekvphy   = is_variable_in_output(var_name='tend_ekvphy')
    ltend_ekvphyvi = is_variable_in_output(var_name='tend_ekvphyvi')
    ltend_ekvphyhi = is_variable_in_output(var_name='tend_ekvphyhi')
    ltend_ekvphyti = is_variable_in_output(var_name='tend_ekvphyti')
    !
    ltend_egpphy   = is_variable_in_output(var_name='tend_egpphy')
    ltend_egpphyvi = is_variable_in_output(var_name='tend_egpphyvi')
    ltend_egpphyhi = is_variable_in_output(var_name='tend_egpphyhi')
    ltend_egpphyti = is_variable_in_output(var_name='tend_egpphyti')
    !
    ltend_etophy   = is_variable_in_output(var_name='tend_etophy')
    ltend_etophyvi = is_variable_in_output(var_name='tend_etophyvi')
    ltend_etophyhi = is_variable_in_output(var_name='tend_etophyhi')
    ltend_etophyti = is_variable_in_output(var_name='tend_etophyti')
    !
    lenergy        = lenergy .OR. ltend_einphy .OR. ltend_einphyvi .OR. ltend_einphyhi .OR. ltend_einphyti &
         &                   .OR. ltend_ekhphy .OR. ltend_ekhphyvi .OR. ltend_ekhphyhi .OR. ltend_ekhphyti &
         &                   .OR. ltend_ekvphy .OR. ltend_ekvphyvi .OR. ltend_ekvphyhi .OR. ltend_ekvphyti &
         &                   .OR. ltend_egpphy .OR. ltend_egpphyvi .OR. ltend_egpphyhi .OR. ltend_egpphyti &
         &                   .OR. ltend_etophy .OR. ltend_etophyvi .OR. ltend_etophyhi .OR. ltend_etophyti

    atm_energy_config%l_atm_energy = lenergy

    WRITE(cg,'(i2.2)') jg
    CALL message('construct_atm_energy_list','create list and allocate memory for jg ='//cg)

    ! number of bits for data representation in grib2
    datatype_grb = MERGE(datatype_pack24, datatype_pack16, gribout_config(jg)%lgribout_24bit)

    ! number of bits for data representation in netcdf
    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
    datatype_int = datatype_int32

    ! array shapes for fields of different dimensionality
    shape0d  = [           1       ] !< a scalar represented as a single value in 1 layer
    shape1d  = [        nlev       ] !< a vertical profile of nlev levels
    shape2d  = [nproma,       nblks] !< a horizontal field stored in an nproma x nblks array
    shape3d  = [nproma, nlev, nblks] !< a 3-dimensional field

    ! define list name
    !
    WRITE(listname,'(a,i2.2)') 'atm_energy_list_D',jg

    ! register atm_energy_list for grid jg
    !
    CALL vlr_add(atm_energy_list            ,&
         &       listname                   ,&
         &       patch_id  = jg             ,&
         &       loutput   = .TRUE.         ,&
         &       lrestart  = .FALSE.        ,&
         &       linitial  = .FALSE.        ,&
         &       model_type=get_my_process_name())

    
    ! add variables to the list


    ! atmosphere energy state
    ! =======================
    !
    ! - internal energy
    !
    !   - 1st stage, for difference ..2.. - ..1.., and for output after the last physics parameterization
    !
    IF (   leinphy   .OR. ltend_eindyn   .OR. ltend_eincld   .OR. ltend_einrad   .OR. ltend_eintmx .OR. &
         & leinphyhi .OR. ltend_eindynhi .OR. ltend_eincldhi .OR. ltend_einradhi .OR. ltend_eintmxhi) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy_density_after_params', 'J m-3', &
            &                'Atmosphere Moist Internal Energy Density After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'einphy', atm_energy%ein1,                  &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein1)
    END IF
    !
    IF (   leinphyvi .OR. ltend_eindynvi .OR. ltend_eincldvi .OR. ltend_einradvi .OR. ltend_eintmxvi .OR. &
         & leinphyti .OR. ltend_eindynti .OR. ltend_eincldti .OR. ltend_einradti .OR. ltend_eintmxti) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy_content_after_params', 'J m-2', &
            &                'Atmosphere Moist Internal Energy Content After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'einphyvi', atm_energy%ein1vi,              &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein1vi)
    END IF
    !
    IF (leinphyhi .OR. ltend_eindynhi .OR. ltend_eincldhi .OR. ltend_einradhi .OR. ltend_eintmxhi) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy_profile_after_params', 'J m-1', &
            &                'Atmosphere Moist Internal Energy Profile After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'einphyhi', atm_energy%ein1hi,              &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein1hi)
    END IF
    !
    IF (leinphyti .OR. ltend_eindynti .OR. ltend_eincldti .OR. ltend_einradti .OR. ltend_eintmxti) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy_after_params', 'J', &
            &                'Atmosphere Moist Internal Energy After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'einphyti', atm_energy%ein1ti,              &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein1ti)
    END IF
    !
    !   - 2nd stage, for difference ..2.. - ..1.. and ..2.. - ..3.., and for output at the end of the time step
    !
    IF (   lein   .OR. leindyn   .OR. leinphy   .OR. &
         & leinhi .OR. leindynhi .OR. leinphyhi .OR. &
         & ltend_eindyn   .OR. ltend_eincld   .OR. ltend_einrad   .OR. ltend_eintmx   .OR. ltend_einphy .OR. &
         & ltend_eindynhi .OR. ltend_eincldhi .OR. ltend_einradhi .OR. ltend_eintmxhi .OR. ltend_einphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy_density', 'J m-3', &
            &                'Atmosphere Moist Internal Energy Density', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ein', atm_energy%ein2,                     &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein2)
    END IF
    !
    IF (   leinvi .OR. leindynvi .OR. leinphyvi .OR. &
         & leinti .OR. leindynti .OR. leinphyti .OR. &
         & ltend_eindynvi .OR. ltend_eincldvi .OR. ltend_einradvi .OR. ltend_eintmxvi .OR. ltend_einphyvi .OR. &
         & ltend_eindynti .OR. ltend_eincldti .OR. ltend_einradti .OR. ltend_eintmxti .OR. ltend_einphyti) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy_content', 'J m-2', &
            &                'Atmosphere Moist Internal Energy Content', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'einvi', atm_energy%ein2vi,                 &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein2vi)
    END IF
    !
    IF (   leinhi .OR. leindynhi .OR. leinphyhi .OR. &
         & ltend_eindynhi .OR. ltend_eincldhi .OR. ltend_einradhi .OR. ltend_eintmxhi .OR. ltend_einphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy_profile', 'J m-1', &
            &                'Atmosphere Moist Internal Energy Profile', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'einhi', atm_energy%ein2hi,                 &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein2hi)
    END IF
    !
    IF (   leinti .OR. leindynti .OR. leinphyti .OR. &
         & ltend_eindynti .OR. ltend_eincldti .OR. ltend_einradti .OR. ltend_eintmxti .OR. ltend_einphyti) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy', 'J', &
            &                'Atmosphere Moist Internal Energy', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'einti', atm_energy%ein2ti,                 &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein2ti)
    END IF
    !
    !   - 3rd stage, for difference ..2.. - ..3.., and for output after dynamics
    !
    IF (   leindyn   .OR. ltend_einphy .OR. &
         & leindynhi .OR. ltend_einphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy_density_after_dynamics', 'J m-3', &
            &                'Atmosphere Moist Internal Energy Density After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'eindyn', atm_energy%ein3,                  &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein3)
    END IF
    !
    IF (   leindynvi .OR. ltend_einphyvi .OR. &
         & leindynti .OR. ltend_einphyti) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy_content_after_dynamics', 'J m-2', &
            &                'Atmosphere Moist Internal Energy Content After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'eindynvi', atm_energy%ein3vi,              &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein3vi)
    END IF
    !
    IF (leindynhi .OR. ltend_einphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy_profile_after_dynamics', 'J m-1', &
            &                'Atmosphere Moist Internal Energy Profile After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'eindynhi', atm_energy%ein3hi,              &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein3hi)
    END IF
    !
    IF (leindynti .OR. ltend_einphyti) THEN
       cf_desc    = t_cf_var('atmosphere_moist_internal_energy_after_dynamics', 'J', &
            &                'Atmosphere Moist Internal Energy After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'eindynti', atm_energy%ein3ti,              &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ein3ti)
    END IF
    !
    !
    ! - horizontal kinetic energy
    !
    !   - 1st stage, for difference ..2.. - ..1.., and for output after the last physics parameterization
    !
    IF (   lekhphy   .OR. ltend_ekhdyn   .OR. ltend_ekhtmx .OR. &
         & lekhphyhi .OR. ltend_ekhdynhi .OR. ltend_ekhtmxhi) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy_density_after_params', 'J m-3', &
            &                'Atmosphere Horizontal Kinetic Energy Density After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekhphy', atm_energy%ekh1,                  &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh1)
    END IF
    !
    IF (   lekhphyvi .OR. ltend_ekhdynvi .OR. ltend_ekhtmxvi .OR. &
         & lekhphyti .OR. ltend_ekhdynti .OR. ltend_ekhtmxti) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy_content_after_params', 'J m-2', &
            &                'Atmosphere Horizontal Kinetic Energy Content After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekhphyvi', atm_energy%ekh1vi,              &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh1vi)
    END IF
    !
    IF (lekhphyhi .OR. ltend_ekhdynhi .OR. ltend_ekhtmxhi) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy_profile_after_params', 'J m-1', &
            &                'Atmosphere Horizontal Kinetic Energy Profile After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekhphyhi', atm_energy%ekh1hi,              &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh1hi)
    END IF
    !
    IF (lekhphyti .OR. ltend_ekhdynti .OR. ltend_ekhtmxti) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy_after_params', 'J', &
            &                'Atmosphere Horizontal Kinetic Energy After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekhphyti', atm_energy%ekh1ti,              &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh1ti)
    END IF
    !
    !   - 2nd stage, for difference ..2.. - ..1.. and ..2.. - ..3.., and for output at the end of the time step
    !
    IF (   lekh   .OR. lekhdyn   .OR. lekhphy   .OR. &
         & lekhhi .OR. lekhdynhi .OR. lekhphyhi .OR. &
         & ltend_ekhdyn   .OR. ltend_ekhtmx   .OR. ltend_ekhphy .OR. &
         & ltend_ekhdynhi .OR. ltend_ekhtmxhi .OR. ltend_ekhphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy_density', 'J m-3', &
            &                'Atmosphere Horizontal Kinetic Energy Density', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekh', atm_energy%ekh2,                     &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh2)
    END IF
    !
    IF (   lekhvi .OR. lekhdynvi .OR. lekhphyvi .OR. &
         & lekhti .OR. lekhdynti .OR. lekhphyti .OR. &
         & ltend_ekhdynvi .OR. ltend_ekhtmxvi .OR. ltend_ekhphyvi .OR. &
         & ltend_ekhdynti .OR. ltend_ekhtmxti .OR. ltend_ekhphyti) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy_content', 'J m-2', &
            &                'Atmosphere Horizontal Kinetic Energy Content', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekhvi', atm_energy%ekh2vi,                 &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh2vi)
    END IF
    !
    IF (   lekhhi .OR. lekhdynhi .OR. lekhphyhi .OR. &
         & ltend_ekhdynhi .OR. ltend_ekhtmxhi .OR. ltend_ekhphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy_profile', 'J m-1', &
            &                'Atmosphere Horizontal Kinetic Energy Profile', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekhhi', atm_energy%ekh2hi,                 &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh2hi)
    END IF
    !
    IF (   lekhti .OR. lekhdynti .OR. lekhphyti .OR. &
         & ltend_ekhdynti .OR. ltend_ekhtmxti .OR. ltend_ekhphyti) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy', 'J', &
            &                'Atmosphere Horizontal Kinetic Energy', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekhti', atm_energy%ekh2ti,                 &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh2ti)
    END IF
    !
    !   - 3rd stage, for difference ..2.. - ..3.., and for output after dynamics
    !
    IF (   lekhdyn   .OR. ltend_ekhphy .OR. &
         & lekhdynhi .OR. ltend_ekhphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy_density_after_dynamics', 'J m-3', &
            &                'Atmosphere Horizontal Kinetic Energy Density After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekhdyn', atm_energy%ekh3,                  &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh3)
    END IF
    !
    IF (   lekhdynvi .OR. ltend_ekhphyvi .OR. &
         & lekhdynti .OR. ltend_ekhphyti) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy_content_after_dynamics', 'J m-2', &
            &                'Atmosphere Horizontal Kinetic Energy Content After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekhdynvi', atm_energy%ekh3vi,              &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh3vi)
    END IF
    !
    IF (lekhdynhi .OR. ltend_ekhphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy_profile_after_dynamics', 'J m-1', &
            &                'Atmosphere Horizontal Kinetic Energy Profile After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekhdynhi', atm_energy%ekh3hi,              &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh3hi)
    END IF
    !
    IF (lekhdynti .OR. ltend_ekhphyti) THEN
       cf_desc    = t_cf_var('atmosphere_horizontal_kinetic_energy_after_dynamics', 'J', &
            &                'Atmosphere Horizontal Kinetic Energy After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekhdynti', atm_energy%ekh3ti,              &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekh3ti)
    END IF
    !
    !
    ! - vertical kinetic energy
    !
    !   - 1st stage, for difference ..2.. - ..1.., and for output after the last physics parameterization
    !
    IF (   lekvphy   .OR. ltend_ekvdyn   .OR. ltend_ekvtmx .OR. &
         & lekvphyhi .OR. ltend_ekvdynhi .OR. ltend_ekvtmxhi) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy_density_after_params', 'J m-3', &
            &                'Atmosphere Vertical Kinetic Energy Density After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekvphy', atm_energy%ekv1,                  &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv1)
    END IF
    !
    IF (   lekvphyvi .OR. ltend_ekvdynvi .OR. ltend_ekvtmxvi .OR. &
         & lekvphyti .OR. ltend_ekvdynti .OR. ltend_ekvtmxti) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy_content_after_params', 'J m-2', &
            &                'Atmosphere Vertical Kinetic Energy Content After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekvphyvi', atm_energy%ekv1vi,              &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv1vi)
    END IF
    !
    IF (lekvphyhi .OR. ltend_ekvdynhi .OR. ltend_ekvtmxhi) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy_profile_after_params', 'J m-1', &
            &                'Atmosphere Vertical Kinetic Energy Profile After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekvphyhi', atm_energy%ekv1hi,              &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv1hi)
    END IF
    !
    IF (lekvphyti .OR. ltend_ekvdynti .OR. ltend_ekvtmxti) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy_after_params', 'J', &
            &                'Atmosphere Vertical Kinetic Energy After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekvphyti', atm_energy%ekv1ti,              &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv1ti)
    END IF
    !
    !   - 2nd stage, for difference ..2.. - ..1.. and ..2.. - ..3.., and for output at the end of the time step
    !
    IF (   lekv   .OR. lekvdyn   .OR. lekvphy   .OR. &
         & lekvhi .OR. lekvdynhi .OR. lekvphyhi .OR. &
         & ltend_ekvdyn   .OR. ltend_ekvtmx   .OR. ltend_ekvphy .OR. &
         & ltend_ekvdynhi .OR. ltend_ekvtmxhi .OR. ltend_ekvphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy_density', 'J m-3', &
            &                'Atmosphere Vertical Kinetic Energy Density', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekv', atm_energy%ekv2,                     &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv2)
    END IF
    !
    IF (   lekvvi .OR. lekvdynvi .OR. lekvphyvi .OR. &
         & lekvti .OR. lekvdynti .OR. lekvphyti .OR. &
         & ltend_ekvdynvi .OR. ltend_ekvtmxvi .OR. ltend_ekvphyvi .OR. &
         & ltend_ekvdynti .OR. ltend_ekvtmxti .OR. ltend_ekvphyti) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy_content', 'J m-2', &
            &                'Atmosphere Vertical Kinetic Energy Content', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekvvi', atm_energy%ekv2vi,                 &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv2vi)
    END IF
    !
    IF (   lekvhi .OR. lekvdynhi .OR. lekvphyhi .OR. &
         & ltend_ekvdynhi .OR. ltend_ekvtmxhi .OR. ltend_ekvphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy_profile', 'J m-1', &
            &                'Atmosphere Vertical Kinetic Energy Profile', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekvhi', atm_energy%ekv2hi,                 &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv2hi)
    END IF
    !
    IF (   lekvti .OR. lekvdynti .OR. lekvphyti .OR. &
         & ltend_ekvdynti .OR. ltend_ekvtmxti .OR. ltend_ekvphyti) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy', 'J', &
            &                'Atmosphere Vertical Kinetic Energy', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekvti', atm_energy%ekv2ti,                 &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv2ti)
    END IF
    !
    !   - 3rd stage, for difference ..2.. - ..3.., and for output after dynamics
    !
    IF (   lekvdyn   .OR. ltend_ekvphy .OR. &
         & lekvdynhi .OR. ltend_ekvphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy_density_after_dynamics', 'J m-3', &
            &                'Atmosphere Vertical Kinetic Energy Density After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekvdyn', atm_energy%ekv3,                  &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv3)
    END IF
    !
    IF (   lekvdynvi .OR. ltend_ekvphyvi .OR. &
         & lekvdynti .OR. ltend_ekvphyti) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy_content_after_dynamics', 'J m-2', &
            &                'Atmosphere Vertical Kinetic Energy Content After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'ekvdynvi', atm_energy%ekv3vi,              &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv3vi)
    END IF
    !
    IF (lekvdynhi .OR. ltend_ekvphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy_profile_after_dynamics', 'J m-1', &
            &                'Atmosphere Vertical Kinetic Energy Profile After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekvdynhi', atm_energy%ekv3hi,              &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv3hi)
    END IF
    !
    IF (lekvdynti .OR. ltend_ekvphyti) THEN
       cf_desc    = t_cf_var('atmosphere_vertical_kinetic_energy_after_dynamics', 'J', &
            &                'Atmosphere Vertical Kinetic Energy After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'ekvdynti', atm_energy%ekv3ti,              &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekv3ti)
    END IF
    !
    !
    ! - geopotential energy
    !
    !   - 1st stage, for difference ..2.. - ..1.., and for output after the last physics parameterization
    !
    IF (   legpphy   .OR. ltend_egpdyn .OR. &
         & legpphyhi .OR. ltend_egpdynhi) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy_density_after_params', 'J m-3', &
            &                'Atmosphere Geopotential Energy Density After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'egpphy', atm_energy%egp1,                  &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp1)
    END IF
    !
    IF (   legpphyvi .OR. ltend_egpdynvi .OR. &
         & legpphyti .OR. ltend_egpdynti) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy_content_after_params', 'J m-2', &
            &                'Atmosphere Geopotential Energy Content After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'egpphyvi', atm_energy%egp1vi,              &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp1vi)
    END IF
    !
    IF (legpphyhi .OR. ltend_egpdynhi) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy_profile_after_params', 'J m-1', &
            &                'Atmosphere Geopotential Energy Profile After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'egpphyhi', atm_energy%egp1hi,              &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp1hi)
    END IF
    !
    IF (legpphyti .OR. ltend_egpdynti) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy_after_params', 'J', &
            &                'Atmosphere Geopotential Energy After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'egpphyti', atm_energy%egp1ti,              &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp1ti)
    END IF
    !
    !   - 2nd stage, for difference ..2.. - ..1.. and ..2.. - ..3.., and for output at the end of the time step
    !
    IF (   legp   .OR. legpdyn   .OR. legpphy   .OR. &
         & legphi .OR. legpdynhi .OR. legpphyhi .OR. &
         & ltend_egpdyn   .OR. ltend_egpphy .OR. &
         & ltend_egpdynhi .OR. ltend_egpphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy_density', 'J m-3', &
            &                'Atmosphere Geopotential Energy Density', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'egp', atm_energy%egp2,                     &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp2)
    END IF
    !
    IF (   legpvi .OR. legpdynvi .OR. legpphyvi .OR. &
         & legpti .OR. legpdynti .OR. legpphyti .OR. &
         & ltend_egpdynvi .OR. ltend_egpphyvi .OR. &
         & ltend_egpdynti .OR. ltend_egpphyti) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy_content', 'J m-2', &
            &                'Atmosphere Geopotential Energy Content', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'egpvi', atm_energy%egp2vi,                 &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp2vi)
    END IF
    !
    IF (   legphi .OR. legpdynhi .OR. legpphyhi .OR. &
         & ltend_egpdynhi .OR. ltend_egpphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy_profile', 'J m-1', &
            &                'Atmosphere Geopotential Energy Profile', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'egphi', atm_energy%egp2hi,                 &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp2hi)
    END IF
    !
    IF (   legpti .OR. legpdynti .OR. legpphyti .OR. &
         & ltend_egpdynti .OR. ltend_egpphyti) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy', 'J', &
            &                'Atmosphere Geopotential Energy', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'egpti', atm_energy%egp2ti,                 &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp2ti)
    END IF
    !
    !   - 3rd stage, for difference ..2.. - ..3.., and for output after dynamics
    !
    IF (   legpdyn   .OR. ltend_egpphy .OR. &
         & legpdynhi .OR. ltend_egpphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy_density_after_dynamics', 'J m-3', &
            &                'Atmosphere Geopotential Energy Density After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'egpdyn', atm_energy%egp3,                  &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp3)
    END IF
    !
    IF (   legpdynvi .OR. ltend_egpphyvi .OR. &
         & legpdynti .OR. ltend_egpphyti) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy_content_after_dynamics', 'J m-2', &
            &                'Atmosphere Geopotential Energy Content After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'egpdynvi', atm_energy%egp3vi,              &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp3vi)
    END IF
    !
    IF (legpdynhi .OR. ltend_egpphyhi) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy_profile_after_dynamics', 'J m-1', &
            &                'Atmosphere Geopotential Energy Profile After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'egpdynhi', atm_energy%egp3hi,              &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp3hi)
    END IF
    !
    IF (legpdynti .OR. ltend_egpphyti) THEN
       cf_desc    = t_cf_var('atmosphere_geopotential_energy_after_dynamics', 'J', &
            &                'Atmosphere Geopotential Energy After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'egpdynti', atm_energy%egp3ti,              &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%egp3ti)
    END IF
    !
    !
    ! - total energy
    !
    !   - 1st stage, for difference ..2.. - ..1.., and for output after the last physics parameterization
    !
    IF (   letophy   .OR. ltend_etodyn .OR. &
         & letophyhi .OR. ltend_etodynhi) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy_density_after_params', 'J m-3', &
            &                'Atmosphere Total Energy Density After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'etophy', atm_energy%eto1,                  &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto1)
    END IF
    !
    IF (   letophyvi .OR. ltend_etodynvi .OR. &
         & letophyti .OR. ltend_etodynti) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy_content_after_params', 'J m-2', &
            &                'Atmosphere Total Energy Content After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'etophyvi', atm_energy%eto1vi,              &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto1vi)
    END IF
    !
    IF (letophyhi .OR. ltend_etodynhi) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy_profile_after_params', 'J m-1', &
            &                'Atmosphere Total Energy Profile After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'etophyhi', atm_energy%eto1hi,              &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto1hi)
    END IF
    !
    IF (letophyti .OR. ltend_etodynti) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy_after_params', 'J', &
            &                'Atmosphere Total Energy After Params', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'etophyti', atm_energy%eto1ti,              &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto1ti)
    END IF
    !
    !   - 2nd stage, for difference ..2.. - ..1.. and ..2.. - ..3.., and for output at the end of the time step
    !
    IF (   leto   .OR. letodyn   .OR. letophy   .OR. &
         & letohi .OR. letodynhi .OR. letophyhi .OR. &
         & ltend_etodyn   .OR. ltend_etophy .OR. &
         & ltend_etodynhi .OR. ltend_etophyhi) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy_density', 'J m-3', &
            &                'Atmosphere Total Energy Density', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'eto', atm_energy%eto2,                     &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto2)
    END IF
    !
    IF (   letovi .OR. letodynvi .OR. letophyvi .OR. &
         & letoti .OR. letodynti .OR. letophyti .OR. &
         & ltend_etodynvi .OR. ltend_etophyvi .OR. &
         & ltend_etodynti .OR. ltend_etophyti) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy_content', 'J m-2', &
            &                'Atmosphere Total Energy Content', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'etovi', atm_energy%eto2vi,                 &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto2vi)
    END IF
    !
    IF (   letohi .OR. letodynhi .OR. letophyhi .OR. &
         & ltend_etodynhi .OR. ltend_etophyhi) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy_profile', 'J m-1', &
            &                'Atmosphere Total Energy Profile', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'etohi', atm_energy%eto2hi,                 &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto2hi)
    END IF
    !
    IF (   letoti .OR. letodynti .OR. letophyti .OR. &
         & ltend_etodynti .OR. ltend_etophyti) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy', 'J', &
            &                'Atmosphere Total Energy', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'etoti', atm_energy%eto2ti,                 &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto2ti)
    END IF
    !
    !   - 3rd stage, for difference ..2.. - ..3.., and for output after dynamics
    !
    IF (   letodyn   .OR. ltend_etophy .OR. &
         & letodynhi .OR. ltend_etophyhi) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy_density_after_dynamics', 'J m-3', &
            &                'Atmosphere Total Energy Density After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'etodyn', atm_energy%eto3,                  &
            &        grid_unstructured_cell, za_reference,                        &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape3d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto3)
    END IF
    !
    IF (   letodynvi .OR. ltend_etophyvi .OR. &
         & letodynti .OR. ltend_etophyti) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy_content_after_dynamics', 'J m-2', &
            &                'Atmosphere Total Energy Content After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_cell)
       CALL add_var( atm_energy_list, 'etodynvi', atm_energy%eto3vi,              &
            &        grid_unstructured_cell, za_surface,                          &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape2d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto3vi)
    END IF
    !
    IF (letodynhi .OR. ltend_etophyhi) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy_profile_after_dynamics', 'J m-1', &
            &                'Atmosphere Total Energy Profile After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'etodynhi', atm_energy%eto3hi,              &
            &        grid_lonlat, za_reference,                                   &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape1d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto3hi)
    END IF
    !
    IF (letodynti .OR. ltend_etophyti) THEN
       cf_desc    = t_cf_var('atmosphere_total_energy_after_dynamics', 'J', &
            &                'Atmosphere Total Energy After Dynamics', datatype_flt)
       grib2_desc = grib2_var(255,255,255, datatype_grb, grid_unstructured, grid_lonlat)
       CALL add_var( atm_energy_list, 'etodynti', atm_energy%eto3ti,              &
            &        grid_lonlat, za_surface,                                     &
            &        cf_desc, grib2_desc,                                         &
            &        ldims=shape0d,                                               &
            &        lrestart = .FALSE.,                                          &
            &        isteptype=tstep_instant,                                     &
            &        lopenacc=.TRUE.)
       __acc_attach(atm_energy%eto3ti)
    END IF
    !

    ! atmosphere energy tendency
    ! ==========================
    !
    ! dyn: due to dynamics
    ! --------------------

    ! - of internal energy 
    !
    IF (ltend_eindyn) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_density_due_to_dynamics', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Density due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_eindyn', atm_energy%eindyn,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eindyn)
    END IF
    !
    !
    IF (ltend_eindynvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_content_due_to_dynamics', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Content due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_eindynvi', atm_energy%eindynvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eindynvi)
    END IF
    !
    !
    IF (ltend_eindynhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_profile_due_to_dynamics', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Profile due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_eindynhi', atm_energy%eindynhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eindynhi)
    END IF
    !
    !
    IF (ltend_eindynti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_due_to_dynamics', 'J s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_eindynti', atm_energy%eindynti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eindynti)
    END IF


    ! - of horizontal kinetic energy
    !
    IF (ltend_ekhdyn) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_density_due_to_dynamics', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy Density due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekhdyn', atm_energy%ekhdyn,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhdyn)
    END IF
    !
    !
    IF (ltend_ekhdynvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_content_due_to_dynamics', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy Content due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekhdynvi', atm_energy%ekhdynvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhdynvi)
    END IF
    !
    !
    IF (ltend_ekhdynhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_profile_due_to_dynamics', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy Profile due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekhdynhi', atm_energy%ekhdynhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhdynhi)
    END IF
    !
    !
    IF (ltend_ekhdynti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_due_to_dynamics', 'J s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekhdynti', atm_energy%ekhdynti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhdynti)
    END IF


    ! - of vertical kinetic energy
    !
    IF (ltend_ekvdyn) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_density_due_to_dynamics', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy Density due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekvdyn', atm_energy%ekvdyn,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvdyn)
    END IF
    !
    !
    IF (ltend_ekvdynvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_content_due_to_dynamics', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy Content due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekvdynvi', atm_energy%ekvdynvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvdynvi)
    END IF
    !
    !
    IF (ltend_ekvdynhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_profile_due_to_dynamics', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy Profile due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekvdynhi', atm_energy%ekvdynhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvdynhi)
    END IF
    !
    !
    IF (ltend_ekvdynti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_due_to_dynamics', 'J s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekvdynti', atm_energy%ekvdynti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvdynti)
    END IF


    ! of geopotential energy
    !
    IF (ltend_egpdyn) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_geopotential_energy_density_due_to_dynamics', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Geopotential Energy Density due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_egpdyn', atm_energy%egpdyn,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%egpdyn)
    END IF
    !
    !
    IF (ltend_egpdynvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_geopotential_energy_content_due_to_dynamics', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Geopotential Energy Content due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_egpdynvi', atm_energy%egpdynvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%egpdynvi)
    END IF
    !
    !
    IF (ltend_egpdynhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_geopotential_energy_profile_due_to_dynamics', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Geopotential Energy Profile due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_egpdynhi', atm_energy%egpdynhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%egpdynhi)
    END IF
    !
    !
    IF (ltend_egpdynti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_geopotential_energy_due_to_dynamics', 'J s-1', &
            &             'Tendency of Atmosphere Geopotential Energy due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_egpdynti', atm_energy%egpdynti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%egpdynti)
    END IF


    ! of total energy
    !
    IF (ltend_etodyn) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_total_energy_density_due_to_dynamics', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Total Energy Density due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_etodyn', atm_energy%etodyn,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%etodyn)
    END IF
    !
    !
    IF (ltend_etodynvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_total_energy_content_due_to_dynamics', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Total Energy Content due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_etodynvi', atm_energy%etodynvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%etodynvi)
    END IF
    !
    !
    IF (ltend_etodynhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_total_energy_profile_due_to_dynamics', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Total Energy Profile due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_etodynhi', atm_energy%etodynhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%etodynhi)
    END IF
    !
    !
    IF (ltend_etodynti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_total_energy_due_to_dynamics', 'J s-1', &
            &             'Tendency of Atmosphere Total Energy due to Dynamics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_etodynti', atm_energy%etodynti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%etodynti)
    END IF


    ! cld: due to cloud microphysics
    ! ------------------------------

    ! - of internal energy 
    !
    IF (ltend_eincld) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_density_due_to_cloud_microphysics', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Density due to Cloud Microphysics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_eincld', atm_energy%eincld,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eincld)
    END IF
    !
    !
    IF (ltend_eincldvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_content_due_to_cloud_microphysics', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Content due to Cloud Microphysics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_eincldvi', atm_energy%eincldvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eincldvi)
    END IF
    !
    !
    IF (ltend_eincldhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_profile_due_to_cloud_microphysics', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Profile due to Cloud Microphysics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_eincldhi', atm_energy%eincldhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eincldhi)
    END IF
    !
    !
    IF (ltend_eincldti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_due_to_cloud_microphysics', 'J s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy due to Cloud Microphysics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_eincldti', atm_energy%eincldti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eincldti)
    END IF


    ! rad: due to radiation
    ! ---------------------

    ! - of internal energy 
    !
    IF (ltend_einrad) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_density_due_to_radiation', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Density due to Radiation', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_einrad', atm_energy%einrad,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%einrad)
    END IF
    !
    !
    IF (ltend_einradvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_content_due_to_radiation', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Content due to Radiation', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_einradvi', atm_energy%einradvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%einradvi)
    END IF
    !
    !
    IF (ltend_einradhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_profile_due_to_radiation', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Profile due to Radiation', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_einradhi', atm_energy%einradhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%einradhi)
    END IF
    !
    !
    IF (ltend_einradti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_due_to_radiation', 'J s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy due to Radiation', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_einradti', atm_energy%einradti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%einradti)
    END IF


    ! tmx: due to turbulent mixing
    ! ----------------------------

    ! - of internal energy 
    !
    IF (ltend_eintmx) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_density_due_to_turbulent_mixing', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Density due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_eintmx', atm_energy%eintmx,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eintmx)
    END IF
    !
    !
    IF (ltend_eintmxvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_content_due_to_turbulent_mixing', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Content due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_eintmxvi', atm_energy%eintmxvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eintmxvi)
    END IF
    !
    !
    IF (ltend_eintmxhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_profile_due_to_turbulent_mixing', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Profile due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_eintmxhi', atm_energy%eintmxhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eintmxhi)
    END IF
    !
    !
    IF (ltend_eintmxti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_due_to_turbulent_mixing', 'J s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_eintmxti', atm_energy%eintmxti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%eintmxti)
    END IF


    ! - of horizontal kinetic energy
    !
    IF (ltend_ekhtmx) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_density_due_to_turbulent_mixing', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy Density due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekhtmx', atm_energy%ekhtmx,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhtmx)
    END IF
    !
    !
    IF (ltend_ekhtmxvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_content_due_to_turbulent_mixing', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy Content due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekhtmxvi', atm_energy%ekhtmxvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhtmxvi)
    END IF
    !
    !
    IF (ltend_ekhtmxhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_profile_due_to_turbulent_mixing', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy Profile due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekhtmxhi', atm_energy%ekhtmxhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhtmxhi)
    END IF
    !
    !
    IF (ltend_ekhtmxti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_due_to_turbulent_mixing', 'J s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekhtmxti', atm_energy%ekhtmxti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhtmxti)
    END IF


    ! - of vertical kinetic energy
    !
    IF (ltend_ekvtmx) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_density_due_to_turbulent_mixing', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy Density due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekvtmx', atm_energy%ekvtmx,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvtmx)
    END IF
    !
    !
    IF (ltend_ekvtmxvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_content_due_to_turbulent_mixing', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy Content due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekvtmxvi', atm_energy%ekvtmxvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvtmxvi)
    END IF
    !
    !
    IF (ltend_ekvtmxhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_profile_due_to_turbulent_mixing', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy Profile due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekvtmxhi', atm_energy%ekvtmxhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvtmxhi)
    END IF
    !
    !
    IF (ltend_ekvtmxti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_due_to_turbulent_mixing', 'J s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy due to Turbulent Mixing', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekvtmxti', atm_energy%ekvtmxti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvtmxti)
    END IF


    ! phy: due to physics
    ! --------------------

    ! - of internal energy 
    !
    IF (ltend_einphy) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_density_due_to_physics', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Density due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_einphy', atm_energy%einphy,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%einphy)
    END IF
    !
    !
    IF (ltend_einphyvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_content_due_to_physics', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Content due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_einphyvi', atm_energy%einphyvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%einphyvi)
    END IF
    !
    !
    IF (ltend_einphyhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_profile_due_to_physics', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy Profile due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_einphyhi', atm_energy%einphyhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%einphyhi)
    END IF
    !
    !
    IF (ltend_einphyti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_moist_internal_energy_due_to_physics', 'J s-1', &
            &             'Tendency of Atmosphere Moist Internal Energy due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_einphyti', atm_energy%einphyti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%einphyti)
    END IF


    ! - of horizontal kinetic energy
    !
    IF (ltend_ekhphy) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_density_due_to_physics', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy Density due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekhphy', atm_energy%ekhphy,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhphy)
    END IF
    !
    !
    IF (ltend_ekhphyvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_content_due_to_physics', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy Content due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekhphyvi', atm_energy%ekhphyvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhphyvi)
    END IF
    !
    !
    IF (ltend_ekhphyhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_profile_due_to_physics', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy Profile due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekhphyhi', atm_energy%ekhphyhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhphyhi)
    END IF
    !
    !
    IF (ltend_ekhphyti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_horizontal_kinetic_energy_due_to_physics', 'J s-1', &
            &             'Tendency of Atmosphere Horizontal Kinetic Energy due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekhphyti', atm_energy%ekhphyti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekhphyti)
    END IF


    ! - of vertical kinetic energy
    !
    IF (ltend_ekvphy) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_density_due_to_physics', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy Density due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekvphy', atm_energy%ekvphy,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvphy)
    END IF
    !
    !
    IF (ltend_ekvphyvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_content_due_to_physics', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy Content due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_ekvphyvi', atm_energy%ekvphyvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvphyvi)
    END IF
    !
    !
    IF (ltend_ekvphyhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_profile_due_to_physics', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy Profile due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekvphyhi', atm_energy%ekvphyhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvphyhi)
    END IF
    !
    !
    IF (ltend_ekvphyti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_vertical_kinetic_energy_due_to_physics', 'J s-1', &
            &             'Tendency of Atmosphere Vertical Kinetic Energy due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_ekvphyti', atm_energy%ekvphyti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%ekvphyti)
    END IF


    ! of geopotential energy
    !
    IF (ltend_egpphy) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_geopotential_energy_density_due_to_physics', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Geopotential Energy Density due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_egpphy', atm_energy%egpphy,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%egpphy)
    END IF
    !
    !
    IF (ltend_egpphyvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_geopotential_energy_content_due_to_physics', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Geopotential Energy Content due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_egpphyvi', atm_energy%egpphyvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%egpphyvi)
    END IF
    !
    !
    IF (ltend_egpphyhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_geopotential_energy_profile_due_to_physics', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Geopotential Energy Profile due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_egpphyhi', atm_energy%egpphyhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%egpphyhi)
    END IF
    !
    !
    IF (ltend_egpphyti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_geopotential_energy_due_to_physics', 'J s-1', &
            &             'Tendency of Atmosphere Geopotential Energy due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_egpphyti', atm_energy%egpphyti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%egpphyti)
    END IF


    ! of total energy
    !
    IF (ltend_etophy) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_total_energy_density_due_to_physics', 'J m-3 s-1', &
            &             'Tendency of Atmosphere Total Energy Density due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_etophy', atm_energy%etophy,     &
            &       grid_unstructured_cell, za_reference,                  &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape3d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%etophy)
    END IF
    !
    !
    IF (ltend_etophyvi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_total_energy_content_due_to_physics', 'J m-2 s-1', &
            &             'Tendency of Atmosphere Total Energy Content due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_cell)
       !
       CALL add_var(atm_energy_list, 'tend_etophyvi', atm_energy%etophyvi, &
            &       grid_unstructured_cell, za_surface,                    &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape2d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%etophyvi)
    END IF
    !
    !
    IF (ltend_etophyhi) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_total_energy_profile_due_to_physics', 'J m-1 s-1', &
            &             'Tendency of Atmosphere Total Energy Profile due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_etophyhi', atm_energy%etophyhi, &
            &       grid_lonlat, za_reference,                             &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape1d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%etophyhi)
    END IF
    !
    !
    IF (ltend_etophyti) THEN
       cf_desc = t_cf_var('tendency_of_atmosphere_total_energy_due_to_physics', 'J s-1', &
            &             'Tendency of Atmosphere Total Energy due to Physics', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, datatype_grb, grid_unstructured, grid_lonlat)
       !
       CALL add_var(atm_energy_list, 'tend_etophyti', atm_energy%etophyti, &
            &       grid_lonlat, za_surface,                               &
            &       cf_desc, grib2_desc,                                   &
            &       ldims=shape0d,                                         &
            &       lrestart = .FALSE.,                                    &
            &       isteptype=tstep_instant,                               &
            &       lopenacc=.TRUE.)
       __acc_attach(atm_energy%etophyti)
    END IF


  END SUBROUTINE construct_atm_energy_list

END MODULE mo_atm_energy_memory
