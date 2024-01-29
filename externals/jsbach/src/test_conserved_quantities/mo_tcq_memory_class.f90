!> tcq (test conserved quantities) memory class
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Contains variables in the memory class of the test conserved quantities process 
!>
MODULE mo_tcq_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_cqt_class,          ONLY: DEAD_CARBON_CQ_TYPE, LIVE_CARBON_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  dsl4jsb_Use_processes TCQ_
  dsl4jsb_Use_config(TCQ_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_tcq_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 10

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_tcq_memory
    TYPE(t_jsb_var_real2d) :: old_sum, conservation_test_field, a_dead_c_tcq, another_dead_c_tcq, &
      & a_veg_c_tcq, an_implicit_scaling_tcq ! some to-be-conserved carbon quantities
    TYPE(t_jsb_var_real2d) :: a_dead_c_ta_tcq, another_dead_c_ta_tcq, &
      & a_veg_c_ta_tcq, an_implicit_scaling_ta_tcq ! corresponding ta variables
  CONTAINS
    PROCEDURE :: Init => Init_tcq_memory
  END TYPE t_tcq_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tcq_memory_class'

CONTAINS

  SUBROUTINE Init_tcq_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, t_jsb_varlist
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_tcq_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(TCQ_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_tcq_memory'
    ! -------------------------------------------------------------------------------------------------- !

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    dsl4jsb_Get_config(TCQ_)

    CALL mem%Add_var('a_veg_c_tcq', mem%a_veg_c_tcq,                                     &
      & hgrid, surface,                                                              &
      & t_cf('a_veg_c_tcq', 'mol(C) m-2(canopy)', 'A C-Pool.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix, output_level=BASIC,                                      &
      & l_conserve_quan=.TRUE., cons_quan_type_id = LIVE_CARBON_CQ_TYPE,     &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

    CALL mem%Add_var('a_dead_c_tcq', mem%a_dead_c_tcq,                                     &
      & hgrid, surface,                                                              &
      & t_cf('a_dead_c_tcq', 'mol(C) m-2(canopy)', 'A C-Pool.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix, output_level=BASIC,                                           &
      & l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,     &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL mem%Add_var('another_dead_c_tcq', mem%another_dead_c_tcq,                                     &
      & hgrid, surface,                                                              &
      & t_cf('another_dead_c_tcq', 'mol(C) m-2(canopy)', 'A C-Pool.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix, output_level=BASIC,                                     &
      & l_conserve_quan=.TRUE., cons_quan_type_id = DEAD_CARBON_CQ_TYPE,     &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL mem%Add_var('an_implicit_scaling_tcq', mem%an_implicit_scaling_tcq,                      &
      & hgrid, surface,                                                              &
      & t_cf('an_implicit_scaling_tcq', 'mol(C) m-2(canopy)', 'A C-Pool.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix, output_level=BASIC,                                                 &
      & l_conserve_quan=.TRUE., cons_quan_type_id = PRODUCT_CARBON_CQ_TYPE,     &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL mem%Add_var('a_veg_c_ta_tcq', mem%a_veg_c_ta_tcq,                                     &
      & hgrid, surface,                                                              &
      & t_cf('a_veg_c_ta_tcq', 'mol(C) m-2', 'A C-Pool.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC,      &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL mem%Add_var('a_dead_c_ta_tcq', mem%a_dead_c_ta_tcq,                                     &
      & hgrid, surface,                                                              &
      & t_cf('a_dead_c_ta_tcq', 'mol(C) m-2', 'A C-Pool.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC,  &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL mem%Add_var('another_dead_c_ta_tcq', mem%another_dead_c_ta_tcq,                                     &
      & hgrid, surface,                                                              &
      & t_cf('another_dead_c_ta_tcq', 'mol(C)', 'A C-Pool.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC,      &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL mem%Add_var('an_implicit_scaling_ta_tcq', mem%an_implicit_scaling_ta_tcq,                      &
      & hgrid, surface,                                                              &
      & t_cf('an_implicit_scaling_ta_tcq', 'mol(C) m-2', 'A C-Pool.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC,      &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL mem%Add_var('conservation_test_field', mem%conservation_test_field,                      &
      & hgrid, surface,                                                              &
      & t_cf('conservation_test_field', 'mol(C) m-2', 'Neither scaled nor relocated -- tcq test.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC,      &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp) 

      CALL mem%Add_var('old_sum', mem%old_sum,                      &
      & hgrid, surface,                                                              &
      & t_cf('old_sum', 'mol(C) m-2', 'Neither scaled nor relocated C-Pool -- old tcq sum'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC,      &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp) 
      

  END SUBROUTINE Init_tcq_memory

#endif
END MODULE mo_tcq_memory_class
