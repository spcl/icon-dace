!> Contains routines for var_lists/streams
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
MODULE mo_jsb_varlist
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp, dp
  USE mo_exception,         ONLY: finish
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d, t_jsb_var_real3d, REAL2D, REAL3D
  USE mo_jsb_varlist_iface, ONLY: t_jsb_varlist => t_var_list,      &
                                  jsb_get_var_list => get_var_list, &
                                  jsb_new_var_list => new_var_list, &
#ifdef __ICON__
                                  add_var_list_element_r1d,         &
#endif
                                  add_var_list_element_r2d, add_var_list_element_r3d, &
                                  is_variable_in_output, &
                                  VARNAME_LEN
  USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2, GRID_UNSTRUCTURED, GRID_UNSTRUCTURED_CELL, GRID_CELL, missval, tables
  USE mo_jsb_io_iface,      ONLY: t_cf_var, t_grib2_var, grib2_var

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)  
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_varlist
#ifdef __ICON__
  PUBLIC :: add_var_list_element_r1d
#endif
  PUBLIC :: add_var_list_element_r2d, add_var_list_element_r3d
  PUBLIC :: jsb_add_var, Get_varlist, Get_var_group_name
  PUBLIC :: VARNAME_LEN
  PUBLIC :: BASIC, MEDIUM, FULL

  INTERFACE jsb_add_var
    MODULE PROCEDURE jsb_add_var_list_element_r2d
  END INTERFACE jsb_add_var

  ENUM, BIND(C)
    ENUMERATOR :: BASIC=1, MEDIUM, FULL
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_varlist'

CONTAINS

  FUNCTION Get_varlist(name, patch_id) RESULT(return_ptr)

    CHARACTER(LEN=*), INTENT(in) :: name
    INTEGER,          INTENT(in) :: patch_id
    TYPE(t_jsb_varlist), POINTER :: return_ptr

    CALL jsb_get_var_list(return_ptr, TRIM(name))

    IF (.NOT. ASSOCIATED(return_ptr)) THEN
      ALLOCATE(return_ptr)
      CALL jsb_new_var_list(return_ptr, TRIM(name), patch_id, &
        & loutput=.TRUE., &
        & lrestart=.TRUE., table=tables(1))
    END IF

  END FUNCTION Get_varlist

  SUBROUTINE jsb_add_var_list_element_r2d(varlist_name, name, ptr, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix,             &
    loutput, lcontainer,                                        &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r,                               &
    groups, verbose                                             &
    )

    CHARACTER(len=*),     INTENT(in)           :: varlist_name
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:)            ! reference to field
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               ! horizontal grid
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               ! vertical grid
    TYPE(t_cf),           INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               ! GRIB1 related metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               ! GRIB2 related metadata
    CHARACTER(len=*),     INTENT(in)           :: prefix              ! Prefix for variable name
    CHARACTER(len=*),     INTENT(in)           :: suffix              ! Suffix for variable name
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    CHARACTER(LEN=VARNAME_LEN), INTENT(in), OPTIONAL :: groups(:)
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information

!!$    REAL(wp), POINTER :: ptr(:,:,:,:,:)
    TYPE(t_jsb_varlist), POINTER :: this_list           ! list
    CHARACTER(len=50) :: name_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(dp)          :: missval_r_loc
    INTEGER           :: ldims(2), gdims(2)
    INTEGER           :: n_groups
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups_loc(:)     ! groups to which this variable belongs to

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsb_add_var_list_element_r2d'

    name_loc = TRIM(name)
    IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
    IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)

    lmiss_loc = .TRUE.
    IF (PRESENT(lmiss)) lmiss_loc = lmiss
    missval_r_loc = missval
    IF (PRESENT(missval_r)) missval_r_loc = missval_r

    cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
    grib2_loc = grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_UNSTRUCTURED, GRID_CELL)

    ldims(1) = hgrid%nproma
    ldims(2) = hgrid%nblks
    gdims(1:2) = hgrid%dims_g(1:2)

    this_list => Get_varlist(TRIM(varlist_name), hgrid%host_patch_id)

    IF (PRESENT(groups)) THEN
      n_groups = SIZE(groups) + 2
    ELSE
      n_groups = 2
    END IF
    ALLOCATE(groups_loc(n_groups))
    groups_loc(1) = 'ALL'
    groups_loc(2) = 'jsb_all'
    IF (PRESENT(groups)) THEN
      groups_loc(3:n_groups) = groups(:)
    END IF

    CALL add_var_list_element_r2d( &
      this_list, TRIM(name_loc), ptr, GRID_UNSTRUCTURED_CELL, vgrid%ZaxisID, &
      cf_loc, grib2_loc, code=grib1%parameter, table=grib1%table, &
      ldims=ldims, gdims=gdims, levelindx=vgrid%echamZaxisIdx, &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
      lmiss=lmiss_loc, missval_r=missval_r_loc, in_groups=groups_loc, verbose=verbose)

  END SUBROUTINE jsb_add_var_list_element_r2d

  FUNCTION Get_var_group_name(id) RESULT(group_name)

    INTEGER, INTENT(in) :: id
    CHARACTER(len=:), ALLOCATABLE :: group_name

    SELECT CASE (id)
    CASE(BASIC)
      group_name = 'basic'
    CASE(MEDIUM)
      group_name = 'medium'
    CASE(FULL)
      group_name = 'full'
    END SELECT

  END FUNCTION

#endif
END MODULE mo_jsb_varlist
