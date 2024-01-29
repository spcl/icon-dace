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

MODULE mo_cf_convention
  !!
  !!----------------------------------------------------------------------------
  !!
  !! File name space:
  !!
  !! title
  !!
  !!    A succint description of what is in the dataset. 
  !!
  !! institution
  !!
  !!    Specifies where the original data was produced. 
  !!
  !! source
  !!
  !!    The method of production of the original data. If it 
  !!    was model-generated, source should name the model and 
  !!    its version, as specifically as could be useful.
  !!
  !! history
  !!
  !!    Provides an audit trail for modifications to the original 
  !!    data. Well-behaved generic netCDF filters will automatically 
  !!    append their name and the parameters with which they were 
  !!    invoked to the global history attribute of an input netCDF file. 
  !!
  !! references
  !!
  !!    Published or web-based references that describe the data or methods 
  !!    used to produce it. 
  !!
  !! comment
  !!
  !!    Miscellaneous information about the data or methods used to produce it. 
  !!
  !!
  !!----------------------------------------------------------------------------
  !!
  !! Variable name space:
  !!
  !! standard name
  !!
  !!    The name used to identify the physical quantity. A standard name 
  !!    contains no whitespace and is case sensitive. 
  !!
  !! canonical units
  !!
  !!    Representative units of the physical quantity. Unless it is 
  !!    dimensionless, a variable with a standard_name attribute must have 
  !!    units which are physically equivalent (not necessarily identical) to 
  !!    the canonical units, possibly modified by an operation specified by 
  !!    either the standard name modifier.
  !!
  !!  description
  !!
  !!    The description is meant to clarify the qualifiers of the fundamental 
  !!    quantities such as which surface a quantity is defined on or what the 
  !!    flux sign conventions are. We don"t attempt to provide precise 
  !!    definitions of fundumental physical quantities (e.g., temperature) 
  !!    which may be found in the literature. 
  !!
  !!----------------------------------------------------------------------------

  IMPLICIT NONE 

  PRIVATE

  TYPE t_cf_global
    CHARACTER(len=128) :: title         = ''
    CHARACTER(len=128) :: institution   = ''
    CHARACTER(len=128) :: source        = ''
    CHARACTER(len=128) :: history       = ''
    CHARACTER(len=128) :: references    = ''
    CHARACTER(len=128) :: comment       = ''
  END TYPE t_cf_global

  TYPE t_cf_var
    CHARACTER(len=128) :: standard_name = ''
    CHARACTER(len=128) :: units         = ''
    CHARACTER(len=256) :: long_name     = ''
    INTEGER            :: datatype      = -1
    CHARACTER(len=256) :: short_name    = ''
  END TYPE t_cf_var

  TYPE t_cf_gridspec
    CHARACTER(len=36)   :: gridspec_coordinates_id   ! uuid
    CHARACTER(len=36)   :: gridspec_data_id          ! uuid
    CHARACTER(len= 9)   :: gridspec_file_type        ! always: grid_file       
    CHARACTER(len=1024) :: gridspec_tile_name        ! URL of grid file
  END type t_cf_gridspec

  PUBLIC :: t_cf_global
  PUBLIC :: t_cf_var
  PUBLIC :: t_cf_gridspec

  PUBLIC :: set_cf_global
  PUBLIC :: set_cf_var
  PUBLIC :: set_cf_gridspec

  TYPE(t_cf_global), SAVE, PUBLIC, PROTECTED :: cf_global_info
  
CONTAINS

  SUBROUTINE set_cf_global(title, institution, source, history, references, comment)
    CHARACTER(len=*), INTENT(in), OPTIONAL :: title
    CHARACTER(len=*), INTENT(in), OPTIONAL :: institution
    CHARACTER(len=*), INTENT(in), OPTIONAL :: source
    CHARACTER(len=*), INTENT(in), OPTIONAL :: history
    CHARACTER(len=*), INTENT(in), OPTIONAL :: references
    CHARACTER(len=*), INTENT(in), OPTIONAL :: comment
   
    IF (PRESENT(title))       cf_global_info%title       = TRIM(title)
    IF (PRESENT(institution)) cf_global_info%institution = TRIM(institution)
    IF (PRESENT(source))      cf_global_info%source      = TRIM(source) 
    IF (PRESENT(history))     cf_global_info%history     = TRIM(history)  
    IF (PRESENT(references))  cf_global_info%references  = TRIM(references)
    IF (PRESENT(comment))     cf_global_info%comment     = TRIM(comment)

  END SUBROUTINE set_cf_global

  FUNCTION set_cf_var(long_name, units, standard_name, datatype) &
       RESULT(cf_var_info)
    TYPE(t_cf_var) :: cf_var_info
    CHARACTER(len=*), INTENT(in), OPTIONAL :: long_name
    CHARACTER(len=*), INTENT(in), OPTIONAL :: units
    CHARACTER(len=*), INTENT(in), OPTIONAL :: standard_name
    INTEGER,          INTENT(in), OPTIONAL :: datatype

    IF (PRESENT(long_name))     cf_var_info%long_name     = long_name
    IF (PRESENT(units))         cf_var_info%units         = units
    IF (PRESENT(standard_name)) cf_var_info%standard_name = standard_name
    IF (PRESENT(datatype))      cf_var_info%datatype      = datatype

  END FUNCTION set_cf_var

  FUNCTION set_cf_gridspec(coordinate_uuid, data_uuid, grid_file) &
       RESULT(cf_gridspec)
    TYPE(t_cf_gridspec) :: cf_gridspec
    CHARACTER(len=*), INTENT(in) :: coordinate_uuid
    CHARACTER(len=*), INTENT(in) :: data_uuid
    CHARACTER(len=*), INTENT(in) :: grid_file

    cf_gridspec%gridspec_coordinates_id = TRIM(coordinate_uuid)
    cf_gridspec%gridspec_data_id        = TRIM(data_uuid)
    cf_gridspec%gridspec_file_type      = 'grid_file'
    cf_gridspec%gridspec_tile_name      = TRIM(grid_file)

  END FUNCTION set_cf_gridspec

END MODULE mo_cf_convention
