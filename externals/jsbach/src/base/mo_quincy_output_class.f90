!> QUINCY model output parameters
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### defines QUINCY model parameters, units and helpers for model output
!>
MODULE mo_quincy_output_class
#ifndef __NO_QUINCY__

  USE mo_kind,              ONLY: wp
  USE mo_io_units,          ONLY: filename_max
  USE mo_exception,         ONLY: finish
  USE mo_util,              ONLY: int2string
  USE mo_jsb_var_class,     ONLY: t_jsb_var_p, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist_iface, ONLY: VARNAME_LEN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_quincy_output, t_quincy_outputfile, t_quincy_out_var
  PUBLIC :: NO_LAYERS, SOIL_LAYERED, CANOPY_LAYERED
  PUBLIC :: NO_FACT, USE_NDAT_ACC, USE_VEG_FLUX_FACT
  PUBLIC :: unit_flux, unit_pool, unit_soillayer_flux, unit_soillayer_pool, unit_cn, unit_nc, unit_np, unit_pn, unit_cp
  PUBLIC :: unitless, unit_hydro_pool, unit_hydro_flux, unit_water_potential, unit_veg_pool_flux, unit_assimi_rates
  PUBLIC :: unit_fraction, unit_energy, unit_temperature, unit_speed, unit_size, unit_bulk_dens, unit_rain
  PUBLIC :: unit_conductivity

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_quincy_output_class'

  ! units
  CHARACTER(len=*), PARAMETER :: unit_flux            = 'micro mol m-2 s-1'
  CHARACTER(len=*), PARAMETER :: unit_pool            = 'mol m-2'
  CHARACTER(len=*), PARAMETER :: unit_soillayer_flux  = 'micro mol m-3 s-1'
  CHARACTER(len=*), PARAMETER :: unit_soillayer_pool  = 'mol m-3'
  !@TODO: unit? mol/mol?
  CHARACTER(len=*), PARAMETER :: unit_cn              = 'mol C mol-1 N'
  CHARACTER(len=*), PARAMETER :: unit_cp              = 'mol C mol-1 P'
  CHARACTER(len=*), PARAMETER :: unit_nc              = 'mol N mol-1 C'
  CHARACTER(len=*), PARAMETER :: unit_np              = 'mol N mol-1 P'
  CHARACTER(len=*), PARAMETER :: unit_pn              = 'mol P mol-1 N'

  CHARACTER(len=*), PARAMETER :: unit_assimi_rates    = 'micro mol CO2 m-2 s-1'

  CHARACTER(len=*), PARAMETER :: unit_bulk_dens       = 'kg m-3'

  CHARACTER(len=*), PARAMETER :: unit_rain            = 'mm s-1'
  CHARACTER(len=*), PARAMETER :: unit_hydro_flux      = 'kg m-2 s-1'
  CHARACTER(len=*), PARAMETER :: unit_hydro_pool      = 'm'
  CHARACTER(len=*), PARAMETER :: unit_water_potential = 'MPa'

  CHARACTER(len=*), PARAMETER :: unit_size            = "m"
  CHARACTER(len=*), PARAMETER :: unit_speed           = "m s-1"
  CHARACTER(len=*), PARAMETER :: unit_temperature     = "K"
  CHARACTER(len=*), PARAMETER :: unit_fraction        = "ppm"
  CHARACTER(len=*), PARAMETER :: unit_energy          = "W m-2"
  CHARACTER(len=*), PARAMETER :: unitless             = '--'
  CHARACTER(len=*), PARAMETER :: unit_conductivity    = "m s-1"

  CHARACTER(len=*), PARAMETER :: unit_veg_pool_flux   = 'mol m-2 timestep-1'

  ! vertical grid layers
  ENUM, BIND(C)
  ENUMERATOR ::      &
    & NO_LAYERS,     &
    & SOIL_LAYERED,  &
    & CANOPY_LAYERED
  END ENUM

  ! factor
  ENUM, BIND(C)
  ENUMERATOR ::   &
    & NO_FACT, &
    & USE_VEG_FLUX_FACT, &
    & USE_NDAT_ACC
  END ENUM

  !> Type used to collect output files
  !>
  TYPE :: t_quincy_output
    INTEGER                             :: nr_of_output_files  = 0
      !< counter for number of output files
    REAL(wp)                            :: ndat_acc            = 0._wp
      !< "timestep counter for averaging" (from mo_qs_output_txt)
    !@TODO: rename or change comment
    TYPE(t_quincy_outputfile), POINTER  :: first_output_file   => NULL()
      !< pointer to first quincy output file
    TYPE(t_quincy_outputfile), POINTER  :: current_output_file => NULL()
      !< pointer to quincy output file in use
    CHARACTER(len=filename_max), POINTER :: sel_output_variables(:) => NULL()
      !< selected output variables
  END TYPE t_quincy_output

  !> Type for an output file
  !>
  TYPE :: t_quincy_outputfile
    CHARACTER(len=200)                  :: file_name_tag
      !< @TODO constant?
    INTEGER                             :: file_id            = -1
      !< assigned and used by netcdf4
    INTEGER                             :: nr_of_dims         = -1
      !<
    LOGICAL                             :: has_soil_layers    = .FALSE.
      !< set to true if file has variables with soil layers
    LOGICAL                             :: has_canopy_layers  = .FALSE.
      !< set to true if file has variables with canopy layers
    INTEGER                             :: current_record     = 1
      !< counter for record dimension (here: 'time')
    CHARACTER(len=filename_max), POINTER :: sel_output_variables(:) => NULL()
      !< Selected input variables
    TYPE(t_quincy_outputfile), POINTER  :: next_output_file   => NULL()
      !< pointer to next quincy output file
    TYPE(t_quincy_outputfile), POINTER  :: prev_output_file   => NULL()
      !< pointer to previous quincy output file
    TYPE(t_quincy_out_var),    POINTER  :: first_var          => NULL()
      !< pointer to the first quincy variable of this file
    TYPE(t_quincy_out_var),    POINTER  :: current_var        => NULL()
      !< pointer to the current quincy variable
  END TYPE t_quincy_outputfile

  !> Type for an output variable
  !>
  TYPE :: t_quincy_out_var
    CHARACTER(len=VARNAME_LEN)      :: var_name
      !< @TODO: maybe too short
    CHARACTER(len=VARNAME_LEN)      :: unit
      !< @TODO: can probably be much shorter! Same for unit_l...
    CHARACTER(len=VARNAME_LEN)      :: unit_sl
      !< only relevant for soil layered variables: unit of 2D sl var
    INTEGER                         :: nr_of_dims       = -1
      !< 1 for 2D and 2 for 3D variables
    INTEGER                         :: var_type         = NO_LAYERS
      !< indicates the var type: if not, soil or canopy layered
    INTEGER                         :: var_id           = -1
      !< assigned and used by netcdf4
    INTEGER                         :: var_add_2d_id    = -1
      !< required for additional file id in case of layered main variable
    TYPE(t_jsb_var_p), ALLOCATABLE  :: var_container(:) !< collects the pointer(s) to the variable(s)
    LOGICAL                         :: l_is_acc_var     = .FALSE.
    LOGICAL                         :: l_avg_sl         = .FALSE.
      !< only relevant for soil layered variables: if false they are averaged, else (default) they are summed
    LOGICAL                         :: l_write_2d       = .TRUE.
      !< only relevant for soil layered variables: write var per sl
    LOGICAL                         :: l_write_1d       = .TRUE.
      !< only relevant for soil layered variables: write sum/avg of sl var
    INTEGER                         :: this_element !TODO: consider renaming
      !< set to the id of an element, if var is an element itself or a collection of compartments
    INTEGER                         :: out_fact_kind    = NO_FACT
      !< if a conversion/averaging factor is to be used upon output writing
    REAL(wp)                        :: accumulated      = 0._wp
    REAL(wp), ALLOCATABLE           :: accumulated_2D(:)

    INTEGER                         :: dim2_len = 0
      !< Size of 2 dim (nsoil_layer or ncanopy_layer)

    ! buffered data to be written only once per year
    REAL(wp), ALLOCATABLE           :: buffer_2D(:)
    REAL(wp), ALLOCATABLE           :: buffer_3D(:,:)

    TYPE(t_quincy_out_var), POINTER :: next_var         => NULL()
      !< pointer to the next output variable in this file
    TYPE(t_quincy_out_var), POINTER :: prev_var         => NULL()
      !< pointer to the previous output variable in this file
  END TYPE t_quincy_out_var

#endif
END MODULE mo_quincy_output_class
