! Auto-generated file by "/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/sdfgs/utils/dace-genfi.py"
module mo_radiation_bindings

  use iso_c_binding

  use mo_exception, only: &
    warning
  use mo_impl_constants, only: &
    MAX_CHAR_LENGTH
  use radiation_aerosol, only: &
    aerosol_type
  use radiation_cloud, only: &
    cloud_type
  use radiation_config, only: &
    config_type
  use radiation_flux, only: &
    flux_type
  use radiation_gas, only: &
    gas_type
  use radiation_single_level, only: &
    single_level_type
  use radiation_thermodynamics, only: &
    thermodynamics_type

  implicit none

  private
  public :: run_radiation
  public :: run_radiation_verification
  public :: verify_radiation
  public :: dace_init_radiation
  public :: dace_exit_radiation
  public :: dace_program_radiation


  type, bind(c) :: dace_aerosol_type

  end type dace_aerosol_type

  type, bind(c) :: dace_cloud_type

  end type dace_cloud_type

  type, bind(c) :: dace_config_type
    integer(kind=c_int) :: f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8
    integer(kind=c_int) :: f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9
    integer(kind=c_int) :: f2dace_SA_i_emiss_from_band_lw_d_0_s_5
    integer(kind=c_int) :: f2dace_SA_sw_albedo_weights_d_0_s_6
    integer(kind=c_int) :: f2dace_SA_sw_albedo_weights_d_1_s_7
    integer(kind=c_int) :: f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_8
    integer(kind=c_int) :: f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9
    integer(kind=c_int) :: f2dace_SOA_i_emiss_from_band_lw_d_0_s_5
    integer(kind=c_int) :: f2dace_SOA_sw_albedo_weights_d_0_s_6
    integer(kind=c_int) :: f2dace_SOA_sw_albedo_weights_d_1_s_7
    type(c_ptr) :: i_band_from_reordered_g_lw
    type(c_ptr) :: i_band_from_reordered_g_sw
    type(c_ptr) :: i_emiss_from_band_lw
    type(c_ptr) :: sw_albedo_weights

  end type dace_config_type

  type, bind(c) :: dace_flux_type

  end type dace_flux_type

  type, bind(c) :: dace_gas_type

  end type dace_gas_type

  type, bind(c) :: dace_single_level_type
    integer(kind=c_int) :: f2dace_SA_lw_emissivity_d_0_s_20
    integer(kind=c_int) :: f2dace_SA_lw_emissivity_d_1_s_21
    integer(kind=c_int) :: f2dace_SA_sw_albedo_d_0_s_16
    integer(kind=c_int) :: f2dace_SA_sw_albedo_d_1_s_17
    integer(kind=c_int) :: f2dace_SA_sw_albedo_direct_d_0_s_18
    integer(kind=c_int) :: f2dace_SA_sw_albedo_direct_d_1_s_19
    integer(kind=c_int) :: f2dace_SOA_lw_emissivity_d_0_s_20
    integer(kind=c_int) :: f2dace_SOA_lw_emissivity_d_1_s_21
    integer(kind=c_int) :: f2dace_SOA_sw_albedo_d_0_s_16
    integer(kind=c_int) :: f2dace_SOA_sw_albedo_d_1_s_17
    integer(kind=c_int) :: f2dace_SOA_sw_albedo_direct_d_0_s_18
    integer(kind=c_int) :: f2dace_SOA_sw_albedo_direct_d_1_s_19
    type(c_ptr) :: lw_emissivity
    type(c_ptr) :: sw_albedo
    type(c_ptr) :: sw_albedo_direct

  end type dace_single_level_type

  type, bind(c) :: dace_thermodynamics_type

  end type dace_thermodynamics_type


  logical :: is_initialized = .false.
  type(c_ptr) :: dace_state = C_NULL_PTR

  type(c_ptr) :: cached_shallow_copy_aerosol = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_cloud = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_config = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_flux = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_gas = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_single_level = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_thermodynamics = C_NULL_PTR

  type(c_ptr) :: verification_deep_copy_aerosol = C_NULL_PTR
  type(c_ptr) :: verification_deep_copy_cloud = C_NULL_PTR
  type(c_ptr) :: verification_deep_copy_config = C_NULL_PTR
  type(c_ptr) :: verification_deep_copy_flux = C_NULL_PTR
  type(c_ptr) :: verification_deep_copy_gas = C_NULL_PTR
  type(c_ptr) :: verification_deep_copy_single_level = C_NULL_PTR
  type(c_ptr) :: verification_deep_copy_thermodynamics = C_NULL_PTR

interface

  type(c_ptr) function malloc(size) &
    bind(c, name="malloc")
    use iso_c_binding
    integer(kind=c_size_t), value :: size
  end function malloc

  subroutine free(ptr) &
    bind(c, name="free")
    use iso_c_binding
    type(c_ptr), value :: ptr
  end subroutine free


  type(c_ptr) function dace_init_radiation( &
    aerosol, &
    cloud, &
    config, &
    flux, &
    gas, &
    single_level, &
    thermodynamics, &
    iendcol, &
    istartcol, &
    ncol, &
    nlev, &
    nulout, &
    sym_iendcol, &
    sym_istartcol &
  ) &
    bind(c, name="__dace_init_radiation")
    use iso_c_binding

    type(c_ptr), value :: aerosol
    type(c_ptr), value :: cloud
    type(c_ptr), value :: config
    type(c_ptr), value :: flux
    type(c_ptr), value :: gas
    type(c_ptr), value :: single_level
    type(c_ptr), value :: thermodynamics
    integer(kind=c_int), value :: iendcol
    integer(kind=c_int), value :: istartcol
    integer(kind=c_int), value :: ncol
    integer(kind=c_int), value :: nlev
    integer(kind=c_int), value :: nulout
    integer(kind=c_int), value :: sym_iendcol
    integer(kind=c_int), value :: sym_istartcol
  end function dace_init_radiation

  integer(c_int) function dace_exit_radiation(state) &
    bind(c, name="__dace_exit_radiation")
    use iso_c_binding

    type(c_ptr), value :: state
  end function dace_exit_radiation

  subroutine dace_program_radiation( &
    state, &
    aerosol, &
    cloud, &
    config, &
    flux, &
    gas, &
    single_level, &
    thermodynamics, &
    iendcol, &
    istartcol, &
    ncol, &
    nlev, &
    nulout, &
    sym_iendcol, &
    sym_istartcol &
  ) &
    bind(c, name="__program_radiation")
    use iso_c_binding

    type(c_ptr), value :: state
    type(c_ptr), value :: aerosol
    type(c_ptr), value :: cloud
    type(c_ptr), value :: config
    type(c_ptr), value :: flux
    type(c_ptr), value :: gas
    type(c_ptr), value :: single_level
    type(c_ptr), value :: thermodynamics
    integer(kind=c_int), value :: iendcol
    integer(kind=c_int), value :: istartcol
    integer(kind=c_int), value :: ncol
    integer(kind=c_int), value :: nlev
    integer(kind=c_int), value :: nulout
    integer(kind=c_int), value :: sym_iendcol
    integer(kind=c_int), value :: sym_istartcol
  end subroutine dace_program_radiation

end interface

interface logical_fix_1d
  module procedure logical_to_int_1d
  module procedure int_to_int_1d
end interface logical_fix_1d


  real(8), parameter :: float32_default_rel_threshold = 1.0e-8
  real(8), parameter :: float32_default_abs_threshold = 0.0

  real(8), parameter :: float64_default_rel_threshold = 1.0e-12
  real(8), parameter :: float64_default_abs_threshold = 0.0

  real(8), parameter :: int32_default_rel_threshold = 1.0e-8
  real(8), parameter :: int32_default_abs_threshold = 0.0

  real(8), parameter :: int64_default_rel_threshold = 1.0e-12
  real(8), parameter :: int64_default_abs_threshold = 0.0


contains

  function copy_in_aerosol_type(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(aerosol_type), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_aerosol_type), pointer :: dace_rich_obj

    dace_obj_ptr = malloc(c_sizeof(dace_rich_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)


  end function copy_in_aerosol_type

  function copy_in_cloud_type(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(cloud_type), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_cloud_type), pointer :: dace_rich_obj

    dace_obj_ptr = malloc(c_sizeof(dace_rich_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)


  end function copy_in_cloud_type

  function copy_in_config_type(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(config_type), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_config_type), pointer :: dace_rich_obj

    dace_obj_ptr = malloc(c_sizeof(dace_rich_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

        dace_rich_obj%f2dace_SA_i_emiss_from_band_lw_d_0_s_5 = 0
            dace_rich_obj%f2dace_SOA_i_emiss_from_band_lw_d_0_s_5 = 0
            dace_rich_obj%f2dace_SA_sw_albedo_weights_d_0_s_6 = 0
            dace_rich_obj%f2dace_SOA_sw_albedo_weights_d_0_s_6 = 0
            dace_rich_obj%f2dace_SA_sw_albedo_weights_d_1_s_7 = 0
            dace_rich_obj%f2dace_SOA_sw_albedo_weights_d_1_s_7 = 0
            dace_rich_obj%f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8 = 0
            dace_rich_obj%f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_8 = 0
            dace_rich_obj%f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9 = 0
            dace_rich_obj%f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9 = 0
            dace_rich_obj%i_emiss_from_band_lw = c_null_ptr
            dace_rich_obj%sw_albedo_weights = c_null_ptr
            dace_rich_obj%i_band_from_reordered_g_lw = c_null_ptr
            dace_rich_obj%i_band_from_reordered_g_sw = c_null_ptr
    
  end function copy_in_config_type

  function copy_in_flux_type(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(flux_type), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_flux_type), pointer :: dace_rich_obj

    dace_obj_ptr = malloc(c_sizeof(dace_rich_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)


  end function copy_in_flux_type

  function copy_in_gas_type(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(gas_type), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_gas_type), pointer :: dace_rich_obj

    dace_obj_ptr = malloc(c_sizeof(dace_rich_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)


  end function copy_in_gas_type

  function copy_in_single_level_type(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(single_level_type), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_single_level_type), pointer :: dace_rich_obj

    dace_obj_ptr = malloc(c_sizeof(dace_rich_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

        dace_rich_obj%f2dace_SA_sw_albedo_d_0_s_16 = 0
            dace_rich_obj%f2dace_SOA_sw_albedo_d_0_s_16 = 0
            dace_rich_obj%f2dace_SA_sw_albedo_d_1_s_17 = 0
            dace_rich_obj%f2dace_SOA_sw_albedo_d_1_s_17 = 0
            dace_rich_obj%f2dace_SA_sw_albedo_direct_d_0_s_18 = 0
            dace_rich_obj%f2dace_SOA_sw_albedo_direct_d_0_s_18 = 0
            dace_rich_obj%f2dace_SA_sw_albedo_direct_d_1_s_19 = 0
            dace_rich_obj%f2dace_SOA_sw_albedo_direct_d_1_s_19 = 0
            dace_rich_obj%f2dace_SA_lw_emissivity_d_0_s_20 = 0
            dace_rich_obj%f2dace_SOA_lw_emissivity_d_0_s_20 = 0
            dace_rich_obj%f2dace_SA_lw_emissivity_d_1_s_21 = 0
            dace_rich_obj%f2dace_SOA_lw_emissivity_d_1_s_21 = 0
            dace_rich_obj%sw_albedo = c_null_ptr
            dace_rich_obj%sw_albedo_direct = c_null_ptr
            dace_rich_obj%lw_emissivity = c_null_ptr
    
  end function copy_in_single_level_type

  function copy_in_thermodynamics_type(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(thermodynamics_type), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_thermodynamics_type), pointer :: dace_rich_obj

    dace_obj_ptr = malloc(c_sizeof(dace_rich_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)


  end function copy_in_thermodynamics_type

  function logical_to_int_1d(inp) result(out)
    logical(4), dimension(:), target :: inp
    integer(kind=c_int), dimension(:), pointer :: out

    call c_f_pointer(c_loc(inp), out, shape=shape(inp))
  end function logical_to_int_1d

  function int_to_int_1d(inp) result(out)
    integer(kind=c_int), dimension(:), target :: inp
    integer(kind=c_int), dimension(:), pointer :: out

    out => inp
  end function int_to_int_1d

  function copy_in_int32_1d_array(fortran_array, steal_arrays, minimal_structs) result(dace_array_ptr)
    integer(kind=c_int), dimension(:), target :: fortran_array
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_array_ptr
    integer(kind=c_int), dimension(:), pointer :: dace_rich_array

    integer :: i0

    if (.not. c_associated(c_loc(fortran_array))) then
      dace_array_ptr = c_null_ptr
      return
    end if


    if (steal_arrays .eqv. .true.) then
      dace_array_ptr = c_loc(fortran_array)
      return
    end if

    dace_array_ptr = malloc(c_sizeof(dace_array_ptr) * size(fortran_array))
    call c_f_pointer(dace_array_ptr, dace_rich_array, shape=shape(fortran_array))

    do i0 = 1, size(fortran_array, dim=1)
      dace_rich_array(i0) = fortran_array(i0)
    end do

  end function copy_in_int32_1d_array

  function copy_in_float64_2d_array(fortran_array, steal_arrays, minimal_structs) result(dace_array_ptr)
    real(kind=c_double), dimension(:,:), target :: fortran_array
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_array_ptr
    real(kind=c_double), dimension(:,:), pointer :: dace_rich_array

    integer :: i0, i1

    if (.not. c_associated(c_loc(fortran_array))) then
      dace_array_ptr = c_null_ptr
      return
    end if


    if (steal_arrays .eqv. .true.) then
      dace_array_ptr = c_loc(fortran_array)
      return
    end if

    dace_array_ptr = malloc(c_sizeof(dace_array_ptr) * size(fortran_array))
    call c_f_pointer(dace_array_ptr, dace_rich_array, shape=shape(fortran_array))

    do i0 = 1, size(fortran_array, dim=1)
      do i1 = 1, size(fortran_array, dim=2)
        dace_rich_array(i0, i1) = fortran_array(i0, i1)
      end do
    end do

  end function copy_in_float64_2d_array


  subroutine copy_back_aerosol_type(fortran_obj, dace_obj_ptr)
    type(aerosol_type), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_aerosol_type), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        call warning( &
          "copy_back_aerosol_type", &
          "Invalid allocation of aerosol_type by DaCe!" &
        )
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



  end subroutine copy_back_aerosol_type
  subroutine copy_back_cloud_type(fortran_obj, dace_obj_ptr)
    type(cloud_type), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_cloud_type), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        call warning( &
          "copy_back_cloud_type", &
          "Invalid allocation of cloud_type by DaCe!" &
        )
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



  end subroutine copy_back_cloud_type
  subroutine copy_back_config_type(fortran_obj, dace_obj_ptr)
    type(config_type), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_config_type), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        call warning( &
          "copy_back_config_type", &
          "Invalid allocation of config_type by DaCe!" &
        )
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



  end subroutine copy_back_config_type
  subroutine copy_back_flux_type(fortran_obj, dace_obj_ptr)
    type(flux_type), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_flux_type), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        call warning( &
          "copy_back_flux_type", &
          "Invalid allocation of flux_type by DaCe!" &
        )
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



  end subroutine copy_back_flux_type
  subroutine copy_back_gas_type(fortran_obj, dace_obj_ptr)
    type(gas_type), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_gas_type), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        call warning( &
          "copy_back_gas_type", &
          "Invalid allocation of gas_type by DaCe!" &
        )
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



  end subroutine copy_back_gas_type
  subroutine copy_back_single_level_type(fortran_obj, dace_obj_ptr)
    type(single_level_type), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_single_level_type), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        call warning( &
          "copy_back_single_level_type", &
          "Invalid allocation of single_level_type by DaCe!" &
        )
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



  end subroutine copy_back_single_level_type
  subroutine copy_back_thermodynamics_type(fortran_obj, dace_obj_ptr)
    type(thermodynamics_type), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_thermodynamics_type), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        call warning( &
          "copy_back_thermodynamics_type", &
          "Invalid allocation of thermodynamics_type by DaCe!" &
        )
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



  end subroutine copy_back_thermodynamics_type

  subroutine check_initializations()
    CHARACTER(len=MAX_CHAR_LENGTH) :: message_text = ''

  end subroutine check_initializations



  subroutine compare_float32_scalar( &
    actual, &
    ref, &
    result, &
    rel_threshold, &
    abs_threshold, &
    scalar_expr &
  )
    real(kind=c_float), intent(in) :: actual, ref
    logical, intent(out) :: result
    real(8), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in), optional :: scalar_expr

    real(8) :: rel_error, abs_error, threshold_ratio
    real(8) :: actual_rel_threshold, actual_abs_threshold
    CHARACTER(len=MAX_CHAR_LENGTH) :: message_text = ''

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = float32_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = float32_default_abs_threshold
    end if

    result = abs(ref - actual) <= max(actual_rel_threshold * abs(ref), actual_abs_threshold)

    if (present(scalar_expr) .and. .not. result) then

      threshold_ratio = real(abs(ref - actual), kind=8) / max(actual_rel_threshold * abs(ref), actual_abs_threshold)
      rel_error = abs(real(ref - actual, kind=8)/ref)
      abs_error = abs(ref - actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,e28.20,a,e28.20)') &
        "Verification failed for scalar '", &
          trim(scalar_expr), &
        "'"//new_line('a')//"    - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
          ", threshold_ratio = ", &
          threshold_ratio, &
        new_line('a')//"    - ref = ", &
          ref, &
          ", actual = ", &
          actual
      call warning("compare_float32_scalar", message_text)

    end if

  end subroutine compare_float32_scalar


  subroutine compare_float64_scalar( &
    actual, &
    ref, &
    result, &
    rel_threshold, &
    abs_threshold, &
    scalar_expr &
  )
    real(kind=c_double), intent(in) :: actual, ref
    logical, intent(out) :: result
    real(8), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in), optional :: scalar_expr

    real(8) :: rel_error, abs_error, threshold_ratio
    real(8) :: actual_rel_threshold, actual_abs_threshold
    CHARACTER(len=MAX_CHAR_LENGTH) :: message_text = ''

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = float64_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = float64_default_abs_threshold
    end if

    result = abs(ref - actual) <= max(actual_rel_threshold * abs(ref), actual_abs_threshold)

    if (present(scalar_expr) .and. .not. result) then

      threshold_ratio = real(abs(ref - actual), kind=8) / max(actual_rel_threshold * abs(ref), actual_abs_threshold)
      rel_error = abs(real(ref - actual, kind=8)/ref)
      abs_error = abs(ref - actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,e28.20,a,e28.20)') &
        "Verification failed for scalar '", &
          trim(scalar_expr), &
        "'"//new_line('a')//"    - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
          ", threshold_ratio = ", &
          threshold_ratio, &
        new_line('a')//"    - ref = ", &
          ref, &
          ", actual = ", &
          actual
      call warning("compare_float64_scalar", message_text)

    end if

  end subroutine compare_float64_scalar


  subroutine compare_int32_scalar( &
    actual, &
    ref, &
    result, &
    rel_threshold, &
    abs_threshold, &
    scalar_expr &
  )
    integer(kind=c_int), intent(in) :: actual, ref
    logical, intent(out) :: result
    real(8), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in), optional :: scalar_expr

    real(8) :: rel_error, abs_error, threshold_ratio
    real(8) :: actual_rel_threshold, actual_abs_threshold
    CHARACTER(len=MAX_CHAR_LENGTH) :: message_text = ''

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = int32_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = int32_default_abs_threshold
    end if

    result = abs(ref - actual) <= max(actual_rel_threshold * abs(ref), actual_abs_threshold)

    if (present(scalar_expr) .and. .not. result) then

      threshold_ratio = real(abs(ref - actual), kind=8) / max(actual_rel_threshold * abs(ref), actual_abs_threshold)
      rel_error = abs(real(ref - actual, kind=8)/ref)
      abs_error = abs(ref - actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,e28.20,a,e28.20)') &
        "Verification failed for scalar '", &
          trim(scalar_expr), &
        "'"//new_line('a')//"    - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
          ", threshold_ratio = ", &
          threshold_ratio, &
        new_line('a')//"    - ref = ", &
          ref, &
          ", actual = ", &
          actual
      call warning("compare_int32_scalar", message_text)

    end if

  end subroutine compare_int32_scalar


  subroutine compare_int64_scalar( &
    actual, &
    ref, &
    result, &
    rel_threshold, &
    abs_threshold, &
    scalar_expr &
  )
    integer(kind=c_long), intent(in) :: actual, ref
    logical, intent(out) :: result
    real(8), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in), optional :: scalar_expr

    real(8) :: rel_error, abs_error, threshold_ratio
    real(8) :: actual_rel_threshold, actual_abs_threshold
    CHARACTER(len=MAX_CHAR_LENGTH) :: message_text = ''

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = int64_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = int64_default_abs_threshold
    end if

    result = abs(ref - actual) <= max(actual_rel_threshold * abs(ref), actual_abs_threshold)

    if (present(scalar_expr) .and. .not. result) then

      threshold_ratio = real(abs(ref - actual), kind=8) / max(actual_rel_threshold * abs(ref), actual_abs_threshold)
      rel_error = abs(real(ref - actual, kind=8)/ref)
      abs_error = abs(ref - actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,e28.20,a,e28.20)') &
        "Verification failed for scalar '", &
          trim(scalar_expr), &
        "'"//new_line('a')//"    - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
          ", threshold_ratio = ", &
          threshold_ratio, &
        new_line('a')//"    - ref = ", &
          ref, &
          ", actual = ", &
          actual
      call warning("compare_int64_scalar", message_text)

    end if

  end subroutine compare_int64_scalar

  subroutine compare_aerosol_type_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(aerosol_type), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=MAX_CHAR_LENGTH) :: member_expr = ''
    type(dace_aerosol_type), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.


    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_aerosol_type_struct

  subroutine compare_cloud_type_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(cloud_type), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=MAX_CHAR_LENGTH) :: member_expr = ''
    type(dace_cloud_type), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.


    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_cloud_type_struct

  subroutine compare_config_type_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(config_type), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=MAX_CHAR_LENGTH) :: member_expr = ''
    type(dace_config_type), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%i_emiss_from_band_lw"

    call compare_int32_1d_array( &
        actual=actual_rich%i_emiss_from_band_lw, &
        ref=logical_fix_1d(ref%i_emiss_from_band_lw), &
        result=local_result, &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%sw_albedo_weights"

    call compare_float64_2d_array( &
        actual=actual_rich%sw_albedo_weights, &
        ref=ref%sw_albedo_weights, &
        result=local_result, &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%i_band_from_reordered_g_lw"

    call compare_int32_1d_array( &
        actual=actual_rich%i_band_from_reordered_g_lw, &
        ref=logical_fix_1d(ref%i_band_from_reordered_g_lw), &
        result=local_result, &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%i_band_from_reordered_g_sw"

    call compare_int32_1d_array( &
        actual=actual_rich%i_band_from_reordered_g_sw, &
        ref=logical_fix_1d(ref%i_band_from_reordered_g_sw), &
        result=local_result, &
        array_expr=member_expr &
    )

    result = result .and. local_result


    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_config_type_struct

  subroutine compare_flux_type_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(flux_type), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=MAX_CHAR_LENGTH) :: member_expr = ''
    type(dace_flux_type), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.


    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_flux_type_struct

  subroutine compare_gas_type_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(gas_type), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=MAX_CHAR_LENGTH) :: member_expr = ''
    type(dace_gas_type), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.


    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_gas_type_struct

  subroutine compare_single_level_type_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(single_level_type), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=MAX_CHAR_LENGTH) :: member_expr = ''
    type(dace_single_level_type), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%sw_albedo"

    call compare_float64_2d_array( &
        actual=actual_rich%sw_albedo, &
        ref=ref%sw_albedo, &
        result=local_result, &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%sw_albedo_direct"

    call compare_float64_2d_array( &
        actual=actual_rich%sw_albedo_direct, &
        ref=ref%sw_albedo_direct, &
        result=local_result, &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%lw_emissivity"

    call compare_float64_2d_array( &
        actual=actual_rich%lw_emissivity, &
        ref=ref%lw_emissivity, &
        result=local_result, &
        array_expr=member_expr &
    )

    result = result .and. local_result


    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_single_level_type_struct

  subroutine compare_thermodynamics_type_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(thermodynamics_type), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=MAX_CHAR_LENGTH) :: member_expr = ''
    type(dace_thermodynamics_type), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.


    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_thermodynamics_type_struct

  subroutine compare_int32_1d_array( &
    actual, &
    ref, &
    result, &
    array_expr, &
    rel_threshold, &
    abs_threshold &
  )
    type(c_ptr) :: actual
    integer(kind=c_int), dimension(:), target, intent(in) :: ref
    logical, intent(out) :: result
    real(kind=c_double), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in) :: array_expr

    real(kind=c_double) :: actual_rel_threshold, actual_abs_threshold
    logical :: local_result
    integer :: i0
    integer(kind=c_int), dimension(:), pointer :: actual_rich
    CHARACTER(len=MAX_CHAR_LENGTH) :: message_text = ''

    integer(kind=c_int) :: error_ref, error_actual
    integer, dimension(0:0) :: max_threshold_ratio_loc
    real(8) :: rel_error, abs_error, threshold_ratio

    integer :: max_threshold_ratio_i0

    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//new_line('a')//"    - ref was NULL, but actual was not!"
        call warning("compare_int32_1d_array", message_text)
      end if

      return
    end if

    result = .true.

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = int32_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = int32_default_abs_threshold
    end if

    call c_f_pointer(actual, actual_rich, shape=shape(ref))


    do i0 = 1, size(ref, dim=1)

    call compare_int32_scalar( &
      actual=actual_rich(i0), &
      ref=transfer(ref(i0), mold=int(1, kind=4)), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

    result = result .and. local_result
    end do


    if (.not. result) then
      max_threshold_ratio_loc = maxloc(abs(ref - actual_rich) / max(actual_rel_threshold * abs(ref), actual_abs_threshold))
      max_threshold_ratio_i0 = max_threshold_ratio_loc(0)

      error_ref = ref(max_threshold_ratio_i0)
      error_actual = actual_rich(max_threshold_ratio_i0)

      threshold_ratio = real(abs(error_ref - error_actual), kind=8) / max(actual_rel_threshold * abs(error_ref), actual_abs_threshold)
      rel_error = abs(real(error_ref - error_actual, kind=8)/error_ref)
      abs_error = abs(error_ref - error_actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,1(i0:,", "),a,e28.20,a,e28.20)') &
        "Verification failed for array '", &
          trim(array_expr), &
        "':"//new_line('a')//"    - max_threshold_ratio = ", &
          threshold_ratio, &
          ", rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
        new_line('a')//"    - at (", &
          [max_threshold_ratio_i0], &
          "), ref = ", &
          error_ref, &
          ", actual = ", &
          error_actual
      call warning("compare_int32_1d_array", message_text)

    end if

    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_int32_1d_array

  subroutine compare_float64_2d_array( &
    actual, &
    ref, &
    result, &
    array_expr, &
    rel_threshold, &
    abs_threshold &
  )
    type(c_ptr) :: actual
    real(kind=c_double), dimension(:,:), target, intent(in) :: ref
    logical, intent(out) :: result
    real(kind=c_double), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in) :: array_expr

    real(kind=c_double) :: actual_rel_threshold, actual_abs_threshold
    logical :: local_result
    integer :: i0, i1
    real(kind=c_double), dimension(:,:), pointer :: actual_rich
    CHARACTER(len=MAX_CHAR_LENGTH) :: message_text = ''

    real(kind=c_double) :: error_ref, error_actual
    integer, dimension(0:1) :: max_threshold_ratio_loc
    real(8) :: rel_error, abs_error, threshold_ratio

    integer :: max_threshold_ratio_i0, max_threshold_ratio_i1

    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//new_line('a')//"    - ref was NULL, but actual was not!"
        call warning("compare_float64_2d_array", message_text)
      end if

      return
    end if

    result = .true.

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = float64_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = float64_default_abs_threshold
    end if

    call c_f_pointer(actual, actual_rich, shape=shape(ref))


    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)

    call compare_float64_scalar( &
      actual=actual_rich(i0, i1), &
      ref=ref(i0, i1), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

    result = result .and. local_result
      end do
    end do


    if (.not. result) then
      max_threshold_ratio_loc = maxloc(abs(ref - actual_rich) / max(actual_rel_threshold * abs(ref), actual_abs_threshold))
      max_threshold_ratio_i0 = max_threshold_ratio_loc(0); max_threshold_ratio_i1 = max_threshold_ratio_loc(1)

      error_ref = ref(max_threshold_ratio_i0, max_threshold_ratio_i1)
      error_actual = actual_rich(max_threshold_ratio_i0, max_threshold_ratio_i1)

      threshold_ratio = real(abs(error_ref - error_actual), kind=8) / max(actual_rel_threshold * abs(error_ref), actual_abs_threshold)
      rel_error = abs(real(error_ref - error_actual, kind=8)/error_ref)
      abs_error = abs(error_ref - error_actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,2(i0:,", "),a,e28.20,a,e28.20)') &
        "Verification failed for array '", &
          trim(array_expr), &
        "':"//new_line('a')//"    - max_threshold_ratio = ", &
          threshold_ratio, &
          ", rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
        new_line('a')//"    - at (", &
          [max_threshold_ratio_i0, max_threshold_ratio_i1], &
          "), ref = ", &
          error_ref, &
          ", actual = ", &
          error_actual
      call warning("compare_float64_2d_array", message_text)

    end if

    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_float64_2d_array


    
  subroutine run_radiation( &
    aerosol, &
    cloud, &
    config, &
    flux, &
    gas, &
    single_level, &
    thermodynamics, &
    iendcol, &
    istartcol, &
    ncol, &
    nlev, &
    nulout, &
    sym_iendcol, &
    sym_istartcol &
  )
    type(aerosol_type), target :: aerosol
    type(cloud_type), target :: cloud
    type(config_type), target :: config
    type(flux_type), target :: flux
    type(gas_type), target :: gas
    type(single_level_type), target :: single_level
    type(thermodynamics_type), target :: thermodynamics
    integer(kind=c_int) :: iendcol
    integer(kind=c_int) :: istartcol
    integer(kind=c_int) :: ncol
    integer(kind=c_int) :: nlev
    integer(kind=c_int) :: nulout
    integer(kind=c_int) :: sym_iendcol
    integer(kind=c_int) :: sym_istartcol




    call check_initializations()

    verification_deep_copy_aerosol = copy_in_aerosol_type( &
    fortran_obj=aerosol, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    verification_deep_copy_cloud = copy_in_cloud_type( &
    fortran_obj=cloud, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    verification_deep_copy_config = copy_in_config_type( &
    fortran_obj=config, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    verification_deep_copy_flux = copy_in_flux_type( &
    fortran_obj=flux, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    verification_deep_copy_gas = copy_in_gas_type( &
    fortran_obj=gas, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    verification_deep_copy_single_level = copy_in_single_level_type( &
    fortran_obj=single_level, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    verification_deep_copy_thermodynamics = copy_in_thermodynamics_type( &
    fortran_obj=thermodynamics, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )


    if (is_initialized .eqv. .false.) then
      is_initialized = .true.
      dace_state = dace_init_radiation( &
      aerosol = verification_deep_copy_aerosol, &
      cloud = verification_deep_copy_cloud, &
      config = verification_deep_copy_config, &
      flux = verification_deep_copy_flux, &
      gas = verification_deep_copy_gas, &
      single_level = verification_deep_copy_single_level, &
      thermodynamics = verification_deep_copy_thermodynamics, &
      iendcol = iendcol, &
      istartcol = istartcol, &
      ncol = ncol, &
      nlev = nlev, &
      nulout = nulout, &
      sym_iendcol = sym_iendcol, &
      sym_istartcol = sym_istartcol &
    )
    end if

    call dace_program_radiation( &
      state = dace_state, &
      aerosol = verification_deep_copy_aerosol, &
      cloud = verification_deep_copy_cloud, &
      config = verification_deep_copy_config, &
      flux = verification_deep_copy_flux, &
      gas = verification_deep_copy_gas, &
      single_level = verification_deep_copy_single_level, &
      thermodynamics = verification_deep_copy_thermodynamics, &
      iendcol = iendcol, &
      istartcol = istartcol, &
      ncol = ncol, &
      nlev = nlev, &
      nulout = nulout, &
      sym_iendcol = sym_iendcol, &
      sym_istartcol = sym_istartcol &
    )

    call copy_back_aerosol_type(aerosol, verification_deep_copy_aerosol)
    call copy_back_cloud_type(cloud, verification_deep_copy_cloud)
    call copy_back_config_type(config, verification_deep_copy_config)
    call copy_back_flux_type(flux, verification_deep_copy_flux)
    call copy_back_gas_type(gas, verification_deep_copy_gas)
    call copy_back_single_level_type(single_level, verification_deep_copy_single_level)
    call copy_back_thermodynamics_type(thermodynamics, verification_deep_copy_thermodynamics)


  end subroutine run_radiation

  subroutine run_radiation_verification( &
    aerosol, &
    cloud, &
    config, &
    flux, &
    gas, &
    single_level, &
    thermodynamics, &
    iendcol, &
    istartcol, &
    ncol, &
    nlev, &
    nulout, &
    sym_iendcol, &
    sym_istartcol &
  )
    type(aerosol_type), target :: aerosol
    type(cloud_type), target :: cloud
    type(config_type), target :: config
    type(flux_type), target :: flux
    type(gas_type), target :: gas
    type(single_level_type), target :: single_level
    type(thermodynamics_type), target :: thermodynamics
    integer(kind=c_int) :: iendcol
    integer(kind=c_int) :: istartcol
    integer(kind=c_int) :: ncol
    integer(kind=c_int) :: nlev
    integer(kind=c_int) :: nulout
    integer(kind=c_int) :: sym_iendcol
    integer(kind=c_int) :: sym_istartcol




    call check_initializations()

    verification_deep_copy_aerosol = copy_in_aerosol_type( &
    fortran_obj=aerosol, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    verification_deep_copy_cloud = copy_in_cloud_type( &
    fortran_obj=cloud, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    verification_deep_copy_config = copy_in_config_type( &
    fortran_obj=config, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    verification_deep_copy_flux = copy_in_flux_type( &
    fortran_obj=flux, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    verification_deep_copy_gas = copy_in_gas_type( &
    fortran_obj=gas, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    verification_deep_copy_single_level = copy_in_single_level_type( &
    fortran_obj=single_level, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    verification_deep_copy_thermodynamics = copy_in_thermodynamics_type( &
    fortran_obj=thermodynamics, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )


    if (is_initialized .eqv. .false.) then
      is_initialized = .true.
      dace_state = dace_init_radiation( &
      aerosol = verification_deep_copy_aerosol, &
      cloud = verification_deep_copy_cloud, &
      config = verification_deep_copy_config, &
      flux = verification_deep_copy_flux, &
      gas = verification_deep_copy_gas, &
      single_level = verification_deep_copy_single_level, &
      thermodynamics = verification_deep_copy_thermodynamics, &
      iendcol = iendcol, &
      istartcol = istartcol, &
      ncol = ncol, &
      nlev = nlev, &
      nulout = nulout, &
      sym_iendcol = sym_iendcol, &
      sym_istartcol = sym_istartcol &
    )
    end if

    call dace_program_radiation( &
      state = dace_state, &
      aerosol = verification_deep_copy_aerosol, &
      cloud = verification_deep_copy_cloud, &
      config = verification_deep_copy_config, &
      flux = verification_deep_copy_flux, &
      gas = verification_deep_copy_gas, &
      single_level = verification_deep_copy_single_level, &
      thermodynamics = verification_deep_copy_thermodynamics, &
      iendcol = iendcol, &
      istartcol = istartcol, &
      ncol = ncol, &
      nlev = nlev, &
      nulout = nulout, &
      sym_iendcol = sym_iendcol, &
      sym_istartcol = sym_istartcol &
    )
  end subroutine run_radiation_verification

  subroutine verify_radiation( &
    aerosol, &
    cloud, &
    config, &
    flux, &
    gas, &
    single_level, &
    thermodynamics, &
    iendcol, &
    istartcol, &
    ncol, &
    nlev, &
    nulout, &
    sym_iendcol, &
    sym_istartcol &
  )
    type(aerosol_type), target :: aerosol
    type(cloud_type), target :: cloud
    type(config_type), target :: config
    type(flux_type), target :: flux
    type(gas_type), target :: gas
    type(single_level_type), target :: single_level
    type(thermodynamics_type), target :: thermodynamics
    integer(kind=c_int) :: iendcol
    integer(kind=c_int) :: istartcol
    integer(kind=c_int) :: ncol
    integer(kind=c_int) :: nlev
    integer(kind=c_int) :: nulout
    integer(kind=c_int) :: sym_iendcol
    integer(kind=c_int) :: sym_istartcol


    logical :: local_result, result = .true.


    if (is_initialized .eqv. .false.) then
      call warning("verify_radiation", "dace state is not initialized")
    end if

    call check_initializations()


    call compare_aerosol_type_struct( &
        actual=verification_deep_copy_aerosol, &
        ref=aerosol, &
        result=local_result, &
        struct_expr="aerosol" &
    )

    result = result .and. local_result

    call compare_cloud_type_struct( &
        actual=verification_deep_copy_cloud, &
        ref=cloud, &
        result=local_result, &
        struct_expr="cloud" &
    )

    result = result .and. local_result

    call compare_config_type_struct( &
        actual=verification_deep_copy_config, &
        ref=config, &
        result=local_result, &
        struct_expr="config" &
    )

    result = result .and. local_result

    call compare_flux_type_struct( &
        actual=verification_deep_copy_flux, &
        ref=flux, &
        result=local_result, &
        struct_expr="flux" &
    )

    result = result .and. local_result

    call compare_gas_type_struct( &
        actual=verification_deep_copy_gas, &
        ref=gas, &
        result=local_result, &
        struct_expr="gas" &
    )

    result = result .and. local_result

    call compare_single_level_type_struct( &
        actual=verification_deep_copy_single_level, &
        ref=single_level, &
        result=local_result, &
        struct_expr="single_level" &
    )

    result = result .and. local_result

    call compare_thermodynamics_type_struct( &
        actual=verification_deep_copy_thermodynamics, &
        ref=thermodynamics, &
        result=local_result, &
        struct_expr="thermodynamics" &
    )

    result = result .and. local_result

    if (.not. result) then
      call warning("verify_radiation", "Failed verification")
    else
      call warning("verify_radiation", "Verification successful :)")
    end if

  end subroutine verify_radiation

end module mo_radiation_bindings
