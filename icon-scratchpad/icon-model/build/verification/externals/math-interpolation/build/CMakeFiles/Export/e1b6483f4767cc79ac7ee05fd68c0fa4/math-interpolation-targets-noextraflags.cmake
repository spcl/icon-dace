#----------------------------------------------------------------
# Generated CMake target import file for configuration "NOEXTRAFLAGS".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "math-interpolation::math-interpolation" for configuration "NOEXTRAFLAGS"
set_property(TARGET math-interpolation::math-interpolation APPEND PROPERTY IMPORTED_CONFIGURATIONS NOEXTRAFLAGS)
set_target_properties(math-interpolation::math-interpolation PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOEXTRAFLAGS "Fortran"
  IMPORTED_LOCATION_NOEXTRAFLAGS "${_IMPORT_PREFIX}/lib/libmath-interpolation.a"
  )

list(APPEND _cmake_import_check_targets math-interpolation::math-interpolation )
list(APPEND _cmake_import_check_files_for_math-interpolation::math-interpolation "${_IMPORT_PREFIX}/lib/libmath-interpolation.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
