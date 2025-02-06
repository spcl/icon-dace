#----------------------------------------------------------------
# Generated CMake target import file for configuration "NOEXTRAFLAGS".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "math-support::math-support" for configuration "NOEXTRAFLAGS"
set_property(TARGET math-support::math-support APPEND PROPERTY IMPORTED_CONFIGURATIONS NOEXTRAFLAGS)
set_target_properties(math-support::math-support PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOEXTRAFLAGS "Fortran"
  IMPORTED_LOCATION_NOEXTRAFLAGS "${_IMPORT_PREFIX}/lib/libmath-support.a"
  )

list(APPEND _cmake_import_check_targets math-support::math-support )
list(APPEND _cmake_import_check_files_for_math-support::math-support "${_IMPORT_PREFIX}/lib/libmath-support.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
