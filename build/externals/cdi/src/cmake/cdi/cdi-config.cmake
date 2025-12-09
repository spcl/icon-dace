# Config file for lib CDI project
# run cmake with CDI_DIR pointing to the directory containing this file.

# set(CDI_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")
# use following line instead of previous if this file can be configured 
set(CDI_ROOT_DIR "/usr/local")

find_path(CDI_INCLUDE_DIRECTORY
  cdi.h
  PATHS "${CDI_ROOT_DIR}/include"
  # use following line instead of previous if this file can be configured 
  # PATHS "${prefix}/include"
  DOC "The CDI include directory"
  )
mark_as_advanced(CDI_INCLUDE_DIRECTORY)

find_library(CDI_LIBRARY
  NAMES cdi
  PATHS "${CDI_ROOT_DIR}/lib"
  # use following line instead of previous if this file can be configured 
  # PATHS "${exec_prefix}/lib"
  DOC "The CDI library"
  )
mark_as_advanced(CDI_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CDI
  REQUIRED_VARS CDI_INCLUDE_DIRECTORY CDI_LIBRARY
  )

if (CDI_FOUND)
  set(CDI_INCLUDE_DIRS "${CDI_INCLUDE_DIRECTORY}")
  set(CDI_LIBRARIES "${CDI_LIBRARY}")

  if (NOT TARGET CDI::CDI)
    add_library(CDI::CDI UNKNOWN IMPORTED)
    set_target_properties(CDI::CDI PROPERTIES
      IMPORTED_LOCATION "${CDI_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${CDI_INCLUDE_DIRS}")
  endif()
endif()
