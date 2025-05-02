# cmake options:
#
#       -DCMAKE_INSTALL_PREFIX=/path/to/install

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( cdi VERSION 2.4.0 LANGUAGES C )
set(CMAKE_C_STANDARD 11)

configure_file (
  "${PROJECT_SOURCE_DIR}/cmake/cdi_config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

include_directories("${PROJECT_BINARY_DIR}")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(NetCDF)

set(HAVE_NETCDF ${netCDF_FOUND} )

if (netCDF_FOUND)
  set(HAVE_NETCDF 1)
endif ()

add_subdirectory( src )
add_subdirectory( app )
#target_link_libraries(cdi PRIVATE NetCDF)
