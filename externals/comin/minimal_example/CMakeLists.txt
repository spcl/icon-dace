find_package(OpenMP REQUIRED)
find_package(NetCDF REQUIRED COMPONENTS Fortran)
find_package(YAC)

add_executable(minimal_example
  minimal_example.F90 descr_data.F90 vars.F90 mo_utilities.F90)
add_executable(ComIn::minimal_example ALIAS minimal_example)

if(YAC_FOUND)
  # We need to link YAC before NetCDF because CMake does not resolve
  # the linking order correctly as both YAC and NetCDF_Fortran depend
  # on NetCDF_C
  target_link_libraries(minimal_example YAC::YAC)
  target_compile_definitions(minimal_example PRIVATE -DENABLE_YAC)
endif()

target_link_libraries(minimal_example
    -Wl,--dynamic-list=${CMAKE_CURRENT_SOURCE_DIR}/dynamic_list.sym)

target_link_libraries(minimal_example
  ComIn MPI::MPI_Fortran OpenMP::OpenMP_Fortran NetCDF::NetCDF_Fortran)

install(
  TARGETS minimal_example
  EXPORT "${PROJECT_NAME}Targets"
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

add_subdirectory(grids)
