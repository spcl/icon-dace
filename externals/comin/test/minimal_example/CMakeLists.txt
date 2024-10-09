find_package(OpenMP REQUIRED)
find_package(NetCDF REQUIRED COMPONENTS Fortran)
find_package(YAC)

add_executable(minimal_example
  minimal_example.F90 descr_data.F90 vars.F90 mo_utilities.F90)
add_executable(${PROJECT_NAME}::minimal_example ALIAS minimal_example)

set_target_properties(minimal_example PROPERTIES
  EXPORT_NAME "${PROJECT_NAME}::minimal_example"
)

if(YAC_FOUND)
  # We need to link YAC before NetCDF because CMake does not resolve
  # the linking order correctly as both YAC and NetCDF_Fortran depend
  # on NetCDF_C
  target_link_libraries(minimal_example YAC::YAC)
  target_compile_definitions(minimal_example PRIVATE -DENABLE_YAC)
endif()

target_link_libraries(minimal_example
  ComIn::ComIn MPI::MPI_Fortran OpenMP::OpenMP_Fortran NetCDF::NetCDF_Fortran)

install(TARGETS minimal_example EXPORT "${PROJECT_NAME}Targets")

add_subdirectory(grids)
