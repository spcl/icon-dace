add_library(utils STATIC utils.F90 netcdf_utils.F90)
target_link_libraries(utils ComIn NetCDF::NetCDF_Fortran)
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
set_property(TARGET utils PROPERTY POSITION_INDEPENDENT_CODE ON)
