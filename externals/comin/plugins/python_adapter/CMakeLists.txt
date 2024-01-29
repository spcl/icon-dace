cmake_minimum_required (VERSION 3.17)
project (python_adapter LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo" FORCE)
  message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified")
endif(NOT CMAKE_BUILD_TYPE)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Python3 REQUIRED COMPONENTS Development Interpreter)
find_package(ComIn REQUIRED)

# TODO find alls cpp files
add_library(python_adapter MODULE python_adapter.cpp callbacks.cpp variables.cpp descrdata.cpp)
target_link_libraries(python_adapter Python3::Python ComIn::ComIn)

add_custom_command( OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/comin.py.h
  COMMAND xxd -i comin.py comin.py.h
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/comin.py
  COMMENT "Compiling comin.py to comin.py.h")

add_custom_target(
    cominpyh
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/comin.py.h
  )

add_dependencies(python_adapter cominpyh)

option(BUILD_TESTING "" OFF) # turn testing off by default as it introduces a lot of dependencies
include(CTest)

enable_testing()
IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
