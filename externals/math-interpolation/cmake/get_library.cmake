# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

# cmake-format: off
# get_libary(<package_target>
#            <tag>
#            [QUIET])
# cmake-format: on
# -----------------------------------------------------------------------------
# <package> and <target> are found from <package_target> which contains the
# package name and the target name separated by '::'. If '::' is not found,
# <package> and <target> are set to <package_target>. If <target> is not given,
# the function will find <package> and use it as <target>. Once the <package>
# and <target> are found, this function will check if <target> is already in the
# TARGETS list. If not, it will find the package using find_package() and if not
# found, it will fetch the package from the git repository using the commit
# <tag>.
#
function(get_library package_target tag)

  cmake_parse_arguments(PARSE_ARGV 1 ARG "QUIET" "" "")

  # Find the position of the double colon '::' in the variable
  string(FIND "${package_target}" "::" double_colon_index)

  # If '::' is found in the variable, extract 'package' and 'target'
  if(NOT ${double_colon_index} EQUAL -1)
    # Extract 'package'
    string(SUBSTRING "${package_target}" 0 ${double_colon_index} package)

    # Extract 'target'
    math(EXPR target_index "${double_colon_index} + 2")
    string(SUBSTRING "${package_target}" ${target_index} -1 target)

    # If 'target' is empty, set it as 'package'
    if("${target}" STREQUAL "")
      set(target "${package}")
    endif()
  else()
    # If '::' is not found, set 'package' and 'target' as what is given
    set(package "${package_target}")
    set(target "${package_target}")
  endif()

  if(NOT TARGET ${target})
    find_package(${package} CONFIG QUIET)
    if(NOT ${package}_FOUND)
      include(FetchContent)

      if(NOT ARG_QUIET)
        message(CHECK_START "${PROJECT_NAME}: Fetching the ${package} library")
      endif()

      FetchContent_Declare(
        lib${package}
        GIT_REPOSITORY git@gitlab.dkrz.de:icon-libraries/lib${package}.git
        GIT_TAG ${tag})

      set(save_BUILD_TESTING "${BUILD_TESTING}")
      set(BUILD_TESTING
          OFF
          CACHE BOOL "Do not need to build testing for ${package}" FORCE)

      FetchContent_MakeAvailable(lib${package})
      if(NOT ARG_QUIET)
        message(CHECK_PASS "done")
      endif()

      set(BUILD_TESTING
          "${save_BUILD_TESTING}"
          CACHE BOOL "Copy back the original value" FORCE)

    endif()
  endif()

endfunction()
