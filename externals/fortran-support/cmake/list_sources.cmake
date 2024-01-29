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

# list_sources(<var>
#              [DIRECTORY <dir>]
#              [INCLUDE_REGEX <include_regex>]
#              [EXCLUDE_GENERATED])
# ------------------------------------------------------------------------------
# Sets <var> to a list of absolute paths to the source files of all targets in
# all subdirectories of <dir> (defaults to the current source directory). The
# duplicates are excluded.
#
# If the INCLUDE_REGEX argument is provided, paths that do not match
# <include_regex> are excluded from the result.
#
# If the EXCLUDE_GENERATED argument is provided, paths to generated source files
# are excluded from the result.
#
function(list_sources var)
  cmake_parse_arguments(PARSE_ARGV 1 ARG
    "EXCLUDE_GENERATED"
    "DIRECTORY;INCLUDE_REGEX"
    ""
  )

  if(ARG_EXCLUDE_GENERATED AND "${CMAKE_VERSION}" VERSION_LESS "3.18")
    message(AUTHOR_WARNING
      "The generated source files can be automatically excluded only with CMake 3.18 or newer.")
    set(ARG_EXCLUDE_GENERATED FALSE)
  endif()

  function(list_sources_recurse var dir)
    get_directory_property(dir_path DIRECTORY ${dir} SOURCE_DIR)
    get_directory_property(targets DIRECTORY ${dir} BUILDSYSTEM_TARGETS)
    set(result)
    foreach(target ${targets})
      if("${CMAKE_VERSION}" VERSION_LESS "3.19")
        # Older versions of CMake fail with an error if the SOURCES property of
        # an interface library is requested:
        get_target_property(type ${target} TYPE)
        if("${type}" STREQUAL "INTERFACE_LIBRARY")
          continue()
        endif()
      endif()
      get_target_property(sources ${target} SOURCES)
      if(sources)
        foreach(source ${sources})
          if(${ARG_EXCLUDE_GENERATED})
            get_source_file_property(source_generated ${source} DIRECTORY ${dir} GENERATED)
            if(source_generated)
              continue()
            endif()
          endif()
          get_filename_component(source ${source} ABSOLUTE BASE_DIR ${dir_path})
          list(APPEND result ${source})
        endforeach()
      endif()
    endforeach()
    get_directory_property(subdirs DIRECTORY ${dir} SUBDIRECTORIES)
    foreach(subdir ${subdirs})
      list_sources_recurse(subdir_sources ${subdir} ${ARG_EXCLUDE_GENERATED})
      list(APPEND result ${subdir_sources})
    endforeach()
    set(${var} ${result} PARENT_SCOPE)
  endfunction()

  if(NOT ARG_DIRECTORY)
    set(ARG_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()

  list_sources_recurse(result ${ARG_DIRECTORY})

  if(ARG_INCLUDE_REGEX)
    list(FILTER result INCLUDE REGEX ${ARG_INCLUDE_REGEX})
  endif()

  list(REMOVE_DUPLICATES result)
  set(${var} ${result} PARENT_SCOPE)
endfunction()
