// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <string.h>

#ifdef YAC_NETCDF_ENABLED
#include <netcdf.h>
#endif

#include "read_woa_data.h"
#include "utils_common.h"
#include "io_utils.h"

/* ------------------------------------------------ */

int yac_open_woa_output(char const * input_file) {

#ifndef YAC_NETCDF_ENABLED

   UNUSED(input_file);
   die(
     "ERROR(yac_open_woa_output): "
     "YAC is built without the NetCDF support");
  return -1;
#else

  int ncid;
  yac_nc_open(input_file, NC_NOWRITE, &ncid);
  return ncid;
#endif
}

/* ------------------------------------------------ */

void yac_close_woa_output(int ncid) {

#ifndef YAC_NETCDF_ENABLED

   UNUSED(ncid);
   die(
     "ERROR(yac_close_woa_output): "
     "YAC is built without the NetCDF support");
#else

  YAC_HANDLE_ERROR(nc_close(ncid));
#endif
}

/* ------------------------------------------------ */

void yac_read_woa_dimensions ( int fileId, char const * fieldName,
                               struct yac_fieldMetadata * fieldInfo ) {

#ifndef YAC_NETCDF_ENABLED

   UNUSED(fileId);
   UNUSED(fieldName);
   UNUSED(fieldInfo);
   die(
     "ERROR(yac_read_woa_dimensions): "
     "YAC is built without the NetCDF support");
#else

  // set defaults
  fieldInfo->nbrLevels    = 1;
  fieldInfo->nbrTimeSteps = 1;
  fieldInfo->nbrLatPoints = 1;
  fieldInfo->nbrLonPoints = 1;
  fieldInfo->timeDimIdx   = -1;
  fieldInfo->latDimIdx    = -1;
  fieldInfo->lonDimIdx    = -1;
  fieldInfo->levelDimIdx  = -1;

  // get id of variable
  int varId;
  yac_nc_inq_varid(fileId, fieldName, &varId);
  fieldInfo->varId   = varId;

  // check type of variable
  nc_type varType;
  YAC_HANDLE_ERROR(nc_inq_vartype(fileId, varId, &varType));
  YAC_ASSERT_F(
    (varType == NC_FLOAT) || (varType == NC_DOUBLE),
    "ERROR(yac_read_woa_dimensions): "
    "unsupported datatype for variable \"%s\" "
    "(has to be either NC_DOUBLE or NC_FLOAT)", fieldName)

  // get number of dimensions
  int varNdims;
  YAC_HANDLE_ERROR(nc_inq_varndims(fileId, varId, &varNdims));
  YAC_ASSERT_F(
    (varNdims > 0) && (varNdims <= 4),
    "ERROR(read_woa_dimensions): invalid number of dimensions for "
    "variable \"%s\" (has to be between 1 and 4)", fieldName)

  // get dimensions of variable
  int varDimids[NC_MAX_VAR_DIMS];   // dimension NetCDF IDs
  YAC_HANDLE_ERROR(nc_inq_vardimid(fileId, varId, varDimids));

  // check dimensions
  for (int i = 0; i < varNdims; ++i) {

    // get details of current dimension
    char dim_name[NC_MAX_NAME + 1];
    size_t dim_len;
    YAC_HANDLE_ERROR(nc_inq_dim(fileId, varDimids[i], dim_name, &dim_len));

    size_t * count = NULL;
    int * idx = NULL;

    if (!strcmp(dim_name, "time")) {
      count = &(fieldInfo->nbrTimeSteps);
      idx = &(fieldInfo->timeDimIdx);
    } else if (!strcmp(dim_name, "depth")) {
      count = &(fieldInfo->nbrLevels);
      idx = &(fieldInfo->levelDimIdx);
    } else if (!strcmp(dim_name, "lon")) {
      count = &(fieldInfo->nbrLonPoints);
      idx = &(fieldInfo->lonDimIdx);
    } else if (!strcmp(dim_name, "lat")) {
      count = &(fieldInfo->nbrLatPoints);
      idx = &(fieldInfo->latDimIdx);
    }

    YAC_ASSERT_F(
      (count != NULL) && (idx != NULL),
      "ERROR(read_woa_dimensions): "
      "supported dimension \"%s\" of variable \"%s\"", dim_name, fieldName)

    *count = dim_len;
    *idx = i;
  }
#endif
}

/* ------------------------------------------------ */

double * yac_get_woa_memory(struct yac_fieldMetadata fieldInfo) {

  return
    xmalloc(fieldInfo.nbrLatPoints * fieldInfo.nbrLonPoints * sizeof(double));
}

/* ------------------------------------------------ */

void yac_free_woa_memory(double * data) {
  free(data);
}

/* ------------------------------------------------ */

void yac_read_woa_timestep_level(
  int fileId, double * cell_data, struct yac_fieldMetadata fieldInfo,
  int timeStep, int level) {

#ifndef YAC_NETCDF_ENABLED

   UNUSED(fileId);
   UNUSED(cell_data);
   UNUSED(fieldInfo);
   UNUSED(timeStep);
   UNUSED(level);
   die(
     "ERROR(yac_read_woa_timestep_level): "
     "YAC is built without the NetCDF support");
#else

  size_t start[4] = {0,0,0,0};
  size_t count[4] = {1,1,1,1};

  /* ... extract one level */
  if ( fieldInfo.levelDimIdx > -1 ) {
    YAC_ASSERT_F(
      ((size_t)level > 0) && ((size_t)level <= fieldInfo.nbrLevels),
      "ERROR(yac_read_woa_timestep_level): "
      "invalid level (has to be between 1 and %d)",
      (int)(fieldInfo.nbrLevels));
    start[fieldInfo.levelDimIdx] = level-1;
    count[fieldInfo.levelDimIdx] = 1;
  }

  /* ... extract one timestep */
  if ( fieldInfo.timeDimIdx > -1 ) {
    YAC_ASSERT_F(
      ((size_t)timeStep > 0) && ((size_t)timeStep <= fieldInfo.nbrTimeSteps),
      "ERROR(yac_read_woa_timestep_level): "
      "invalid time step (has to be between 1 and %d)",
      (int)(fieldInfo.nbrTimeSteps));
    start[fieldInfo.timeDimIdx] = timeStep-1;
    count[fieldInfo.timeDimIdx] = 1;
  }

  /* ... extract one horizontal data set */
  if ( fieldInfo.latDimIdx > -1 ) {
    count[fieldInfo.latDimIdx] = fieldInfo.nbrLatPoints;
  }
  if ( fieldInfo.lonDimIdx > -1 ) {
    count[fieldInfo.lonDimIdx] = fieldInfo.nbrLonPoints;
  }

  /* ... get data from netcdf file */
  YAC_HANDLE_ERROR(
    nc_get_vara_double(fileId, fieldInfo.varId, start, count, cell_data));
#endif
}
