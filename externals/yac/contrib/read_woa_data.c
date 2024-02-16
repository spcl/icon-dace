/**
 * @file read_woa_data.c
 *
 * @copyright Copyright  (C)  2022 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <stdio.h>
#include <string.h>
#include <netcdf.h>

#include "read_woa_data.h"
#include "utils.h"
#include "io_utils.h"

/* ------------------------------------------------ */

int open_woa_output(char const * input_file) {

  int ncid;
  yac_nc_open(input_file, NC_NOWRITE, &ncid);
  return ncid;
}

/* ------------------------------------------------ */

void close_woa_output(int ncid) {
  HANDLE_ERROR(nc_close(ncid));
}

/* ------------------------------------------------ */

void read_woa_dimensions ( int fileId, char const * fieldName,
                            struct fieldMetadata * fieldInfo ) {

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
  HANDLE_ERROR(nc_inq_vartype(fileId, varId, &varType));
  YAC_ASSERT_F(
    (varType == NC_FLOAT) || (varType == NC_DOUBLE),
    "ERROR(read_woa_dimensions): unsupported datatype for variable \"%s\" "
    "(has to be either NC_DOUBLE or NC_FLOAT)", fieldName)

  // get number of dimensions
  int varNdims;
  HANDLE_ERROR(nc_inq_varndims(fileId, varId, &varNdims));
  YAC_ASSERT_F(
    (varNdims > 0) && (varNdims <= 4),
    "ERROR(read_woa_dimensions): invalid number of dimensions for "
    "variable \"%s\" (has to be between 1 and 4)", fieldName)

  // get dimensions of variable
  int varDimids[NC_MAX_VAR_DIMS];   // dimension NetCDF IDs
  HANDLE_ERROR(nc_inq_vardimid(fileId, varId, varDimids));

  // check dimensions
  for (int i = 0; i < varNdims; ++i) {

    // get details of current dimension
    char dim_name[NC_MAX_NAME + 1];
    size_t dim_len;
    HANDLE_ERROR(nc_inq_dim(fileId, varDimids[i], dim_name, &dim_len));

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
}

/* ------------------------------------------------ */

double * get_woa_memory(struct fieldMetadata fieldInfo) {

  return
    xmalloc(fieldInfo.nbrLatPoints * fieldInfo.nbrLonPoints * sizeof(double));
}

/* ------------------------------------------------ */

void free_woa_memory(double * data) {
  free(data);
}

/* ------------------------------------------------ */

void read_woa_timestep_level(
  int fileId, double * cell_data, struct fieldMetadata fieldInfo,
  int timeStep, int level) {

  size_t start[4] = {0,0,0,0};
  size_t count[4] = {1,1,1,1};

  /* ... extract one level */
  if ( fieldInfo.levelDimIdx > -1 ) {
    YAC_ASSERT_F(
      (level > 0) && (level <= fieldInfo.nbrLevels),
      "ERROR(read_woa_timestep_level): "
      "invalid level (has to be between 1 and %d)",
      (int)(fieldInfo.nbrLevels));
    start[fieldInfo.levelDimIdx] = level-1;
    count[fieldInfo.levelDimIdx] = 1;
  }

  /* ... extract one timestep */
  if ( fieldInfo.timeDimIdx > -1 ) {
    YAC_ASSERT_F(
      (timeStep > 0) && (timeStep <= fieldInfo.nbrTimeSteps),
      "ERROR(read_woa_timestep_level): "
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
  HANDLE_ERROR(
    nc_get_vara_double(fileId, fieldInfo.varId, start, count, cell_data));
}
