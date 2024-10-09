// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef WOA_DATA_H
#define WOA_DATA_H

#include <stdlib.h>

// YAC PUBLIC HEADER START

/** \example test_read_woa_data.c
 * This contains examples for read_woa_data.
 */

/** \file read_woa_data.h
  * \brief general routines for reading WOA output NetCDF files
  *
  * These routines should be called in a specific order:
  *
  * -# \ref yac_open_woa_output
  * -# \ref yac_read_woa_dimensions
  * -# \ref yac_get_woa_memory
  * -# \ref yac_read_woa_timestep_level
  * -# \ref yac_free_woa_memory
  * -# \ref yac_close_woa_output
  *
  * \remark These routines are adapted to read the current WOA NetCDF output.
  * \remark It is assumed that the horizonal dimension is available in a 1d array cell,
  * \remark the vertical information can be accessed via depth, and time via time.
 **/



struct yac_fieldMetadata {
  int varId;            //!< NetCDF variable ID
  size_t nbrTimeSteps;  //!< number of timesteps contained in the NetCDF file
  size_t nbrLevels;     //!< number of vertical levels/layers contained in the NetCDF file
  size_t nbrLatPoints;  //!< number of latitude cells contained in the NetCDF file
  size_t nbrLonPoints;  //!< number of longitude cells contained in the NetCDF file
  int timeDimIdx;       //!< dimension index from NetCDF containing the time
  int levelDimIdx;      //!< dimension index from NetCDF containing the vertical
  int latDimIdx;        //!< dimension index from NetCDF containing the latitude
  int lonDimIdx;        //!< dimension index from NetCDF containing the longitude
};

/**
 * To open a NetCDF file
 * @param[in] input_file file name of the NetCDF input file including the file extension
 * The function returns an integer value (the NetCDF ID) which has to be used by subsequent calls
 * when this file is accessed.
 */
int yac_open_woa_output ( char const * input_file );

/**
 * To close a NetCDF file
 * @param[in] fileId NetCDF file ID as it was returned by yac_open_woa_output
 */
void yac_close_woa_output ( int fileId );

/**
 * To read in the dimensions for arrays stored in the file
 * @param[in]  fileId    NetCDF file ID as it was returned by yac_open_woa_output
 * @param[in]  fieldName name of the array that shall be read in
 * @param[out] fieldInfo metadata information for the array
 */
void yac_read_woa_dimensions ( int fileId, char const * fieldName,
                               struct yac_fieldMetadata * fieldInfo );

/**
 * To get the appropriate memory for the data to be read in
 * @param[in]  fieldInfo metadata information for the array
 * @return pointer to memory that is big enough to hold a single level of a
 *         single time step of the field associated to fieldInfo
 */
double * yac_get_woa_memory(struct yac_fieldMetadata fieldInfo );

/**
 * To release allocated memory
 * @param[in]  data memory for data from NetCDF file
 */
void yac_free_woa_memory(double * data);

/**
 * To read in on timestep of one particular field
 * @param[in]  fileId    NetCDF file ID as it was returned by yac_open_woa_output
 * @param[out] cell_data one level of the requested field
 * @param[out] fieldInfo metadata information for the array
 * @param[in]  timestep  time step to be read from the NetCDF file
 * @param[in]  level     vertical level to be read from the NetCDF file
 */
void yac_read_woa_timestep_level(
  int fileId, double * cell_data, struct yac_fieldMetadata fieldInfo,
  int timestep, int level);

// YAC PUBLIC HEADER STOP

#endif // WOA_DATA_H

