/**
 * @file read_woa_data.h
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
#ifndef WOA_DATA_H
#define WOA_DATA_H

#include <stdlib.h>

/** \example test_read_woa_data.c
 * This contains examples for read_woa_data.
 */

/** \file read_woa_data.h
  * \brief general routines for reading WOA output NetCDF files
  *
  * These routines should be called in a specific order:
  *
  * -# \ref open_woa_output
  * -# \ref read_woa_dimensions
  * -# \ref get_woa_memory
  * -# \ref read_woa_timestep_level
  * -# \ref free_woa_memory
  * -# \ref close_woa_output
  *
  * \remark These routines are adapted to read the current WOA NetCDF output.
  * \remark It is assumed that the horizonal dimension is available in a 1d array cell,
  * \remark the vertical information can be accessed via depth, and time via time.
 **/



struct fieldMetadata {
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
int open_woa_output ( char const * input_file );

/**
 * To close a NetCDF file
 * @param[in] fileId NetCDF file ID as it was returned by open_woa_output
 */
void close_woa_output ( int fileId );

/**
 * To read in the dimensions for arrays stored in the file
 * @param[in]  fileId    NetCDF file ID as it was returned by open_woa_output
 * @param[in]  fieldName name of the array that shall be read in
 * @param[out] fieldInfo metadata information for the array
 */
void read_woa_dimensions ( int fileId, char const * fieldName,
                           struct fieldMetadata * fieldInfo );

/**
 * To get the appropriate memory for the data to be read in
 * @param[in]  fieldInfo metadata information for the array
 * @return pointer to memory that is big enough to hold a single level of a
 *         single time step of the field associated to fieldInfo
 */
double * get_woa_memory(struct fieldMetadata fieldInfo );

/**
 * To release allocated memory
 * @param[in]  data memory for data from NetCDF file
 */
void free_woa_memory(double * data);

/**
 * To read in on timestep of one particular field
 * @param[in]  fileId    NetCDF file ID as it was returned by open_woa_output
 * @param[out] cell_data one level of the requested field
 * @param[out] fieldInfo metadata information for the array
 * @param[in]  timestep  time step to be read from the NetCDF file
 * @param[in]  level     vertical level to be read from the NetCDF file
 */
void read_woa_timestep_level(
  int fileId, double * cell_data, struct fieldMetadata fieldInfo,
  int timestep, int level);

#endif // WOA_DATA_H

