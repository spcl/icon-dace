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

#include <netcdf.h>
#include <string.h>

#include "grid.h"
#include "utils.h"
#include "io_utils.h"

static size_t check_dimension(int ncid, int varids[2]) {

  for (int i = 0; i < 2; ++i) {
    int ndims;
    HANDLE_ERROR(nc_inq_varndims(ncid, varids[i], &ndims));
    YAC_ASSERT(
      ndims == 1,
      "ERROR(check_dimension): coordinate array has more than one dimension")
  }

  int dimids[2];
  for (int i = 0; i < 2; ++i)
    HANDLE_ERROR(nc_inq_vardimid(ncid, varids[i], &(dimids[i])));

  YAC_ASSERT(
    dimids[0] == dimids[1],
    "ERROR(check_dimension): "
    "lon lat coordinate arrays have differing dimensions")

  size_t dimlen;
  HANDLE_ERROR(nc_inq_dimlen(ncid, dimids[0], &dimlen));
  return dimlen;
}

int check_coord_units(int ncid, int varid) {

  int is_degree = 0;
  nc_type att_type;
  size_t att_len;
  int status = nc_inq_att(ncid, varid, "units", &att_type, &att_len);
  // if the status is not "attribute not found"
  if (status != NC_ENOTATT) {
    HANDLE_ERROR(status);
    // if the attribute is not a string or too long
    YAC_ASSERT(
      (att_type == NC_CHAR) && (att_len <= 8),
      "ERROR(read_coord): invalid units type or len")
    char units[8];
    memset(units, 0, 8 * sizeof(units[0]));
    HANDLE_ERROR(nc_get_att_text(ncid, varid, "units", units));
    is_degree = !strcmp(units, "degree");
    YAC_ASSERT(
      is_degree || !strcmp(units, "radian"),
      "ERROR(read_coord): unsupported units type")
  }
  return is_degree;
}

static double * read_coord(int ncid, int varid, size_t varlen) {

  int is_degree = check_coord_units(ncid, varid);

  double * coord = xmalloc(varlen * sizeof(*coord));
  HANDLE_ERROR(nc_get_var_double (ncid, varid, coord));

  // convert to radiant if necessary
  if (is_degree) for (size_t i = 0; i < varlen; ++i) coord[i] *= YAC_RAD;

  return coord;
}

void read_coords(
  int ncid, char const * lon_name, char const * lat_name,
  double ** lon, double ** lat, size_t * len) {

  int vlonid, vlatid;
  yac_nc_inq_varid(ncid, lon_name, &vlonid);
  yac_nc_inq_varid(ncid, lat_name, &vlatid);

  size_t varlen = (*len = check_dimension(ncid, (int[]){vlonid, vlatid}));
  *lon = read_coord(ncid, vlonid, varlen);
  *lat = read_coord(ncid, vlatid, varlen);
}
