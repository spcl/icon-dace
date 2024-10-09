// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include <netcdf.h>

#include "yac_utils.h"

// redefine YAC assert macros
#undef YAC_ASSERT
#undef YAC_ASSERT_F

#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

static char const * cmd;
#define STR_USAGE "Usage: %s -n grid_name -f grid_filename -o vtk_filename\n"

#define YAC_ASSERT(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, cmd); \
      exit(EXIT_FAILURE); \
    } \
  }

#define YAC_ASSERT_F(exp, format, ...) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: " format "\n" STR_USAGE, __VA_ARGS__, cmd); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(int argc, char ** argv,
                            char const ** grid_filename,
                            char const ** grid_name,
                            char const ** vtk_filename);
static void read_grid_file(
  char const * grid_filename, char const * grid_name,
  size_t * nv, size_t * nc,
  double ** clon, double ** clat, double ** vlon, double ** vlat,
  int ** gid, int ** cmk, int ** rnk);
static void write_vtk_file(
  char const * vtk_filename, char const * grid_name,
  size_t nv, size_t nc,
  double * clon, double * clat, double * vlon, double * vlat,
  int * gid, int * cmk, int * rnk);

int main(int argc, char ** argv) {

  cmd = argv[0];

  char const * grid_name, * grid_filename, * vtk_filename;
  parse_arguments(argc, argv, &grid_filename, &grid_name, &vtk_filename);

  size_t nv, nc;
  double * clon, * clat, * lon, * lat;
  int * gid, * cmk, * rnk;

  read_grid_file(
    grid_filename, grid_name,
    &nv, &nc, &clon, &clat, &lon, &lat, &gid, &cmk, &rnk);

  write_vtk_file(
    vtk_filename, grid_name,
    nv, nc, clon, clat, lon, lat, gid, cmk, rnk);

  free(clon);
  free(clat);
  free(lon);
  free(lat);
  free(gid);
  free(cmk);
  free(rnk);
}

static double * read_coords(int ncid, int varid, size_t varlen) {

  int is_degree = yac_check_coord_units(ncid, varid);

  double * coord = malloc(varlen * sizeof(*coord));
  YAC_HANDLE_ERROR(nc_get_var_double (ncid, varid, coord));

  // convert to radiant if necessary
  if (is_degree) for (size_t i = 0; i < varlen; ++i) coord[i] *= YAC_RAD;

  return coord;
}

static int * read_int(int ncid, int varid, size_t varlen) {

  int * array = malloc(varlen * sizeof(*array));
  YAC_HANDLE_ERROR(nc_get_var_int(ncid, varid, array));
  return array;
}

static void read_grid_file(
  char const * grid_filename, char const * grid_name,
  size_t * nv, size_t * nc,
  double ** clon, double ** clat, double ** lon, double ** lat,
  int ** gid, int ** cmk, int ** rnk) {

  int status;

  size_t grid_name_len = strlen(grid_name) + 1;
  char nv_dim_name[3 + grid_name_len];
  char nc_dim_name[3 + grid_name_len];
  char cla_var_name[4 + grid_name_len];
  char clo_var_name[4 + grid_name_len];
  char lat_var_name[4 + grid_name_len];
  char lon_var_name[4 + grid_name_len];
  char gid_var_name[4 + grid_name_len];
  char cmk_var_name[4 + grid_name_len];
  char rnk_var_name[4 + grid_name_len];

  snprintf(nv_dim_name, 3 + grid_name_len, "nv_%s", grid_name);
  snprintf(nc_dim_name, 3 + grid_name_len, "nc_%s", grid_name);
  snprintf(cla_var_name, 4 + grid_name_len, "%s.cla", grid_name);
  snprintf(clo_var_name, 4 + grid_name_len, "%s.clo", grid_name);
  snprintf(lat_var_name, 4 + grid_name_len, "%s.lat", grid_name);
  snprintf(lon_var_name, 4 + grid_name_len, "%s.lon", grid_name);
  snprintf(gid_var_name, 4 + grid_name_len, "%s.gid", grid_name);
  snprintf(cmk_var_name, 4 + grid_name_len, "%s.cmk", grid_name);
  snprintf(rnk_var_name, 4 + grid_name_len, "%s.rnk", grid_name);

  YAC_ASSERT_F(
    yac_file_exists(grid_filename), "File %s does not exist.", grid_filename)

  int ncid, dimid, var_id;
  yac_nc_open(grid_filename, NC_NOWRITE, &ncid);

  yac_nc_inq_dimid(ncid, nv_dim_name, &dimid);
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, nv));
  yac_nc_inq_dimid(ncid, nc_dim_name, &dimid);
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, nc));

  yac_nc_inq_varid(ncid, cla_var_name, &var_id);
  *clat = read_coords(ncid, var_id, *nv * *nc);
  yac_nc_inq_varid(ncid, clo_var_name, &var_id);
  *clon = read_coords(ncid, var_id, *nv * *nc);

  status = nc_inq_varid(ncid, lat_var_name, &var_id);
  if (status == NC_NOERR) {
    *lat = read_coords(ncid, var_id, *nc);
  } else if (status == NC_ENOTVAR) {
    *lat = NULL;
    status = NC_NOERR;
  }
  YAC_HANDLE_ERROR(status);
  status = nc_inq_varid(ncid, lon_var_name, &var_id);
  if (status == NC_NOERR) {
    *lon = read_coords(ncid, var_id, *nc);
  } else if (status == NC_ENOTVAR) {
    *lon = NULL;
    status = NC_NOERR;
  }
  YAC_HANDLE_ERROR(status);

  status = nc_inq_varid(ncid, gid_var_name, &var_id);
  if (status == NC_NOERR) {
    *gid = read_int(ncid, var_id, *nc);
  } else if (status == NC_ENOTVAR) {
    *gid = NULL;
    status = NC_NOERR;
  }
  YAC_HANDLE_ERROR(status);

  status = nc_inq_varid(ncid, cmk_var_name, &var_id);
  if (status == NC_NOERR) {
    *cmk = read_int(ncid, var_id, *nc);
  } else if (status == NC_ENOTVAR) {
    *cmk = NULL;
    status = NC_NOERR;
  }
  YAC_HANDLE_ERROR(status);

  YAC_HANDLE_ERROR(nc_inq_varid(ncid, rnk_var_name, &var_id));
  *rnk = read_int(ncid, var_id, *nc);

  YAC_HANDLE_ERROR(nc_close(ncid));
}

static void LLtoXYZ(double lon, double lat, double p_out[]) {

   if ((fabs(lon) > 8.0 * M_PI) || (fabs(lat) > 2.0 * M_PI)) {
     p_out[0] = 0.0;
     p_out[1] = 0.0;
     p_out[2] = 0.0;
     return;
   }

   while (lon < -M_PI) lon += 2.0 * M_PI;
   while (lon >= M_PI) lon -= 2.0 * M_PI;

   double cos_lat = cos(lat);
   p_out[0] = cos_lat * cos(lon);
   p_out[1] = cos_lat * sin(lon);
   p_out[2] = sin(lat);
}

static void write_vtk_file(
  char const * vtk_filename, char const * grid_name,
  size_t nv, size_t nc,
  double * clon, double * clat, double * lon, double * lat,
  int * gid, int * cmk, int * rnk) {

  yac_coordinate_pointer points =
    malloc((nv * nc + ((lon && lat)?nc:0)) * sizeof(*points));

  for (size_t i = 0; i < nv * nc; ++i) LLtoXYZ(clon[i], clat[i], points[i]);
  if (lon && lat)
    for (size_t i = 0; i < nc; ++i)
      LLtoXYZ(lon[i], lat[i], points[nv * nc + i]);

  unsigned * num_points_per_polygon =
    malloc(nc * sizeof(*num_points_per_polygon));
  for (size_t i = 0; i < nc; ++i) num_points_per_polygon[i] = (unsigned)nv;

  unsigned * polygon_data = malloc(nc * nv * sizeof(*polygon_data));
  for (size_t i = 0; i < nv * nc; ++i) polygon_data[i] = (unsigned)i;

  YAC_VTK_FILE * file = yac_vtk_open(vtk_filename, grid_name);
  yac_vtk_write_point_data(file, &(points[0][0]), nv * nc + ((lon && lat)?nc:0));
  yac_vtk_write_cell_data(
    file, polygon_data, num_points_per_polygon, nc);
  if (gid) yac_vtk_write_cell_scalars_int(file, gid, nc, "gid");
  if (cmk) yac_vtk_write_cell_scalars_int(file, cmk, nc, "cmk");
  if (rnk) yac_vtk_write_cell_scalars_int(file, rnk, nc, "rnk");
  yac_vtk_close(file);

  free(polygon_data);
  free(num_points_per_polygon);
  free(points);
}

static void parse_arguments(int argc, char ** argv,
                            char const ** grid_filename,
                            char const ** grid_name,
                            char const ** vtk_filename) {

  *grid_name = NULL;
  *grid_filename = NULL;
  *vtk_filename = NULL;;

  int opt;
  while ((opt = getopt(argc, argv, "n:f:o:")) != -1) {
    YAC_ASSERT(
      (opt == 'n') ||
      (opt == 'f') ||
      (opt == 'o'), "invalid command argument")
    switch (opt) {
      default:
      case 'n':
        *grid_name = optarg;
        break;
      case 'f':
        *grid_filename = optarg;
        break;
      case 'o':
        *vtk_filename = optarg;
        break;
    }
  }
  YAC_ASSERT_F(optind >= argc, "non-option ARGV-element: \"%s\"", argv[optind])
  YAC_ASSERT(argc != 1, "too few arguments")
  YAC_ASSERT(*grid_name != NULL,  "grid_name argument is missing")
  YAC_ASSERT(*grid_filename != NULL,  "grid_filename argument is missing")
  YAC_ASSERT(*vtk_filename != NULL, "vtk_filename argument is missing")
}
