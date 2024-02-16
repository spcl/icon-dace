/**
 * @file weights2vtk.c
 *
 * @copyright Copyright  (C)  2015 DKRZ, MPI-M
 *
 * @version 1.0
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

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include <netcdf.h>

#include "read_icon_grid.h"
#include "read_mpiom_grid.h"
#include "read_scrip_grid.h"
#include "utils.h"
#include "geometry.h"
#include "vtk_output.h"
#include "generate_cubed_sphere.h"
#include "io_utils.h"

#define YAC_WEIGHT_FILE_VERSION_1_0_STRING "yac weight file 1.0"

// redefine YAC assert macros
#undef YAC_ASSERT
#undef YAC_ASSERT_F

static char const * cmd;
#define STR_USAGE \
  "Usage: %s -S src_grid_type -T tgt_grid_type " \
  "-s src_filename -t tgt_filename -w weight_filename " \
  "-o vtk_filename\n\n" \
  "   grid_type can have the following values:\n" \
  "      'c', 'C': cubed sphere grid (instead of a " \
  "filename provide N, where\n" \
  "                N = sqrt(n/6), with n being the " \
  "total number of cells\n" \
  "      'm', 'M': curvilinear grid\n" \
  "      'i', 'I': unstructured grid\n" \
  "      'g', 'G': gaussian grid (instead of a filename provide the grid\n" \
  "                configuration \"x1,y1,x2,y2,nx,ny\", where:\n"\
  "                x1,y1: are the coordinates of the first grid corner\n"\
  "                x2,y2: are the coordinates of the last grid corner\n"\
  "                nx,ny: are the number of cells in each direction\n"\
  "                example: 0.0,-90.0,360.0,90.0,360,180)\n"\
  "      's', 'S': unstructured grid in scrip file format\n" \
  "                (in addition to the grid filename also provide \n" \
  "                the mask filename and the grid name \"g,m,n\", where:\n" \
  "                g: grid filename\n" \
  "                m: mask filename\n" \
  "                n: grid name)\n" \
  "      - a lower case grid type indicates, that the global ids\n" \
  "        of the respective grid were \"0\"-based in the model\n" \
  "        component that generated the weight file\n" \
  "      - an upper case grid type indicates \"1\"-based global\n" \
  "        ids\n" \
  "\n" \
  "Example:\n" \
  "  Configuration:\n" \
  "    source grid type: unstructured\n" \
  "    source grid file: src_grid.nc\n" \
  "    target grid type: cubed sphere\n" \
  "    target grid file: 100 (N instead of filename)\n" \
  "    weight file name: weight_file.nc\n" \
  "    vtk file name:    weights.vtk\n" \
  "  Command:\n" \
  "     %s -S i -T c -s src_grid.nc -t 100 -w weight_file.nc " \
  "-o weights.vtk\n"

#define YAC_ASSERT(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, cmd, cmd); \
      exit(EXIT_FAILURE); \
    } \
  }

#define YAC_ASSERT_F(exp, format, ...) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: " format "\n" STR_USAGE, __VA_ARGS__, cmd, cmd); \
      exit(EXIT_FAILURE); \
    } \
  }

/** \example test_weights2vtk.c
 * A test for weights2vtk utility program.
 */

struct link_data {
  int src_idx, tgt_idx, points_idx;
  double weight;
};

enum grid_type {
  NONE,
  CUBE,
  CURVE,
  UNSTRUCT,
  GAUSS,
  SCRIP,
};

struct grid_config {
  union {
    struct {
      unsigned n;
    } cube;
    struct {
      char const * filename;
    } curve, unstruct;
    struct {
      double corners[2][2];
      size_t num_cells[2];
    } gauss;
    struct {
      char const * grid_filename;
      char const * mask_filename;
      char const * grid_name;
    } scrip;
  } config;
  int address_offset;
  enum grid_type type;
};

static void parse_arguments(int argc, char ** argv,
                            struct grid_config * src_grid_config,
                            struct grid_config * tgt_grid_config,
                            char const ** weight_filename,
                            char const ** vtk_filename);
static void read_link_data(char const * weight_filename, int src_address_offset,
                           int tgt_address_offset, struct link_data ** links,
                           unsigned * num_links, enum yac_location ** src_locations,
                           enum yac_location * tgt_location, unsigned * max_src_idx,
                           unsigned * max_tgt_idx);
static void write_data_to_file(char const * filename,
                               struct basic_grid_data src_grid,
                               struct basic_grid_data tgt_grid,
                               struct link_data * links,
                               unsigned num_links,
                               enum yac_location * src_locations,
                               enum yac_location tgt_location);
static struct basic_grid_data create_grid(struct grid_config grid_config);

int main(int argc, char ** argv) {

  cmd = argv[0];

  struct grid_config src_grid_config, tgt_grid_config;
  char const * weight_filename, * vtk_filename;

  parse_arguments(argc, argv, &src_grid_config, &tgt_grid_config,
                  &weight_filename, &vtk_filename);

  //-------------------------------------
  // read/generate source and target grid
  //-------------------------------------
  struct basic_grid_data src_grid = create_grid(src_grid_config);
  struct basic_grid_data tgt_grid = create_grid(tgt_grid_config);

  //-----------------
  // read weight file
  //-----------------
  YAC_ASSERT_F(
    yac_file_exists(weight_filename),
    "File %s does not exist.", weight_filename)

  struct link_data * links;
  unsigned num_links, max_src_idx, max_tgt_idx;
  enum yac_location * src_locations;
  enum yac_location tgt_location;
  read_link_data(weight_filename, src_grid_config.address_offset,
                 tgt_grid_config.address_offset, &links, &num_links,
                 &src_locations, &tgt_location, &max_src_idx, &max_tgt_idx);

  YAC_ASSERT(
    (max_src_idx < src_grid.num_cells) && (max_tgt_idx < tgt_grid.num_cells),
     "weight file does not match with source and target grid")

  //---------------
  // write vtk file
  //---------------
  write_data_to_file(
    vtk_filename, src_grid, tgt_grid, links, num_links,
    src_locations, tgt_location);

  //---------
  // clean up
  //---------

  free(src_locations);
  free(links);
  yac_basic_grid_data_free(tgt_grid);
  yac_basic_grid_data_free(src_grid);

  return EXIT_SUCCESS;
}

static void read_link_data(
  char const * weight_filename, int src_address_offset,
  int tgt_address_offset, struct link_data ** links,
  unsigned * num_links, enum yac_location ** src_locations,
  enum yac_location * tgt_location, unsigned * max_src_idx,
  unsigned * max_tgt_idx) {

  int ncid, dimid, var_id;
  yac_nc_open(weight_filename, NC_NOWRITE, &ncid);

  // global attributes
  size_t version_len;
  char * version;
  HANDLE_ERROR(nc_inq_attlen(ncid, NC_GLOBAL, "version", &version_len));
  version = xmalloc(version_len + 1);
  version[version_len] = '\0';
  HANDLE_ERROR(nc_get_att_text(ncid, NC_GLOBAL, "version", version));
  YAC_ASSERT(
    (strlen(YAC_WEIGHT_FILE_VERSION_1_0_STRING) == version_len) &&
    !strncmp(version, YAC_WEIGHT_FILE_VERSION_1_0_STRING, version_len),
    "version string from weight file is not \""
    YAC_WEIGHT_FILE_VERSION_1_0_STRING "\"");
  free(version);

  size_t str_link_len;
  char * str_link;
  var_id = NC_GLOBAL;
  HANDLE_ERROR(nc_inq_attlen(ncid, var_id, "contains_links", &str_link_len));
  str_link = xmalloc(str_link_len + 1);
  str_link[str_link_len] = '\0';
  HANDLE_ERROR(nc_get_att_text(ncid, var_id, "contains_links", str_link));

  int contains_links =
    (strlen("TRUE") == str_link_len) &&
    !strncmp("TRUE", str_link, str_link_len);
  YAC_ASSERT(
    contains_links ||
    ((strlen("FALSE") == str_link_len) &&
    !strncmp("FALSE", str_link, str_link_len)),
    "invalid global attribute contains_links")

  if (!contains_links) {
    *links = NULL;
    *num_links = 0;
    *src_locations = NULL;
    *max_src_idx = 0;
    *max_tgt_idx = 0;
    return;
  }

  size_t num_wgts = 0;
  HANDLE_ERROR(nc_inq_dimid(ncid, "num_wgts", &dimid));
  HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_wgts));
  YAC_ASSERT(
    num_wgts == 1, "YAC only supports num_wgts == 1")

  // get number of links from file
  size_t num_links_in_file = 0;
  HANDLE_ERROR(nc_inq_dimid(ncid, "num_links", &dimid));
  HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_links_in_file));
  YAC_ASSERT(num_links_in_file > 0, "no links defined")
  *num_links = num_links_in_file;
  *links = xmalloc(num_links_in_file * sizeof(**links));

  // get number of pointsets
  size_t num_pointsets = 0;
  HANDLE_ERROR(nc_inq_dimid(ncid, "num_src_fields", &dimid));
  HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_pointsets));
  YAC_ASSERT(num_pointsets > 0, "no point sets in file")

  // get max location string length from file
  size_t max_loc_str_len;
  HANDLE_ERROR(nc_inq_dimid(ncid, "max_loc_str_len", &dimid));
  HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &max_loc_str_len));
  YAC_ASSERT(
    max_loc_str_len == YAC_MAX_LOC_STR_LEN,
    "wrong max location string length in weight file")

  // get source locations
  *src_locations = xmalloc(num_pointsets * sizeof(**src_locations));
  yac_nc_inq_varid(ncid, "src_locations", &var_id);

  for (unsigned i = 0; i < num_pointsets; ++i) {

    char loc_str[YAC_MAX_LOC_STR_LEN];

    size_t str_start[2] = {i, 0};
    size_t str_count[2] = {1, YAC_MAX_LOC_STR_LEN};

    HANDLE_ERROR(nc_get_vara_text(ncid, var_id, str_start, str_count, loc_str));

    (*src_locations)[i] = yac_str2loc(loc_str);
  }

  // get target location
  yac_nc_inq_varid(ncid, "dst_location", &var_id);
  {
    char loc_str[YAC_MAX_LOC_STR_LEN];
    HANDLE_ERROR(nc_get_var_text(ncid, var_id, loc_str));
    *tgt_location = yac_str2loc(loc_str);
  }

  // get number of links per pointset
  unsigned * num_links_per_pointset =
    xmalloc(num_pointsets * sizeof(*num_links_per_pointset));
  yac_nc_inq_varid(ncid, "num_links_per_src_field", &var_id);
  HANDLE_ERROR(nc_get_var_uint(ncid, var_id, num_links_per_pointset));

  for (unsigned points_idx = 0, link_idx = 0; points_idx < num_pointsets;
       ++points_idx)
    for (unsigned i = 0; i < num_links_per_pointset[points_idx];
         ++i, ++link_idx)
      (*links)[link_idx].points_idx = points_idx;
  free(num_links_per_pointset);

  // get links
  int * address = xmalloc(num_links_in_file * sizeof(*address));
  yac_nc_inq_varid(ncid, "src_address", &var_id);
  HANDLE_ERROR(nc_get_var_int(ncid, var_id, address));
  *max_src_idx = 0;
  for (unsigned i = 0; i < num_links_in_file; ++i) {
    int idx = address[i] + src_address_offset;
    YAC_ASSERT_F(idx >= 0, "invalid src index (%d)", idx);
    if ((unsigned)idx > *max_src_idx) *max_src_idx = (unsigned)idx;
    (*links)[i].src_idx = (unsigned)idx;
  }

  yac_nc_inq_varid(ncid, "dst_address", &var_id);
  HANDLE_ERROR(nc_get_var_int(ncid, var_id, address));
  *max_tgt_idx = 0;
  for (unsigned i = 0; i < num_links_in_file; ++i) {
    int idx = address[i] + tgt_address_offset;
    YAC_ASSERT_F(idx >= 0, "invalid tgt index (%d)", idx);
    if ((unsigned)idx > *max_tgt_idx) *max_tgt_idx = (unsigned)idx;
    (*links)[i].tgt_idx = (unsigned)idx;
  }
  free(address);

  double * weights = xmalloc(num_links_in_file * sizeof(*weights));
  yac_nc_inq_varid(ncid, "remap_matrix", &var_id);
  HANDLE_ERROR(nc_get_var_double(ncid, var_id, weights));
  for (size_t i = 0; i < num_links_in_file; ++i)
    (*links)[i].weight = weights[i];
  free(weights);

  HANDLE_ERROR(nc_close(ncid));
}

static void get_cell_middle_point(
  struct basic_grid_data * grid, size_t cell_index, double * point) {

  int num_cell_corners = grid->num_vertices_per_cell[cell_index];
  size_t * cell_to_vertex = grid->cell_to_vertex +
                            grid->cell_to_vertex_offsets[cell_index];
  coordinate_pointer vertex_coordinates = grid->vertex_coordinates;

  point[0] = 0;
  point[1] = 0;
  point[2] = 0;

  for (int i = 0; i < num_cell_corners; ++i)
    for (int j = 0; j < 3; ++j)
      point[j] += vertex_coordinates[cell_to_vertex[i]][j];
  normalise_vector(point);
}

static void get_edge_middle_point(
  struct basic_grid_data * grid, size_t edge_index, double * point) {

  size_t * edge_to_vertex = grid->edge_to_vertex[edge_index];
  coordinate_pointer vertex_coordinates = grid->vertex_coordinates;

  point[0] = vertex_coordinates[edge_to_vertex[0]][0] +
             vertex_coordinates[edge_to_vertex[1]][0];
  point[1] = vertex_coordinates[edge_to_vertex[0]][1] +
             vertex_coordinates[edge_to_vertex[1]][1];
  point[2] = vertex_coordinates[edge_to_vertex[0]][2] +
             vertex_coordinates[edge_to_vertex[1]][2];

  normalise_vector(point);
}

static void get_point_coordinates(
  struct basic_grid_data * grid, size_t point_index,
  enum yac_location location, double * point) {

  YAC_ASSERT(
    (location == CELL) ||
    (location == CORNER) ||
    (location == EDGE), "unsupported point location")
  switch (location) {
    case(CELL):
      get_cell_middle_point(grid, point_index, point);
      break;
    case(CORNER):
      memcpy(point, grid->vertex_coordinates[point_index], 3 * sizeof(*point));
      break;
    default:
    case(EDGE):
      get_edge_middle_point(grid, point_index, point);
      break;
  };
}

static void get_link_xyz_coordinates(
  struct link_data * links, unsigned num_links,
  struct basic_grid_data * src_grid, struct basic_grid_data * tgt_grid,
  enum yac_location * src_locations, enum yac_location tgt_location,
  coordinate_pointer points) {

  for (unsigned i = 0; i < num_links; ++i) {

    get_point_coordinates(
      src_grid, (size_t)(links[i].src_idx), src_locations[links[i].points_idx],
      points[2*i+0]);
    get_point_coordinates(
      tgt_grid, (size_t)(links[i].tgt_idx), tgt_location, points[2*i+1]);
  }
}

static void get_grid_cell_data(
  struct basic_grid_data * grid, unsigned * cell_data, unsigned offset) {

   for (size_t i = 0, k = 0; i < grid->num_cells; ++i) {

      size_t * curr_cell_corners =
        grid->cell_to_vertex + grid->cell_to_vertex_offsets[i];
      int num_cell_corners = grid->num_vertices_per_cell[i];

      for (int j = 0; j < num_cell_corners; ++j)
         cell_data[k++] = curr_cell_corners[j] + offset;
   }
}

static void get_link_address_data(
  unsigned num_links, unsigned * polygon_data, unsigned offset) {

  for (unsigned i = 0; i < 2 * num_links; ++i)
    polygon_data[i] = i + offset;
}

static void write_data_to_file(char const * filename,
                               struct basic_grid_data src_grid,
                               struct basic_grid_data tgt_grid,
                               struct link_data * links, unsigned num_links,
                               enum yac_location * src_locations,
                               enum yac_location tgt_location) {

  //------------------------------------
  // generate cell data for the vtk file
  //------------------------------------

  size_t num_grid_corners[2] =
    {src_grid.num_vertices, tgt_grid.num_vertices};
  size_t total_num_points =
    num_grid_corners[0] + num_grid_corners[1] + 2 * (size_t)num_links;
  size_t num_polygons[3] =
    {src_grid.num_cells, tgt_grid.num_cells, num_links};
  size_t total_num_polygons =
    num_polygons[0] + num_polygons[1] + num_polygons[2];

  coordinate_pointer points = xmalloc(total_num_points * sizeof(*points));

  // get the point data of the grids
  memcpy(points, src_grid.vertex_coordinates,
         src_grid.num_vertices * sizeof(*points));
  memcpy(points + src_grid.num_vertices, tgt_grid.vertex_coordinates,
         tgt_grid.num_vertices * sizeof(*points));
  get_link_xyz_coordinates(
    links, num_links, &src_grid, &tgt_grid, src_locations, tgt_location,
    points + num_grid_corners[0] + num_grid_corners[1]);

  unsigned * num_points_per_polygon =
    xmalloc(total_num_polygons * sizeof(*num_points_per_polygon));
  unsigned num_points_per_polygon_sum[3] = {0, 0, 0};
  for (size_t i = 0; i < num_polygons[0]; ++i) {
    num_points_per_polygon_sum[0] +=
      (num_points_per_polygon[i] =
        (unsigned)(src_grid.num_vertices_per_cell[i]));
  }
  for (size_t i = 0; i < num_polygons[1]; ++i) {
   num_points_per_polygon_sum[1] +=
     (num_points_per_polygon[num_polygons[0] + i] =
       (unsigned)(tgt_grid.num_vertices_per_cell[i]));
  }
  for (size_t i = 0; i < num_polygons[2]; ++i)
    num_points_per_polygon[num_polygons[0] + num_polygons[1] + i] = 2;
  num_points_per_polygon_sum[2] = 2 * num_polygons[2];

  unsigned * polygon_data =
    xmalloc((num_points_per_polygon_sum[0] + num_points_per_polygon_sum[1] +
             num_points_per_polygon_sum[2]) * sizeof(*polygon_data));
  get_grid_cell_data(&src_grid, polygon_data, 0);
  get_grid_cell_data(&tgt_grid, polygon_data + num_points_per_polygon_sum[0],
                     num_grid_corners[0]);
  get_link_address_data(
    num_links, polygon_data + num_points_per_polygon_sum[0] +
    num_points_per_polygon_sum[1], num_grid_corners[0] + num_grid_corners[1]);

  //------------------------------------
  // generate scalar data
  //------------------------------------

  unsigned * polygon_type = xmalloc(total_num_polygons * sizeof(*polygon_type));
  for (size_t i = 0; i < num_polygons[0]; ++i) polygon_type[i] = 0;
  for (size_t i = 0; i < num_polygons[1]; ++i)
    polygon_type[num_polygons[0] + i] = 1;
  for (size_t i = 0; i < num_polygons[2]; ++i)
    polygon_type[num_polygons[0] + num_polygons[1] + i] = 2;

  double * weights = NULL;
  if (num_links > 0) {
    weights = xmalloc(total_num_polygons * sizeof(*weights));
    for (size_t i = 0; i < num_polygons[0] + num_polygons[1]; ++i)
      weights[i] = -1;
    for (size_t i = 0; i < num_polygons[2]; ++i)
     weights[i + num_polygons[0] + num_polygons[1]] = links[i].weight;
  }

  int * cell_ids = xmalloc(total_num_polygons * sizeof(*cell_ids));
  if (src_grid.cell_ids != NULL) {
    for (size_t i = 0; i < num_polygons[0]; ++i)
      cell_ids[i] = (int)(src_grid.cell_ids[i]);
  } else {
    for (size_t i = 0; i < num_polygons[0]; ++i) cell_ids[i] = (int)i;
  }
  if (tgt_grid.cell_ids != NULL) {
    for (size_t i = 0; i < num_polygons[1]; ++i)
      cell_ids[num_polygons[0] + i] = (int)(tgt_grid.cell_ids[i]);
  } else {
    for (size_t i = 0; i < num_polygons[1]; ++i)
      cell_ids[num_polygons[0] + i] = (int)i;
  }
  for (size_t i = 0; i < num_polygons[2]; ++i)
    cell_ids[num_polygons[0] + num_polygons[1] + i] = (int)i;

  int * vertex_ids = xmalloc(total_num_points * sizeof(*vertex_ids));
  if (src_grid.vertex_ids != NULL) {
    for (size_t i = 0; i < src_grid.num_vertices; ++i)
      vertex_ids[i] = (int)(src_grid.vertex_ids[i]);
  } else {
    for (size_t i = 0; i < src_grid.num_vertices; ++i)
      vertex_ids[i] = (int)i;
  }
  if (tgt_grid.vertex_ids != NULL) {
    for (size_t i = 0; i < tgt_grid.num_vertices; ++i)
      vertex_ids[src_grid.num_vertices + i] = (int)(tgt_grid.vertex_ids[i]);
  } else {
    for (size_t i = 0; i < tgt_grid.num_vertices; ++i)
      vertex_ids[src_grid.num_vertices + i] = (int)i;
  }
  for (size_t i = 0; i < 2 * (size_t)num_links; ++i)
    vertex_ids[src_grid.num_vertices + tgt_grid.num_vertices + i] =
      (int)(i >> 1);

  //------------------------------------
  // generate the actual vtk file
  //------------------------------------

  VTK_FILE * file = vtk_open(filename, "grid data");
  vtk_write_point_data(file, &(points[0][0]), total_num_points);
  vtk_write_cell_data(file, polygon_data, num_points_per_polygon, total_num_polygons);
  vtk_write_cell_scalars_int(file, cell_ids, total_num_polygons, "cell_ids");
  vtk_write_cell_scalars_uint(file, polygon_type, total_num_polygons, "polygon_type");
  if (num_links > 0)
    vtk_write_cell_scalars_double(file, weights, total_num_polygons, "weights");
  vtk_write_point_scalars_int(file, vertex_ids, total_num_points, "vertex_ids");
  vtk_close(file);

  //------------------------------------
  // some final cleanup
  //------------------------------------

  free(weights);
  free(polygon_type);
  free(cell_ids);
  free(vertex_ids);
  free(polygon_data);
  free(num_points_per_polygon);
  free(points);
}

static void interpret_grid_arg(
  struct grid_config * grid_config, char * arg, char * str) {

  YAC_ASSERT_F(arg != NULL, "-%c argument is missing", str[0])

  YAC_ASSERT_F(
    (grid_config->type == CUBE) ||
    (grid_config->type == CURVE) ||
    (grid_config->type == UNSTRUCT) ||
    (grid_config->type == GAUSS) ||
    (grid_config->type == SCRIP),
    "invalid %s grid type\n", str);

  switch (grid_config->type) {
    default:
    case CUBE:
      grid_config->config.cube.n = atoi(arg);
      YAC_ASSERT_F(
        grid_config->config.cube.n > 0,
        "invalid N for cubed sphere %s grid\n", str)
      break;
    case CURVE:
      grid_config->config.curve.filename = arg;
      break;
    case UNSTRUCT:
      grid_config->config.unstruct.filename = arg;
      break;
    case GAUSS:
      YAC_ASSERT_F(
        sscanf(
          arg, "%lf,%lf,%lf,%lf,%zu,%zu",
          &grid_config->config.gauss.corners[0][0],
          &grid_config->config.gauss.corners[0][1],
          &grid_config->config.gauss.corners[1][0],
          &grid_config->config.gauss.corners[1][1],
          &grid_config->config.gauss.num_cells[0],
          &grid_config->config.gauss.num_cells[1]) == 6,
        "invalid %s grid configuration (gauss grid)", str);
      YAC_ASSERT_F(
        (grid_config->config.gauss.num_cells[0] > 0) &&
        (grid_config->config.gauss.num_cells[1] > 0),
        "invalid %s grid configuration "
        "(gauss grid has invalid number of cells)", str)
      break;
    case SCRIP: {
      char * arg_copy = strdup(arg);
      grid_config->config.scrip.grid_filename = strtok(arg_copy, ",");
      grid_config->config.scrip.mask_filename = strtok(NULL, ",");
      grid_config->config.scrip.grid_name = strtok(NULL, ",");
      break;
    }
  }
}

void parse_arguments(int argc, char ** argv,
                     struct grid_config * src_grid_config,
                     struct grid_config * tgt_grid_config,
                     char const ** weight_filename,
                     char const ** vtk_filename) {

  src_grid_config->type = NONE;
  tgt_grid_config->type = NONE;
  *weight_filename = NULL;
  *vtk_filename = NULL;

  char * src_arg = NULL;
  char * tgt_arg = NULL;

  int opt;
  while ((opt = getopt(argc, argv, "S:T:s:t:w:o:")) != -1) {
    YAC_ASSERT(
      (opt == 'S') ||
      (opt == 'T') ||
      (opt == 's') ||
      (opt == 't') ||
      (opt == 'w') ||
      (opt == 'o'), "invalid command argument")
    switch (opt) {
      default:
      case 'S':
      case 'T':
      {
        struct grid_config * grid_config =
          (opt == 'S')?src_grid_config:tgt_grid_config;

        YAC_ASSERT_F(
          strlen(optarg) == 1, "invalid grid type for argument %c", (char)opt);

        switch (optarg[0]) {
          case 'c':
          case 'C':
            grid_config->type = CUBE;
            break;
          case 'm':
          case 'M':
            grid_config->type = CURVE;
            break;
          case 'i':
          case 'I':
            grid_config->type = UNSTRUCT;
            break;
          case 'g':
          case 'G':
            grid_config->type = GAUSS;
            break;
          case 's':
          case 'S':
              grid_config->type = SCRIP;
        };
        grid_config->address_offset = (optarg[0] < 'a')?-2:-1;
        break;
      }
      case 's':
        src_arg = optarg;
        break;
      case 't':
        tgt_arg = optarg;
        break;
      case 'w':
        *weight_filename = optarg;
        break;
      case 'o':
        *vtk_filename = optarg;
        break;
    }
  }
  YAC_ASSERT_F(optind >= argc, "non-option ARGV-element: \"%s\"", argv[optind])
  YAC_ASSERT(argc != 1, "too few arguments")
  YAC_ASSERT(*weight_filename != NULL,  "weight_filename argument is missing")
  YAC_ASSERT(*vtk_filename != NULL, "vtk_filename argument is missing")

  interpret_grid_arg(src_grid_config, src_arg, "source");
  interpret_grid_arg(tgt_grid_config, tgt_arg, "target");
}

static double * generate_vertices(double start, double end, size_t count) {

  double * vertices = xmalloc(count * sizeof(*vertices));

  double d = (end - start) / (double)(count - 1);

  for (size_t i = 0; i < count; ++i)
    vertices[i] = start + d * (double)i;
  vertices[count - 1] = end;

  return vertices;
}

static struct basic_grid_data generate_gauss_grid(
  double * first_corner, double * last_corner, size_t * num_cells) {

  size_t num_vertices[2] = {num_cells[0] + 1, num_cells[1] + 1};
  int cyclic[2] = {0, 0};
  double * lon_vertices =
    generate_vertices(first_corner[0], last_corner[0], num_vertices[0]);
  double * lat_vertices =
    generate_vertices(first_corner[1], last_corner[1], num_vertices[1]);

  struct basic_grid_data grid_data =
    yac_generate_basic_grid_data_reg_2d_deg(
      num_vertices, cyclic, lon_vertices, lat_vertices);

  free(lon_vertices);
  free(lat_vertices);

  return grid_data;
}

static struct basic_grid_data create_grid(struct grid_config grid_config) {

  YAC_ASSERT(
    (grid_config.type == CUBE) ||
    (grid_config.type == CURVE) ||
    (grid_config.type == UNSTRUCT) ||
    (grid_config.type == GAUSS) ||
    (grid_config.type == SCRIP), "invalid grid configuration")
  switch(grid_config.type) {
    default:
    case CUBE:
      return generate_cubed_sphere_grid(grid_config.config.cube.n);
    case CURVE:
      YAC_ASSERT_F(
        yac_file_exists(grid_config.config.curve.filename),
        "File %s does not exist.", grid_config.config.curve.filename);
      return read_mpiom_grid((char*)(grid_config.config.curve.filename));
    case UNSTRUCT:
      YAC_ASSERT_F(
        yac_file_exists(grid_config.config.unstruct.filename),
        "File %s does not exist.", grid_config.config.unstruct.filename);
      return read_icon_grid((char*)(grid_config.config.unstruct.filename));
    case GAUSS:
      return
        generate_gauss_grid(
          grid_config.config.gauss.corners[0],
          grid_config.config.gauss.corners[1],
          grid_config.config.gauss.num_cells);
    case SCRIP:
      YAC_ASSERT_F(
        yac_file_exists(grid_config.config.scrip.grid_filename),
        "File %s does not exist.", grid_config.config.scrip.grid_filename);
      YAC_ASSERT_F(
        yac_file_exists(grid_config.config.scrip.mask_filename),
        "File %s does not exist.", grid_config.config.scrip.mask_filename);
      return
        read_scrip_grid(
          (char*)(grid_config.config.scrip.grid_filename),
          (char*)(grid_config.config.scrip.mask_filename),
          (char*)(grid_config.config.scrip.grid_name), 0);
  }
}

