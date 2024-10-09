// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <string.h>

#include "basic_grid.h"
#include "grid_cell.h"
#include "field_data_set.h"
#include "io_utils.h"
#include "yac_mpi_internal.h"
#include "time.h"
#include "yac_config.h"
#include "geometry.h"

struct yac_basic_grid {
  char * name;
  int is_empty;
  struct yac_field_data_set * field_data_set;
  struct yac_basic_grid_data data;
};

static struct yac_basic_grid_data yac_basic_grid_data_empty = {
  .vertex_coordinates = NULL,
  .cell_ids = NULL,
  .vertex_ids = NULL,
  .edge_ids = NULL,
  .num_cells = 0,
  .num_vertices = 0,
  .num_edges = 0,
  .core_cell_mask = NULL,
  .core_vertex_mask = NULL,
  .core_edge_mask = NULL,
  .num_vertices_per_cell = NULL,
  .num_cells_per_vertex = NULL,
  .cell_to_vertex = NULL,
  .cell_to_vertex_offsets = NULL,
  .cell_to_edge = NULL,
  .cell_to_edge_offsets = NULL,
  .vertex_to_cell = NULL,
  .vertex_to_cell_offsets = NULL,
  .edge_to_vertex = NULL,
  .edge_type = NULL,
  .num_total_cells = 0,
  .num_total_vertices = 0,
  .num_total_edges = 0
};

struct yac_basic_grid * yac_basic_grid_new(
  char const * name, struct yac_basic_grid_data grid_data) {

  struct yac_basic_grid * grid = xmalloc(1 * sizeof(*grid));

  grid->name = strdup(name);
  grid->is_empty = 0;
  grid->field_data_set = yac_field_data_set_empty_new();
  grid->data = grid_data;

  return grid;
}

struct yac_basic_grid * yac_basic_grid_empty_new(char const * name) {
  struct yac_basic_grid * grid =
    yac_basic_grid_new(name, yac_basic_grid_data_empty);
  grid->is_empty = 1;
  return grid;
}

void yac_basic_grid_delete(struct yac_basic_grid * grid) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_delete): "
    "NULL is not a valid value for argument grid")
  free(grid->name);
  yac_field_data_set_delete(grid->field_data_set);
  yac_basic_grid_data_free(grid->data);
  free(grid);
}

yac_const_coordinate_pointer yac_basic_grid_get_field_coordinates(
  struct yac_basic_grid * grid, struct yac_interp_field field) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_get_field_coordinates): "
    "NULL is not a valid value for argument grid")

  yac_const_coordinate_pointer coords =
    (field.coordinates_idx != SIZE_MAX)?
      yac_field_data_get_coordinates_data(
        yac_basic_grid_get_field_data(grid, field.location),
        field.coordinates_idx):NULL;

  // if no field coordinates are defined, but the location is at the corners of
  // of the grid cells, return coordinates of them
  return
    ((coords != NULL) || (field.location != YAC_LOC_CORNER))?
       coords:((yac_const_coordinate_pointer)(grid->data.vertex_coordinates));
}

int const * yac_basic_grid_get_core_mask(
  struct yac_basic_grid * grid, enum yac_location location) {

  YAC_ASSERT(
    (location == YAC_LOC_CELL) ||
    (location == YAC_LOC_CORNER) ||
    (location == YAC_LOC_EDGE),
    "ERROR(yac_basic_grid_get_core_mask): invalid location")

  switch (location) {
    default:
    case(YAC_LOC_CELL): return grid->data.core_cell_mask;
    case(YAC_LOC_CORNER): return grid->data.core_vertex_mask;
    case(YAC_LOC_EDGE): return grid->data.core_edge_mask;
  };
}

int const * yac_basic_grid_get_field_mask(
  struct yac_basic_grid * grid, struct yac_interp_field field) {

  if (field.masks_idx == SIZE_MAX) return NULL;

  return
    yac_field_data_get_mask_data(
      yac_basic_grid_get_field_data(grid, field.location), field.masks_idx);
}

char const * yac_basic_grid_get_name(struct yac_basic_grid * grid) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_get_name): "
    "NULL is not a valid value for argument grid")

  return grid->name;
}

struct yac_basic_grid_data * yac_basic_grid_get_data(
  struct yac_basic_grid * grid) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_get_data): "
    "NULL is not a valid value for argument grid")

  return &(grid->data);
}

size_t yac_basic_grid_get_data_size(
  struct yac_basic_grid * grid, enum yac_location location) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_get_data_size): "
    "NULL is not a valid value for argument grid")
  YAC_ASSERT(
    (location == YAC_LOC_CELL) ||
    (location == YAC_LOC_CORNER) ||
    (location == YAC_LOC_EDGE),
    "ERROR(yac_basic_grid_get_data_size): invalid location")

  switch (location) {
    default:
    case (YAC_LOC_CELL):
      return grid->data.num_cells;
    case (YAC_LOC_CORNER):
      return grid->data.num_vertices;
    case (YAC_LOC_EDGE):
      return grid->data.num_edges;
  };
}

size_t yac_basic_grid_get_data_size_f2c(
  struct yac_basic_grid * grid, int location) {

  return
    yac_basic_grid_get_data_size(grid, yac_get_location(location));
}

size_t yac_basic_grid_get_named_mask_idx(
  struct yac_basic_grid * grid, enum yac_location location,
  char const * mask_name) {

  if (mask_name == NULL) return SIZE_MAX;

  struct yac_field_data * data =
    yac_basic_grid_get_field_data(grid, location);

  if (data == NULL) return SIZE_MAX;

  size_t mask_idx = SIZE_MAX;
  size_t masks_count = yac_field_data_get_masks_count(data);

  for (size_t i = 0; (i < masks_count) && (mask_idx == SIZE_MAX); ++i) {
    char const * curr_mask_name =
      yac_field_data_get_mask_name(data, i);
    if ((curr_mask_name != NULL) &&
        (!strcmp(mask_name, curr_mask_name)))
      mask_idx = i;
  }

  YAC_ASSERT_F(
    mask_idx != SIZE_MAX,
    "ERROR(yac_basic_grid_get_named_mask_idx): grid \"%s\" does not contain "
    "%s-mask with the name \"%s\"",
    grid->name, yac_loc2str(location), mask_name)

  return mask_idx;
}

size_t yac_basic_grid_add_coordinates_nocpy(
  struct yac_basic_grid * grid,
  enum yac_location location, yac_coordinate_pointer coordinates) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_add_coordinates_nocpy): "
    "NULL is not a valid value for argument grid")
  YAC_ASSERT_F(
    !grid->is_empty, "ERROR(yac_basic_grid_add_coordinates_nocpy): "
    "grid \"%s\" is an empty grid", grid->name)

  return
    yac_field_data_set_add_coordinates_nocpy(
      grid->field_data_set, location, coordinates);
}

size_t yac_basic_grid_add_coordinates_nocpy_f2c(
  struct yac_basic_grid * grid, int location, double * coordinates) {

  return
    yac_basic_grid_add_coordinates_nocpy(
      grid, yac_get_location(location), (yac_coordinate_pointer)coordinates);
}

size_t yac_basic_grid_add_coordinates(
  struct yac_basic_grid * grid, enum yac_location location,
  yac_coordinate_pointer coordinates, size_t count) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_add_coordinates): "
    "NULL is not a valid value for argument grid")
  YAC_ASSERT_F(
    !grid->is_empty, "ERROR(yac_basic_grid_add_coordinates): "
    "grid \"%s\" is an empty grid", grid->name)

  return
    yac_field_data_set_add_coordinates(
      grid->field_data_set, location, coordinates, count);
}

size_t yac_basic_grid_add_coordinates_f2c(
  struct yac_basic_grid * grid, int location,
  double * coordinates, size_t count) {

  return
    yac_basic_grid_add_coordinates(
      grid, yac_get_location(location),
      (yac_coordinate_pointer)coordinates, count);
}

size_t yac_basic_grid_add_mask_nocpy(
  struct yac_basic_grid * grid,
  enum yac_location location, int const * mask,
  char const * mask_name) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_add_mask_nocpy): "
    "NULL is not a valid value for argument grid")
  YAC_ASSERT_F(
    !grid->is_empty, "ERROR(yac_basic_grid_add_mask_nocpy): "
    "grid \"%s\" is an empty grid", grid->name)

  return
    yac_field_data_set_add_mask_nocpy(
      grid->field_data_set, location, mask, mask_name);
}

size_t yac_basic_grid_add_mask_nocpy_f2c(
  struct yac_basic_grid * grid, int location,
  int const * mask, char const * mask_name) {

  return
    yac_basic_grid_add_mask_nocpy(
      grid, yac_get_location(location), mask, mask_name);
}

size_t yac_basic_grid_add_mask(
  struct yac_basic_grid * grid, enum yac_location location,
  int const * mask, size_t count, char const * mask_name) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_add_mask): "
    "NULL is not a valid value for argument grid")
  YAC_ASSERT_F(
    !grid->is_empty, "ERROR(yac_basic_grid_add_mask): "
    "grid \"%s\" is an empty grid", grid->name)

  return
    yac_field_data_set_add_mask(
      grid->field_data_set, location, mask, count, mask_name);
}

size_t yac_basic_grid_add_mask_f2c(
  struct yac_basic_grid * grid, int location,
  int const * mask, size_t count, char const * mask_name) {

  return
    yac_basic_grid_add_mask(
      grid, yac_get_location(location), mask, count, mask_name);
}

struct yac_field_data * yac_basic_grid_get_field_data(
  struct yac_basic_grid * grid, enum yac_location location) {

  YAC_ASSERT(
    grid, "ERROR(yac_basic_grid_get_field_data): "
    "NULL is not a valid value for argument grid")

  if (grid->is_empty)
    return NULL;
  else
    return
      yac_field_data_set_get_field_data(
        grid->field_data_set, location);
}

struct yac_basic_grid * yac_basic_grid_reg_2d_new(
  char const * name, size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_basic_grid_data_reg_2d(
        nbr_vertices, cyclic, lon_vertices, lat_vertices));
}

struct yac_basic_grid * yac_basic_grid_reg_2d_deg_new(
  char const * name, size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_basic_grid_data_reg_2d_deg(
        nbr_vertices, cyclic, lon_vertices, lat_vertices));
}

struct yac_basic_grid * yac_basic_grid_curve_2d_new(
  char const * name, size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_basic_grid_data_curve_2d(
        nbr_vertices, cyclic, lon_vertices, lat_vertices));
}

struct yac_basic_grid * yac_basic_grid_curve_2d_deg_new(
  char const * name, size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_basic_grid_data_curve_2d_deg(
        nbr_vertices, cyclic, lon_vertices, lat_vertices));
}

struct yac_basic_grid * yac_basic_grid_unstruct_new(
  char const * name, size_t nbr_vertices, size_t nbr_cells,
  int *num_vertices_per_cell, double *x_vertices, double *y_vertices,
  int *cell_to_vertex) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_basic_grid_data_unstruct(
        nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex));
}

struct yac_basic_grid * yac_basic_grid_unstruct_deg_new(
  char const * name, size_t nbr_vertices, size_t nbr_cells,
  int *num_vertices_per_cell, double *x_vertices, double *y_vertices,
  int *cell_to_vertex) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_basic_grid_data_unstruct_deg(
        nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex));
}

struct yac_basic_grid * yac_basic_grid_unstruct_ll_new(
  char const * name, size_t nbr_vertices, size_t nbr_cells,
  int *num_vertices_per_cell, double *x_vertices, double *y_vertices,
  int *cell_to_vertex) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_basic_grid_data_unstruct_ll(
        nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex));
}

struct yac_basic_grid * yac_basic_grid_unstruct_ll_deg_new(
  char const * name, size_t nbr_vertices, size_t nbr_cells,
  int *num_vertices_per_cell, double *x_vertices, double *y_vertices,
  int *cell_to_vertex) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_basic_grid_data_unstruct_ll_deg(
        nbr_vertices, nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex));
}

struct yac_basic_grid * yac_basic_grid_cloud_new(
  char const * name, size_t nbr_points, double *x_points, double *y_points) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_basic_grid_data_cloud(nbr_points, x_points, y_points));
}

struct yac_basic_grid * yac_basic_grid_cloud_deg_new(
  char const * name, size_t nbr_points, double *x_points, double *y_points) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_basic_grid_data_cloud_deg(nbr_points, x_points, y_points));
}

#ifdef YAC_NETCDF_ENABLED

static int def_dim(
  char const * filename, int ncid, char const * name, size_t dim_len,
  int file_is_new) {

  int create_dim = file_is_new;
  int dim_id;

  // if the netcdf file already exists
  if (!file_is_new) {

    // check whether the dimension already exists in the file
    int status = nc_inq_dimid(ncid, name, &dim_id);

    // if the dimension already exists, check for a consistent dimension size
    if (!((create_dim = (status == NC_EBADDIM)))) {

      YAC_HANDLE_ERROR(status);

      size_t temp_dim_len;
      YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &temp_dim_len));

      YAC_ASSERT_F(
        dim_len == temp_dim_len,
        "ERROR(def_dim): file \"%s\" already contains dimension \"%s\", "
        "but it has a different size %zu != %zu",
        filename, name, dim_len, temp_dim_len);
    }
  }

  if (create_dim) YAC_HANDLE_ERROR(nc_def_dim(ncid, name, dim_len, &dim_id));

  return dim_id;
}

static void def_var(
  char const * filename, int ncid, char const * name, nc_type type,
  int ndims, const int * dimids, char const * att_name, char const * att_value,
  int file_is_new) {

  int create_var = file_is_new;
  int var_id;

  // if the netcdf file already exists
  if (!file_is_new) {

    // check whether the variable already exists in the file
    int status = nc_inq_varid(ncid, name, &var_id);

    // if the variable already exists
    if (!((create_var = (status == NC_ENOTVAR)))) {

      YAC_HANDLE_ERROR(status);

      // check the consistency of the variable definition
      nc_type temp_type;
      int temp_ndims;
      int temp_dimids[NC_MAX_VAR_DIMS];
      YAC_HANDLE_ERROR(
        nc_inq_var(
          ncid, var_id, NULL, &temp_type, &temp_ndims, temp_dimids, NULL));
      YAC_ASSERT_F(
        type == temp_type,
        "ERROR(def_var): file \"%s\" already contains variable \"%s\", "
        "but it has a different data type %d != %d",
        filename, name, type, temp_type);
      YAC_ASSERT_F(
        ndims == temp_ndims,
        "ERROR(def_var): file \"%s\" already contains variable \"%s\", "
        "but it has a different number of dimensions %d != %d",
        filename, name, ndims, temp_ndims);
      for (int i = 0; i < ndims; ++i) {
        YAC_ASSERT_F(
          dimids[i] == temp_dimids[i],
          "ERROR(def_var): file \"%s\" already contains variable \"%s\", "
          "but it has a different dimensions dimid[%d] %d != %d",
          filename, name, i, dimids[i], temp_dimids[i])
      }

      if (att_name != NULL) {

        size_t att_len;
        YAC_HANDLE_ERROR(nc_inq_attlen(ncid, var_id, att_name, &att_len));

        char * temp_att_value = xcalloc((att_len + 1), sizeof(temp_att_value));
        YAC_HANDLE_ERROR(nc_get_att_text(ncid, var_id, att_name, temp_att_value));

        YAC_ASSERT_F(
          !strcmp(att_value, temp_att_value),
          "ERROR(def_var): file \"%s\" already contains variable \"%s\", "
          "but it has a different attribute value (name: \"%s\") "
          "\"%s\" != \"%s\"",
          filename, name, att_name, att_value, temp_att_value);
        free(temp_att_value);
      }
    }
  }

  if (create_var) {
    YAC_HANDLE_ERROR(nc_def_var(ncid, name, type, ndims, dimids, &var_id));
    if (att_name != NULL)
      YAC_HANDLE_ERROR(
        nc_put_att_text(ncid, var_id, att_name, strlen(att_value), att_value));
  }
}

static void create_grid_file(
  char const * filename, char const * grid_name, size_t num_cells,
  int num_vertices_per_cell, int cell_center_coords_avaiable,
  int cell_global_ids_available, int core_cell_mask_available,
  int vertex_global_ids_available, int core_vertex_mask_available,
  int edge_global_ids_available, int core_edge_mask_available) {

  size_t grid_name_len = strlen(grid_name) + 1;

  char nv_dim_name[3 + grid_name_len];
  char nc_dim_name[3 + grid_name_len];

  char cla_var_name[4 + grid_name_len];
  char clo_var_name[4 + grid_name_len];
  char lat_var_name[4 + grid_name_len];
  char lon_var_name[4 + grid_name_len];
  char rnk_var_name[4 + grid_name_len];
  char cmk_var_name[4 + grid_name_len];
  char gid_var_name[4 + grid_name_len];
  char vcmk_var_name[5 + grid_name_len];
  char vgid_var_name[5 + grid_name_len];
  char ecmk_var_name[5 + grid_name_len];
  char egid_var_name[5 + grid_name_len];
  char const * unit_att = "units";
  char const * coord_unit = "degree";

  snprintf(nv_dim_name, 3 + grid_name_len, "nv_%s", grid_name);
  snprintf(nc_dim_name, 3 + grid_name_len, "nc_%s", grid_name);

  snprintf(cla_var_name, 4 + grid_name_len, "%s.cla", grid_name);
  snprintf(clo_var_name, 4 + grid_name_len, "%s.clo", grid_name);
  snprintf(lat_var_name, 4 + grid_name_len, "%s.lat", grid_name);
  snprintf(lon_var_name, 4 + grid_name_len, "%s.lon", grid_name);
  snprintf(rnk_var_name, 4 + grid_name_len, "%s.rnk", grid_name);
  snprintf(cmk_var_name, 4 + grid_name_len, "%s.cmk", grid_name);
  snprintf(gid_var_name, 4 + grid_name_len, "%s.gid", grid_name);
  snprintf(vcmk_var_name, 5 + grid_name_len, "%s.vcmk", grid_name);
  snprintf(vgid_var_name, 5 + grid_name_len, "%s.vgid", grid_name);
  snprintf(ecmk_var_name, 5 + grid_name_len, "%s.ecmk", grid_name);
  snprintf(egid_var_name, 5 + grid_name_len, "%s.egid", grid_name);

  int ncid;

  // create/open file
  int file_is_new = !yac_file_exists(filename);
  if (file_is_new) {
    yac_nc_create(filename, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
  } else {
    yac_nc_open(filename, NC_WRITE | NC_SHARE, &ncid);
    nc_redef(ncid);
  }

  // define dimensions
  int nv_dim_id =
    def_dim(filename, ncid, nv_dim_name, (size_t)num_vertices_per_cell, file_is_new);
  int nc_dim_id =
    def_dim(filename, ncid, nc_dim_name, (size_t)num_cells, file_is_new);

  // define variables
  int corner_dims[2] = {nc_dim_id, nv_dim_id};
  int cell_dims[1] = {nc_dim_id};
  def_var(
    filename, ncid, cla_var_name, NC_DOUBLE, 2, corner_dims,
    unit_att, coord_unit, file_is_new);
  def_var(
    filename, ncid, clo_var_name, NC_DOUBLE, 2, corner_dims,
    unit_att, coord_unit, file_is_new);
  if (cell_center_coords_avaiable) {
    def_var(
      filename, ncid, lat_var_name, NC_DOUBLE, 1, cell_dims,
      unit_att, coord_unit, file_is_new);
    def_var(
      filename, ncid, lon_var_name, NC_DOUBLE, 1, cell_dims,
      unit_att, coord_unit, file_is_new);
  }
  if (cell_global_ids_available)
    def_var(
      filename, ncid, gid_var_name, NC_INT, 1, cell_dims, NULL, NULL, file_is_new);
  if (core_cell_mask_available)
    def_var(
      filename, ncid, cmk_var_name, NC_INT, 1, cell_dims, NULL, NULL, file_is_new);
  if (vertex_global_ids_available)
    def_var(
      filename, ncid, vgid_var_name, NC_INT, 2, corner_dims, NULL, NULL, file_is_new);
  if (core_vertex_mask_available)
    def_var(
      filename, ncid, vcmk_var_name, NC_INT, 2, corner_dims, NULL, NULL, file_is_new);
  if (edge_global_ids_available)
    def_var(
      filename, ncid, egid_var_name, NC_INT, 2, corner_dims, NULL, NULL, file_is_new);
  if (core_edge_mask_available)
    def_var(
      filename, ncid, ecmk_var_name, NC_INT, 2, corner_dims, NULL, NULL, file_is_new);
  def_var(
    filename, ncid, rnk_var_name, NC_INT, 1, cell_dims, NULL, NULL, file_is_new);

  time_t now = time(NULL);
  char str_now_UTC[32];
  strftime(
    str_now_UTC, sizeof(str_now_UTC), "%Y-%m-%dT%H:%M:%SZ",
    gmtime(&now));
  char const * yac_version = YAC_VERSION;
  char const * created_by = "Created by YAC";
  char const * grid_type = "curvilinear";
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, NC_GLOBAL, "YAC", strlen(yac_version), yac_version));
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, NC_GLOBAL, "title", strlen(created_by), created_by));
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, NC_GLOBAL, "description", strlen(created_by), created_by));
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, NC_GLOBAL, "grid", strlen(grid_type), grid_type));
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, NC_GLOBAL, "timeStamp", strlen(str_now_UTC), str_now_UTC));

  // end definition
  YAC_HANDLE_ERROR(nc_enddef(ncid));

  // close file
  YAC_HANDLE_ERROR(nc_close(ncid));
}

static void put_vara(
  int ncid, char const * grid_name, char const * var_ext,
  size_t * start, size_t * count, void const * buffer) {

  if (count[0] == 0) return;

  char * var_name =
    xmalloc((strlen(grid_name) + strlen(var_ext) + 2) * sizeof(*var_name));

  sprintf(var_name, "%s.%s", grid_name, var_ext);

  int var_id;

  yac_nc_inq_varid(ncid, var_name, &var_id);

  free(var_name);

  YAC_HANDLE_ERROR(
    nc_put_vara(ncid, var_id, start, count, buffer));
}

#endif // YAC_NETCDF_ENABLED

void yac_basic_grid_to_file_parallel(
  struct yac_basic_grid * grid, char const * filename, MPI_Comm comm) {

#ifndef YAC_NETCDF_ENABLED

  UNUSED(grid);
  UNUSED(filename);
  UNUSED(comm);

  die(
    "ERROR(yac_basic_grid_to_file_parallel): "
    "YAC is built without the NetCDF support");
#else

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  { // consistency check
    int filename_len = (int)strlen(filename) + 1;
    yac_mpi_call(MPI_Bcast(&filename_len, 1, MPI_INT, 0, comm), comm);
    char * filename_buffer =
      xmalloc((size_t)filename_len * sizeof(*filename_buffer));
    if (comm_rank == 0) strcpy(filename_buffer, filename);
    yac_mpi_call(
      MPI_Bcast(filename_buffer, filename_len, MPI_CHAR, 0, comm), comm);

    YAC_ASSERT_F(
      !strcmp(filename, filename_buffer),
      "ERROR(yac_basic_grid_to_file_parallel): "
      "inconsistent filename (\"%s\" on rank 0 != \"%s\" on rank %d)",
      filename_buffer, filename, comm_rank);
    free(filename_buffer);
  }

  struct yac_field_data * cell_field =
    yac_basic_grid_get_field_data(grid, YAC_LOC_CELL);
  yac_const_coordinate_pointer cell_center_coords =
    ((cell_field != NULL) &&
     (yac_field_data_get_coordinates_count(cell_field) > 0))?
      yac_field_data_get_coordinates_data(cell_field, 0):NULL;
  struct yac_basic_grid_data * grid_data = yac_basic_grid_get_data(grid);

  uint64_t local_cell_count = (uint64_t)grid->data.num_cells;
  uint64_t * global_cell_counts =
    xmalloc((size_t)comm_size * sizeof(*global_cell_counts));
  int max_num_vertices_per_cell = 0;
  int cell_center_coords_avaiable =
    (grid_data->num_cells == 0) || (cell_center_coords != NULL);
  int cell_global_ids_available =
    (grid_data->num_cells == 0) || (grid_data->cell_ids != NULL);
  int core_cell_mask_available =
    (grid_data->num_cells == 0) || (grid_data->core_cell_mask != NULL);
  int vertex_global_ids_available =
    (grid_data->num_vertices == 0) || (grid_data->vertex_ids != NULL);
  int core_vertex_mask_available =
    (grid_data->num_vertices == 0) || (grid_data->core_vertex_mask != NULL);
  int edge_global_ids_available =
    (grid_data->num_edges == 0) || (grid_data->edge_ids != NULL);
  int core_edge_mask_available =
    (grid_data->num_edges == 0) || (grid_data->core_edge_mask != NULL);

  for (size_t i = 0; i < grid->data.num_cells; ++i)
    if (grid->data.num_vertices_per_cell[i] > max_num_vertices_per_cell)
      max_num_vertices_per_cell = grid->data.num_vertices_per_cell[i];

  yac_mpi_call(
    MPI_Allgather(
      &local_cell_count, 1, MPI_UINT64_T,
      global_cell_counts, 1, MPI_UINT64_T, comm), comm);
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &max_num_vertices_per_cell, 1, MPI_INT, MPI_MAX, comm),
    comm);
  int ints[7] = {cell_center_coords_avaiable,
                 cell_global_ids_available,
                 core_cell_mask_available,
                 vertex_global_ids_available,
                 core_vertex_mask_available,
                 edge_global_ids_available,
                 core_edge_mask_available};
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, ints, sizeof(ints)/sizeof(ints[0]),
      MPI_INT, MPI_MIN, comm), comm);
  cell_center_coords_avaiable = ints[0];
  cell_global_ids_available = ints[1];
  core_cell_mask_available = ints[2];
  vertex_global_ids_available = ints[3];
  core_vertex_mask_available = ints[4];
  edge_global_ids_available = ints[5];
  core_edge_mask_available = ints[6];

  size_t global_cell_count = 0;
  for (int i = 0; i < comm_size; ++i)
    global_cell_count += (size_t)(global_cell_counts[i]);

  YAC_ASSERT_F(
    global_cell_count > 0,
    "ERROR(yac_basic_grid_to_file_parallel): grid \"%s\" has no cells",
    grid->name);

  // create the grid file
  if (comm_rank == 0)
    create_grid_file(
      filename, grid->name, global_cell_count, max_num_vertices_per_cell,
      cell_center_coords_avaiable,
      cell_global_ids_available, core_cell_mask_available,
      vertex_global_ids_available, core_vertex_mask_available,
      edge_global_ids_available, core_edge_mask_available);
  yac_mpi_call(MPI_Barrier(comm), comm);

  // determine processes that will do output
  int io_flag;
  int * io_ranks;
  int num_io_ranks;
  yac_get_io_ranks(comm, &io_flag, &io_ranks, &num_io_ranks);

  size_t recv_count = 0;
  size_t io_start_idx = SIZE_MAX;

  if (io_flag) {

    int io_rank_idx;
    for (io_rank_idx = 0; io_rank_idx < num_io_ranks; ++io_rank_idx)
      if (io_ranks[io_rank_idx] == comm_rank) break;

    io_start_idx =
      (size_t)(
        ((long long)(io_rank_idx) * (long long)global_cell_count)/
        (long long)num_io_ranks);
    recv_count =
      (size_t)(
        ((long long)(io_rank_idx + 1) * (long long)global_cell_count)/
        (long long)num_io_ranks - (long long)io_start_idx);
  }

  // generate transfer positions
  struct Xt_com_pos * com_pos_buffer =
    xmalloc(5 * (size_t)comm_size * sizeof(*com_pos_buffer));
  struct Xt_com_pos * cell_src_com = com_pos_buffer;
  struct Xt_com_pos * cell_dst_com = com_pos_buffer + (size_t)comm_size;
  struct Xt_com_pos * vertex_src_com = com_pos_buffer + 2 * (size_t)comm_size;
  struct Xt_com_pos * vertex_dst_com = com_pos_buffer + 3 * (size_t)comm_size;
  struct Xt_com_pos * edge_src_com = com_pos_buffer + 4 * (size_t)comm_size;
  struct Xt_com_pos * edge_dst_com = vertex_dst_com;
  int num_src_msg = 0, num_dst_msg = 0;
  int * transfer_pos_buffer =
    xmalloc(
      (size_t)(2 * max_num_vertices_per_cell + 1) *
      (grid->data.num_cells + recv_count) * sizeof(*transfer_pos_buffer));
  int * cell_send_pos = transfer_pos_buffer;
  int * cell_recv_pos = transfer_pos_buffer + grid->data.num_cells;
  int * vertex_send_pos =
    transfer_pos_buffer + grid->data.num_cells + recv_count;
  int * vertex_recv_pos =
    transfer_pos_buffer + grid->data.num_cells + recv_count +
    (size_t)max_num_vertices_per_cell * grid->data.num_cells;
  int * edge_send_pos =
    transfer_pos_buffer + grid->data.num_cells + recv_count +
    (size_t)max_num_vertices_per_cell * grid->data.num_cells +
    (size_t)max_num_vertices_per_cell * recv_count;

  for (size_t i = 0; i < grid->data.num_cells; ++i) {

    size_t * cell_to_vertex =
      grid->data.cell_to_vertex + grid->data.cell_to_vertex_offsets[i];
    size_t * cell_to_edge =
      grid->data.cell_to_edge + grid->data.cell_to_edge_offsets[i];
    int num_vertices = grid->data.num_vertices_per_cell[i];

    cell_send_pos[i] = (int)i;
    for (int j = 0; j < num_vertices; ++j)
      vertex_send_pos[i * (size_t)max_num_vertices_per_cell + (size_t)j] =
        cell_to_vertex[j];
    for (int j = num_vertices; j < max_num_vertices_per_cell; ++j)
      vertex_send_pos[i * (size_t)max_num_vertices_per_cell + (size_t)j] = 0;
    for (int j = 0; j < num_vertices; ++j)
      edge_send_pos[i * (size_t)max_num_vertices_per_cell + (size_t)j] =
        cell_to_edge[j];
    for (int j = num_vertices; j < max_num_vertices_per_cell; ++j)
      edge_send_pos[i * (size_t)max_num_vertices_per_cell + (size_t)j] = 0;
  }

  for (size_t i = 0, k = 0; i < recv_count; ++i) {
    cell_recv_pos[i] = (int)i;
    for (int j = 0; j < max_num_vertices_per_cell; ++j, ++k)
      vertex_recv_pos[k] = k;
  }

  size_t io_end_idx = 0;
  size_t global_count_accu = 0;
  for (int io_rank_idx = 0, curr_rank = 0; io_rank_idx < num_io_ranks;
       ++io_rank_idx) {

    int curr_io_rank = io_ranks[io_rank_idx];
    size_t io_start_idx = io_end_idx;
    io_end_idx =
      (size_t)(
        ((long long)(io_rank_idx + 1) * (long long)global_cell_count)/
        (long long)num_io_ranks);
    size_t io_size = io_end_idx - io_start_idx;
    size_t curr_recv_size = 0;

    while(global_count_accu < io_end_idx) {

      size_t count =
        MIN(
          (global_count_accu + global_cell_counts[curr_rank]) -
          (io_start_idx + curr_recv_size), io_size - curr_recv_size);

      if (curr_rank == comm_rank) {
        cell_src_com[num_src_msg].transfer_pos = cell_send_pos;
        cell_src_com[num_src_msg].num_transfer_pos = (int)count;
        cell_src_com[num_src_msg].rank = curr_io_rank;
        cell_send_pos += count;

        vertex_src_com[num_src_msg].transfer_pos = vertex_send_pos;
        vertex_src_com[num_src_msg].num_transfer_pos =
          (int)count * max_num_vertices_per_cell;
        vertex_src_com[num_src_msg].rank = curr_io_rank;
        vertex_send_pos +=
          count * (size_t)max_num_vertices_per_cell;

        edge_src_com[num_src_msg].transfer_pos = edge_send_pos;
        edge_src_com[num_src_msg].num_transfer_pos =
          (int)count * max_num_vertices_per_cell;
        edge_src_com[num_src_msg].rank = curr_io_rank;
        edge_send_pos +=
          count * (size_t)max_num_vertices_per_cell;

        num_src_msg++;
      }

      if (curr_io_rank == comm_rank) {
        cell_dst_com[num_dst_msg].transfer_pos = cell_recv_pos;
        cell_dst_com[num_dst_msg].num_transfer_pos = (int)count;
        cell_dst_com[num_dst_msg].rank = curr_rank;
        cell_recv_pos += count;

        vertex_dst_com[num_dst_msg].transfer_pos = vertex_recv_pos;
        vertex_dst_com[num_dst_msg].num_transfer_pos =
          (int)count * max_num_vertices_per_cell;
        vertex_dst_com[num_dst_msg].rank = curr_rank;
        vertex_recv_pos += count * (size_t)max_num_vertices_per_cell;

        num_dst_msg++;
      }

      if ((global_count_accu + global_cell_counts[curr_rank]) <=
          io_end_idx) {
        global_count_accu += global_cell_counts[curr_rank];
        curr_rank++;
      }

      curr_recv_size += count;

      if (curr_recv_size >= io_size) break;
    }
  }

  free(global_cell_counts);

  // open grid file
  int ncid;
  size_t start[2], count[2];
  if (io_flag) {

    int io_rank_idx;
    for (io_rank_idx = 0; io_rank_idx < num_io_ranks; ++io_rank_idx)
      if (io_ranks[io_rank_idx] == comm_rank) break;

    start[0] = io_start_idx;
    start[1] = 0;
    count[0] = recv_count;
    count[1] = max_num_vertices_per_cell;

    yac_nc_open(filename, NC_WRITE | NC_SHARE, &ncid);
  } else {
    count[0] = 0;
    count[1] = 0;
  }

  void * recv_buffer =
    xmalloc(
      (size_t)max_num_vertices_per_cell * recv_count *
      MAX(MAX(sizeof(double), sizeof(int)), sizeof(yac_int)));

  //----------------------------------------------------------------------------
  // redistribute cell based data
  //----------------------------------------------------------------------------

  Xt_xmap cell_xmap =
    xt_xmap_intersection_pos_new(
      num_src_msg, cell_src_com, num_dst_msg, cell_dst_com, comm);

  Xt_redist cell_redist_int = xt_redist_p2p_new(cell_xmap, MPI_INT);

  // number of vertices per cell
  int * num_vertices_per_cell =
    xmalloc(recv_count * sizeof(*num_vertices_per_cell));
  xt_redist_s_exchange1(
    cell_redist_int, grid->data.num_vertices_per_cell, num_vertices_per_cell);

  // cell center coordinates (if available)
  if (cell_center_coords_avaiable) {

    // generate cell center lon/lat data
    double * send_buffer_dble =
      xmalloc(2 * grid->data.num_cells * sizeof(*send_buffer_dble));;
    double * send_buffer_lon = send_buffer_dble;
    double * send_buffer_lat = send_buffer_dble + grid->data.num_cells;
    for (size_t i = 0; i < grid->data.num_cells; ++i) {
      XYZtoLL(
        cell_center_coords[i], send_buffer_lon + i, send_buffer_lat + i);
      send_buffer_lon[i] /= YAC_RAD;
      send_buffer_lat[i] /= YAC_RAD;
    }

    // exchange cell center lon/lat data and write it to file
    Xt_redist cell_redist_dble = xt_redist_p2p_new(cell_xmap, MPI_DOUBLE);
    xt_redist_s_exchange1(cell_redist_dble, send_buffer_lon, recv_buffer);
    put_vara(ncid, grid->name, "lon", start, count, recv_buffer);
    xt_redist_s_exchange1(cell_redist_dble, send_buffer_lat, recv_buffer);
    put_vara(ncid, grid->name, "lat", start, count, recv_buffer);
    xt_redist_delete(cell_redist_dble);
    free(send_buffer_dble);
  }

  int * cell_send_buffer_int =
    xmalloc(grid->data.num_cells * sizeof(*cell_send_buffer_int));

  // cell ranks
  {
    // generate cell owner rank data
    for (size_t i = 0; i < grid->data.num_cells; ++i)
      cell_send_buffer_int[i] = comm_rank;

    // exchange cell owner ranks and write it to file
    xt_redist_s_exchange1(cell_redist_int, cell_send_buffer_int, recv_buffer);
    put_vara(ncid, grid->name, "rnk", start, count, recv_buffer);
  }

  // cell core mask (if available)
  if (core_cell_mask_available) {

    // exchange core cell mask and write to file
    xt_redist_s_exchange1(
      cell_redist_int, grid->data.core_cell_mask, recv_buffer);
    put_vara(ncid, grid->name, "cmk", start, count, recv_buffer);
  }

  // cell global ids (if available)
  if (cell_global_ids_available) {

    // convert global ids to int
    for (size_t i = 0; i < grid->data.num_cells; ++i)
      cell_send_buffer_int[i] = (int)(grid->data.cell_ids[i]);

    // exchange cell global ids and write to file
    xt_redist_s_exchange1(
      cell_redist_int, cell_send_buffer_int, recv_buffer);
    put_vara(ncid, grid->name, "gid", start, count, recv_buffer);
  }

  xt_redist_delete(cell_redist_int);
  free(cell_send_buffer_int);

  xt_xmap_delete(cell_xmap);

  //----------------------------------------------------------------------------
  // redistribute vertex based data
  //----------------------------------------------------------------------------

  Xt_xmap vertex_xmap =
    xt_xmap_intersection_pos_new(
      num_src_msg, vertex_src_com, num_dst_msg, vertex_dst_com, comm);

  // vertex coordinates
  {
    // generate vertex lon/lat data
    double * send_buffer_dble =
      xmalloc(2 * grid->data.num_vertices * sizeof(*send_buffer_dble));
    double * send_buffer_lon = send_buffer_dble;
    double * send_buffer_lat = send_buffer_dble + grid->data.num_vertices;
    for (size_t i = 0; i < grid->data.num_vertices; ++i) {
      XYZtoLL(
        grid->data.vertex_coordinates[i],
        send_buffer_lon + i, send_buffer_lat + i);
      send_buffer_lon[i] /= YAC_RAD;
      send_buffer_lat[i] /= YAC_RAD;
    }

    // exchange vertex lon/lat data and write it to file
    Xt_redist vertex_redist_dble = xt_redist_p2p_new(vertex_xmap, MPI_DOUBLE);
    xt_redist_s_exchange1(vertex_redist_dble, send_buffer_lon, recv_buffer);
    for (size_t i = 0, k = 0; i < recv_count;
         ++i, k += (size_t)max_num_vertices_per_cell)
      for (size_t j = (size_t)num_vertices_per_cell[i];
           j < (size_t)max_num_vertices_per_cell; ++j)
        ((double*)recv_buffer)[k+j] = DBL_MAX;
    put_vara(ncid, grid->name, "clo", start, count, recv_buffer);
    xt_redist_s_exchange1(vertex_redist_dble, send_buffer_lat, recv_buffer);
    free(send_buffer_dble);
    for (size_t i = 0, k = 0; i < recv_count;
         ++i, k += (size_t)max_num_vertices_per_cell)
      for (size_t j = (size_t)num_vertices_per_cell[i];
           j < (size_t)max_num_vertices_per_cell; ++j)
        ((double*)recv_buffer)[k+j] = DBL_MAX;
    put_vara(ncid, grid->name, "cla", start, count, recv_buffer);
    xt_redist_delete(vertex_redist_dble);
  }

  Xt_redist vertex_redist_int = xt_redist_p2p_new(vertex_xmap, MPI_INT);
  int * vertex_send_buffer_int =
    xmalloc(grid->data.num_vertices * sizeof(*vertex_send_buffer_int));

  // core vertex mask (if available)
  if (core_vertex_mask_available) {

    // exchange core vertex mask and write to file
    xt_redist_s_exchange1(
      vertex_redist_int, grid->data.core_vertex_mask, recv_buffer);
    for (size_t i = 0, k = 0; i < recv_count;
         ++i, k += (size_t)max_num_vertices_per_cell)
      for (size_t j = (size_t)num_vertices_per_cell[i];
           j < (size_t)max_num_vertices_per_cell; ++j)
        ((int*)recv_buffer)[k+j] = INT_MAX;
    put_vara(ncid, grid->name, "vcmk", start, count, recv_buffer);
  }

  // vertex global ids (if available)
  if (vertex_global_ids_available) {

    // convert global ids to int
    for (size_t i = 0; i < grid->data.num_vertices; ++i)
      vertex_send_buffer_int[i] = (int)(grid->data.vertex_ids[i]);

    // exchange vertex global ids and write to file
    xt_redist_s_exchange1(
      vertex_redist_int, vertex_send_buffer_int, recv_buffer);
    for (size_t i = 0, k = 0; i < recv_count;
         ++i, k += (size_t)max_num_vertices_per_cell)
      for (size_t j = (size_t)num_vertices_per_cell[i];
           j < (size_t)max_num_vertices_per_cell; ++j)
        ((int*)recv_buffer)[k+j] = INT_MAX;
    put_vara(ncid, grid->name, "vgid", start, count, recv_buffer);
  }

  free(vertex_send_buffer_int);
  xt_redist_delete(vertex_redist_int);

  xt_xmap_delete(vertex_xmap);

  //----------------------------------------------------------------------------
  // redistribute edge based data
  //----------------------------------------------------------------------------

  if (core_edge_mask_available || edge_global_ids_available) {

    Xt_xmap edge_xmap =
      xt_xmap_intersection_pos_new(
        num_src_msg, edge_src_com, num_dst_msg, edge_dst_com, comm);

    Xt_redist edge_redist_int = xt_redist_p2p_new(edge_xmap, MPI_INT);
    int * edge_send_buffer_int =
      xmalloc(grid->data.num_edges * sizeof(*edge_send_buffer_int));

    // core edge mask (if available)
    if (core_edge_mask_available) {

      // exchange core edge mask and write to file
      xt_redist_s_exchange1(
        edge_redist_int, grid->data.core_edge_mask, recv_buffer);
      for (size_t i = 0, k = 0; i < recv_count;
          ++i, k += (size_t)max_num_vertices_per_cell)
        for (size_t j = (size_t)num_vertices_per_cell[i];
            j < (size_t)max_num_vertices_per_cell; ++j)
          ((int*)recv_buffer)[k+j] = INT_MAX;
      put_vara(ncid, grid->name, "ecmk", start, count, recv_buffer);
    }

    // edge global ids (if available)
    if (edge_global_ids_available) {

      // convert global ids to int
      for (size_t i = 0; i < grid->data.num_edges; ++i)
        edge_send_buffer_int[i] = (int)(grid->data.edge_ids[i]);

      // exchange edge global ids and write to file
      xt_redist_s_exchange1(
        edge_redist_int, edge_send_buffer_int, recv_buffer);
      for (size_t i = 0, k = 0; i < recv_count;
          ++i, k += (size_t)max_num_vertices_per_cell)
        for (size_t j = (size_t)num_vertices_per_cell[i];
            j < (size_t)max_num_vertices_per_cell; ++j)
          ((int*)recv_buffer)[k+j] = INT_MAX;
      put_vara(ncid, grid->name, "egid", start, count, recv_buffer);
    }

    free(edge_send_buffer_int);
    xt_redist_delete(edge_redist_int);

    xt_xmap_delete(edge_xmap);
  }

  free(num_vertices_per_cell);
  free(recv_buffer);
  free(com_pos_buffer);
  free(transfer_pos_buffer);
  free(io_ranks);

  // ensure that the writing of the weight file is complete
  yac_mpi_call(MPI_Barrier(comm), comm);

#endif // YAC_NETCDF_ENABLED
}

void yac_basic_grid_compute_cell_areas(
  struct yac_basic_grid * grid, double * cell_areas) {

  yac_basic_grid_data_compute_cell_areas(
    grid->data, cell_areas);
}
