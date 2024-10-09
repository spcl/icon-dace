// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

// #define VERBOSE

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <string.h>
#include <stdbool.h>
#include "yac.h"
#include "yac_utils.h"

enum interp_stack_type {
  INTERP_STACK_CONSERV_DESTAREA,
  INTERP_STACK_CONSERV_FRACAREA,
  INTERP_STACK_CONS2ND_FRACAREA,
  INTERP_STACK_DISTWGT_4,
  INTERP_STACK_DISTWGT_1,
  INTERP_STACK_HCSBB,
  INTERP_STACK_AVG_DIST,
  NUM_INTERP_STACK_TYPES
};

enum experiment_type {
    GEN_MASK,
    NOICOH,
    NOGT_ICOH,
    ICOS_ICOH,
    NUM_EXPERIMENT,
};

struct experiment_configuration {
  bool write_weights;
  char const **component_names;
  size_t num_components;
  enum interp_stack_type * interpolations;
  size_t num_interpolations;
} experiment_configs[] =
  {// GEN_MASK
   {.write_weights = 0,
    .component_names =
      (char const*[]){"torc","bggd","nogt","sse7","icos","icoh"},
    .num_components = 6,
    .interpolations =
      (enum interp_stack_type[]){INTERP_STACK_CONSERV_DESTAREA},
    .num_interpolations = 1},
   // NOICOH
   {.write_weights = 1,
    .component_names =
      (char const*[]){"torc","bggd","nogt","sse7","icos"},
    .num_components = 5,
    .interpolations =
      (enum interp_stack_type[]){INTERP_STACK_CONSERV_DESTAREA,
                                 INTERP_STACK_CONSERV_FRACAREA,
                                 INTERP_STACK_CONS2ND_FRACAREA,
                                 INTERP_STACK_DISTWGT_4,
                                 INTERP_STACK_DISTWGT_1,
                                 INTERP_STACK_HCSBB,
                                 INTERP_STACK_AVG_DIST},
    .num_interpolations = 7},
   // NOGT_ICOH
   {.write_weights = 1,
    .component_names = (char const*[]){"nogt","icoh"},
    .num_components = 2,
    .interpolations =
      (enum interp_stack_type[]){INTERP_STACK_CONSERV_DESTAREA,
                                 INTERP_STACK_CONSERV_FRACAREA,
                                 INTERP_STACK_CONS2ND_FRACAREA,
                                 INTERP_STACK_DISTWGT_4,
                                 INTERP_STACK_DISTWGT_1,
                                 INTERP_STACK_HCSBB,
                                 INTERP_STACK_AVG_DIST},
    .num_interpolations = 7},
   // ICOS_ICOH
   {.write_weights = 1,
    .component_names = (char const*[]){"icos","icoh"},
    .num_components = 2,
    .interpolations =
      (enum interp_stack_type[]){INTERP_STACK_CONSERV_DESTAREA,
                                 INTERP_STACK_CONSERV_FRACAREA,
                                 INTERP_STACK_CONS2ND_FRACAREA,
                                 INTERP_STACK_DISTWGT_4,
                                 INTERP_STACK_DISTWGT_1,
                                 INTERP_STACK_HCSBB,
                                 INTERP_STACK_AVG_DIST},
    .num_interpolations = 7}};

struct {
  char const * name;
  int id;
} interp_stacks[NUM_INTERP_STACK_TYPES];

struct {
  double (*p)(double lon, double lat);
  char const * name;
} test_functions[] =
  {{.p = yac_test_func,
    .name = "yac"},
    {.p = yac_test_ana_fcos,
    .name = "fcos"},
    {.p = yac_test_ana_fcossin,
    .name = "fcossin"},
    {.p = yac_test_one,
    .name = "one"},
    {.p = yac_test_gulfstream,
    .name = "gulfstream"},
    {.p = yac_test_harmonic,
    .name = "harmonic"},
    {.p = yac_test_vortex,
    .name = "vortex"}};
enum{NUM_TEST_FUNCTIONS = sizeof(test_functions)/sizeof(test_functions[0])};

static void generate_interp_stacks(void);
static void free_interp_stacks(void);
static void LLtoXYZ(double lon, double lat, double p_out[]);

static int * read_mask(char const * filename, char const * grid_name, int const valid_mask_value);

#define STR_USAGE "Usage: %s -e experimentName\n"
#define YAC_ASSERT_ARGS(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(
  int argc, char ** argv, enum experiment_type * experiment);

int main (int argc, char *argv[]) {

  // Initialisation of MPI

  MPI_Init (0, NULL);

  /* Available experiment configurations
     "gen_mask" : all the grids, only interpolate the constant one
                  function with conservative 1st ord destarea
     "noicoh" : torc, nogt, bggd, sse7, icos grids, all the interpolations
     "nogt_icoh" : nogt, icoh grids, all the interpolations
     "icos_icoh" : icos, icoh grids, all the interpolations
  */

  // get experiment configuration
  enum experiment_type experiment = NOICOH; // default experiment
  parse_arguments(argc, argv, &experiment);
  struct experiment_configuration exp_config = experiment_configs[experiment];

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  YAC_ASSERT_F(
    (size_t)size == exp_config.num_components,
    "Wrong number of processes (has to be %zu)\n", exp_config.num_components);

  yac_cinit();
  yac_cdef_datetime("2008-03-09T16:05:07", "2008-03-10T16:05:07");
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);

  // generate interpolation stacks
  generate_interp_stacks();

  // configure interpolations
  for (size_t src_comp_idx = 0; src_comp_idx < exp_config.num_components;
       ++src_comp_idx) {
    char const * src_comp_name = exp_config.component_names[src_comp_idx];
    for (size_t tgt_comp_idx = src_comp_idx + 1;
         tgt_comp_idx < exp_config.num_components; ++tgt_comp_idx) {
      char const * tgt_comp_name = exp_config.component_names[tgt_comp_idx];
      for (size_t interp_idx = 0; interp_idx < exp_config.num_interpolations;
           ++interp_idx) {
        enum interp_stack_type interp_type =
          exp_config.interpolations[interp_idx];
        char field_name[256];
        sprintf(field_name, "%s_%s_out",
                src_comp_name, interp_stacks[interp_type].name);
        char weight_filename[1024];
        snprintf(
          weight_filename, 1024, "./output/rmp_%s_to_%s_yac_%s.nc",
          src_comp_name,
          exp_config.component_names[tgt_comp_idx],
          interp_stacks[interp_type].name);
        int ext_couple_config_id;
        yac_cget_ext_couple_config(&ext_couple_config_id);
        if (exp_config.write_weights)
          yac_cset_ext_couple_config_weight_file(
            ext_couple_config_id, weight_filename);
        yac_cdef_couple_custom(
          src_comp_name,                 // src_comp_name
          src_comp_name,                 // src_grid_name
          field_name,                    // src_field_name
          tgt_comp_name,                 // tgt_comp_name
          tgt_comp_name,                 // tgt_grid_name
          field_name,                    // tgt_field_name
          "10",                          // coupling period
          YAC_TIME_UNIT_SECOND,          // time unit
          YAC_REDUCTION_TIME_ACCUMULATE, // time reduction operation
          interp_stacks[interp_type].id, // interpolation stack id
          0, 0,                          // source and target lag
          ext_couple_config_id);         // additional configuration parameters
        yac_cfree_ext_couple_config(ext_couple_config_id);
      }
    }
  }

  // register local component
  int comp_id;
  char const * comp_name = exp_config.component_names[rank];
  yac_cdef_comp(comp_name, &comp_id);

  // read grid information
  char * grid_filename = "./input/grids.nc";
  char * mask_filename = "./input/masks_no_atm.nc";
  int valid_mask_value = 0;
  size_t num_vertices;
  size_t num_cells;
  int * num_vertices_per_cell;
  double * x_vertices;
  double * y_vertices;
  double * x_cell;
  double * y_cell;
  int * cell_to_vertex;
  int * cell_core_mask;
  yac_read_scrip_grid_information(
    grid_filename, mask_filename, comp_name, valid_mask_value,
    &num_vertices, &num_cells, &num_vertices_per_cell,
    &x_vertices, &y_vertices, &x_cell, &y_cell, &cell_to_vertex,
    &cell_core_mask, NULL, NULL, NULL);

  // register local grid
  int grid_id;
  yac_cdef_grid_unstruct(
    comp_name, (int)num_vertices, (int)num_cells, num_vertices_per_cell,
    x_vertices, y_vertices, cell_to_vertex, &grid_id);

  // set global ids and core masks
  yac_int * global_cell_ids = malloc(num_cells * sizeof(*global_cell_ids));
  yac_int * global_vertex_ids =
    malloc(num_vertices * sizeof(*global_vertex_ids));
  for (size_t i = 0; i < num_cells; ++i) global_cell_ids[i] = i;
  for (size_t i = 0; i < num_vertices; ++i) global_vertex_ids[i] = i;
  yac_cset_global_index(global_cell_ids, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);
  yac_cset_global_index(global_vertex_ids, YAC_LOCATION_CORNER, grid_id);

  int * nogt_mask = read_mask("./input/masks_nogt_yac.nc", comp_name, 0);
  int * torc_mask = read_mask("./input/masks_torc_yac.nc", comp_name, 0);

  // register cell middle points and mask
  int cell_point_id_no_mask,
      cell_point_id_nogt_mask,
      cell_point_id_torc_mask;
  yac_cdef_points_unstruct(
    grid_id, (int)num_cells, YAC_LOCATION_CELL, x_cell, y_cell,
    &cell_point_id_no_mask);
  if (nogt_mask != NULL) {
    yac_cdef_points_unstruct(
      grid_id, (int)num_cells, YAC_LOCATION_CELL, x_cell, y_cell,
      &cell_point_id_nogt_mask);
    yac_cset_mask(nogt_mask, cell_point_id_nogt_mask);
  } else cell_point_id_nogt_mask = cell_point_id_no_mask;
  if (torc_mask != NULL) {
    yac_cdef_points_unstruct(
      grid_id, (int)num_cells, YAC_LOCATION_CELL, x_cell, y_cell,
      &cell_point_id_torc_mask);
    yac_cset_mask(torc_mask, cell_point_id_torc_mask);
  } else cell_point_id_torc_mask = cell_point_id_no_mask;

  // register in fields
  int field_ids[exp_config.num_components][exp_config.num_interpolations];
  for (size_t src_comp_idx = 0; src_comp_idx < exp_config.num_components;
       ++src_comp_idx) {
    int cell_point_id = cell_point_id_no_mask;
    if (!strcmp("nogt", exp_config.component_names[src_comp_idx]))
      cell_point_id = cell_point_id_nogt_mask;
    if (!strcmp("torc", exp_config.component_names[src_comp_idx]))
      cell_point_id = cell_point_id_torc_mask;
    for (size_t interp_idx = 0; interp_idx < exp_config.num_interpolations;
         ++interp_idx) {
      enum interp_stack_type interp_type =
        exp_config.interpolations[interp_idx];
      char field_name[1024];
      snprintf(
        field_name, 1024, "%s_%s_out",
        exp_config.component_names[src_comp_idx],
        interp_stacks[interp_type].name);
      yac_cdef_field(
        field_name, comp_id, &cell_point_id, 1, 1, "10", YAC_TIME_UNIT_SECOND,
        &field_ids[src_comp_idx][interp_idx]);
    }
  }

  // initialise coupling
  yac_cenddef( );

  // initialise vtk output file
  YAC_VTK_FILE *vtk_file;
  {
    char vtk_filename[1024];
    sprintf(vtk_filename, "./output/%s.vtk", comp_name);

    double (*point_data)[3] = malloc(num_vertices * sizeof(*point_data));
    for (size_t i = 0; i < num_vertices; ++i)
      LLtoXYZ(x_vertices[i], y_vertices[i], point_data[i]);

    vtk_file = yac_vtk_open(vtk_filename, comp_name);
    yac_vtk_write_point_data(vtk_file, &point_data[0][0], (int)num_vertices);
    yac_vtk_write_cell_data(
      vtk_file, (unsigned *)cell_to_vertex,
      (unsigned*)num_vertices_per_cell, (int)num_cells);
    yac_vtk_write_point_scalars_int(
      vtk_file, global_vertex_ids, (int)num_vertices, "global_vertex_id");
    yac_vtk_write_cell_scalars_int(
      vtk_file, cell_core_mask, (int)num_cells, "cell_core_mask");
    yac_vtk_write_cell_scalars_int(
      vtk_file, global_cell_ids, (int)num_cells, "global_cell_id");
    if (nogt_mask != NULL)
      yac_vtk_write_cell_scalars_int(
        vtk_file, nogt_mask, (int)num_cells, "nogt_mask");
    if (torc_mask != NULL)
      yac_vtk_write_cell_scalars_int(
        vtk_file, torc_mask, (int)num_cells, "torc_mask");
    free(point_data);
  }

  if (nogt_mask) free(nogt_mask);
  if (torc_mask) free(torc_mask);

  // allocate memory for field data
  double * out_data = malloc(num_cells * sizeof(*out_data));
  double * in_data = malloc(num_cells * sizeof(*in_data));

  // "time loop"
  for (size_t test_idx = 0; test_idx < NUM_TEST_FUNCTIONS; ++test_idx) {

    int info, err;

    // generate out data
    for (size_t j = 0; j < num_cells; ++j)
      out_data[j] = test_functions[test_idx].p(x_cell[j], y_cell[j]);

    // write out field to vtk file
    {
      char field_name[1024];
      snprintf(field_name, 1024, "test_%s_out", test_functions[test_idx].name);
      yac_vtk_write_cell_scalars_double(
        vtk_file, out_data, (int)num_cells, field_name);
    }

    // put out fields
    double *point_set_data[1] = {out_data};
    double **collection_data[1] = {point_set_data};
    for (size_t tgt_comp_idx = 0; tgt_comp_idx < exp_config.num_interpolations;
         ++tgt_comp_idx)
      yac_cput(field_ids[rank][tgt_comp_idx], 1, collection_data, &info, &err);

    // get in fields
    for (size_t src_comp_idx = 0; src_comp_idx < exp_config.num_components;
         ++src_comp_idx) {
      if ((int)src_comp_idx == rank) continue;
      for (size_t interp_idx = 0; interp_idx < exp_config.num_interpolations;
           ++interp_idx) {
        enum interp_stack_type interp_type =
          exp_config.interpolations[interp_idx];
        for (size_t l = 0; l < num_cells; ++l) in_data[l] = 0.0; //-10.0;
        double *collection_data[1] = {in_data};
        yac_cget(
          field_ids[src_comp_idx][interp_idx], 1, collection_data, &info, &err);

        // write in field to vtk file
        {
          char field_name[1024];
          snprintf(
            field_name, 1024, "test_%s_%s_%s",
            test_functions[test_idx].name,
            exp_config.component_names[src_comp_idx],
            interp_stacks[interp_type].name);
          yac_vtk_write_cell_scalars_double(
            vtk_file, in_data, (int)num_cells, field_name);
        }
      }
    }
  }

  // free interpolation stacks
  free_interp_stacks();

  // close vtk file
  {
    yac_vtk_close(vtk_file);
  }

  // finalize YAC and MPI
  yac_cfinalize();
  MPI_Finalize();

  free(in_data);
  free(out_data);
  free(num_vertices_per_cell);
  free(x_vertices);
  free(y_vertices);
  free(x_cell);
  free(y_cell);
  free(cell_to_vertex);
  free(cell_core_mask);
  free(global_vertex_ids);
  free(global_cell_ids);

  return EXIT_SUCCESS;
}

static int * read_mask(
  char const * filename, char const * grid_name, int const valid_mask_value) {

  size_t grid_name_len = strlen(grid_name) + 1;
  char x_dim_name[2 + grid_name_len];
  char y_dim_name[2 + grid_name_len];
  char mask_var_name[4 + grid_name_len];

  snprintf(x_dim_name, 2 + grid_name_len, "x_%s", grid_name);
  snprintf(y_dim_name, 2 + grid_name_len, "y_%s", grid_name);
  snprintf(mask_var_name, 4 + grid_name_len, "%s.msk", grid_name);

  if (!yac_file_exists(filename)) return NULL;

  int ncid, status;
  yac_nc_open(filename, 0, &ncid);

  int mask_var_id;
  status = nc_inq_varid(ncid, mask_var_name, &mask_var_id);

  if (status == NC_ENOTVAR) {
    YAC_HANDLE_ERROR(nc_close(ncid));
    return NULL;
  } else if (status != NC_NOERR)
    YAC_HANDLE_ERROR(status);

  int x_dim_id;
  int y_dim_id;
  yac_nc_inq_dimid(ncid, x_dim_name, &x_dim_id);
  yac_nc_inq_dimid(ncid, y_dim_name, &y_dim_id);

  size_t x_dim_len;
  size_t y_dim_len;
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, x_dim_id, &x_dim_len));
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, y_dim_id, &y_dim_len));

  size_t num_cells = x_dim_len * y_dim_len;

  int * cell_mask = malloc(num_cells * sizeof(*cell_mask));
  YAC_HANDLE_ERROR(nc_get_var_int(ncid, mask_var_id, cell_mask));
  for (size_t i = 0; i < num_cells; ++i)
    cell_mask[i] = cell_mask[i] == valid_mask_value;


  YAC_HANDLE_ERROR(nc_close(ncid));

  return cell_mask;
}

static void parse_arguments(
  int argc, char ** argv, enum experiment_type * experiment) {

  int opt;
  while ((opt = getopt(argc, argv, "e:")) != -1) {
    YAC_ASSERT_ARGS((opt == 'e'), "invalid command argument")
    switch (opt) {
      default:
      case 'e':
        if (!strcmp(optarg, "GEN_MASK")) *experiment = GEN_MASK;
        else if (!strcmp(optarg, "NOICOH")) *experiment = NOICOH;
        else if (!strcmp(optarg, "NOGT_ICOH")) *experiment = NOGT_ICOH;
        else if (!strcmp(optarg, "ICOS_ICOH")) *experiment = ICOS_ICOH;
        else *experiment = NUM_EXPERIMENT;
        YAC_ASSERT_ARGS(
          *experiment != NUM_EXPERIMENT, "invalid experiment name")
        break;
    }
  }
}

static void generate_interp_stacks(void) {
  for (int interp_stack_type = 0; interp_stack_type < NUM_INTERP_STACK_TYPES;
       ++interp_stack_type) {
    switch (interp_stack_type) {
      default:
      case(INTERP_STACK_CONSERV_DESTAREA):
        interp_stacks[interp_stack_type].name = strdup("CONSERV_DESTAREA");
        yac_cget_interp_stack_config(&(interp_stacks[interp_stack_type].id));
        yac_cadd_interp_stack_config_conservative(
          interp_stacks[interp_stack_type].id, 1, 1, 1, YAC_CONSERV_DESTAREA);
        yac_cadd_interp_stack_config_fixed(
          interp_stacks[interp_stack_type].id, 0.0);
        break;
      case(INTERP_STACK_CONSERV_FRACAREA):
        interp_stacks[interp_stack_type].name = strdup("CONSERV_FRACAREA");
        yac_cget_interp_stack_config(&(interp_stacks[interp_stack_type].id));
        yac_cadd_interp_stack_config_conservative(
          interp_stacks[interp_stack_type].id, 1, 1, 1, YAC_CONSERV_FRACAREA);
        yac_cadd_interp_stack_config_fixed(
          interp_stacks[interp_stack_type].id, 0.0);
        break;
      case(INTERP_STACK_CONS2ND_FRACAREA):
        interp_stacks[interp_stack_type].name = strdup("CONS2ND_FRACAREA");
        yac_cget_interp_stack_config(&(interp_stacks[interp_stack_type].id));
        yac_cadd_interp_stack_config_conservative(
          interp_stacks[interp_stack_type].id, 2, 0, 1, YAC_CONSERV_FRACAREA);
        yac_cadd_interp_stack_config_fixed(
          interp_stacks[interp_stack_type].id, 0.0);
        break;
      case(INTERP_STACK_DISTWGT_4):
        interp_stacks[interp_stack_type].name = strdup("DISTWGT_4");
        yac_cget_interp_stack_config(&(interp_stacks[interp_stack_type].id));
        yac_cadd_interp_stack_config_nnn(
          interp_stacks[interp_stack_type].id, YAC_NNN_DIST, 4,
          YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT, -1.0);
        break;
      case(INTERP_STACK_DISTWGT_1):
        interp_stacks[interp_stack_type].name = strdup("DISTWGT_1");
        yac_cget_interp_stack_config(&(interp_stacks[interp_stack_type].id));
        yac_cadd_interp_stack_config_nnn(
          interp_stacks[interp_stack_type].id, YAC_NNN_AVG, 1,
          YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT, -1.0);
        break;
      case(INTERP_STACK_HCSBB):
        interp_stacks[interp_stack_type].name = strdup("HCSBB");
        yac_cget_interp_stack_config(&(interp_stacks[interp_stack_type].id));
        yac_cadd_interp_stack_config_hcsbb(interp_stacks[interp_stack_type].id);
        yac_cadd_interp_stack_config_nnn(
          interp_stacks[interp_stack_type].id, YAC_NNN_DIST, 4,
          YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT, -1.0);
        yac_cadd_interp_stack_config_fixed(
          interp_stacks[interp_stack_type].id, 0.0);
        break;
      case(INTERP_STACK_AVG_DIST):
        interp_stacks[interp_stack_type].name = strdup("AVG_DIST");
        yac_cget_interp_stack_config(&(interp_stacks[interp_stack_type].id));
        yac_cadd_interp_stack_config_average(
          interp_stacks[interp_stack_type].id, YAC_AVG_ARITHMETIC, 1);
        yac_cadd_interp_stack_config_nnn(
          interp_stacks[interp_stack_type].id, YAC_NNN_DIST, 2,
          YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT, -1.0);
        break;
    }
  }
}

static void free_interp_stacks(void) {
  for (int interp_stack_type = 0; interp_stack_type < NUM_INTERP_STACK_TYPES;
       ++interp_stack_type) {
    free((void*)(interp_stacks[interp_stack_type].name));
    yac_cfree_interp_stack_config(interp_stacks[interp_stack_type].id);
  }
}


static void LLtoXYZ(double lon, double lat, double p_out[]) {

   while (lon < -M_PI) lon += 2.0 * M_PI;
   while (lon >= M_PI) lon -= 2.0 * M_PI;

   double cos_lat = cos(lat);
   p_out[0] = cos_lat * cos(lon);
   p_out[1] = cos_lat * sin(lon);
   p_out[2] = sin(lat);
}
