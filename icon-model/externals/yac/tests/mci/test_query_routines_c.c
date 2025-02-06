// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

#include "tests.h"

#include "yac.h"
#include "yaxt.h"

void yac_cget_comp_size_c2py(int comp_id, int* size);
void yac_cget_comp_rank_c2py(int comp_id, int* rank);

int main(void) {

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);

  int yac_id;
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cinit();
  yac_cinit_instance(&yac_id);

  int global_rank, global_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &global_size);

  //--------------------------------------------------
  // definition of components, grid, points and fields
  //--------------------------------------------------

  int rank_type = global_rank%5;

  if (rank_type == 0) {

    // ranks of type 0 define no data
    yac_cdef_comps(NULL, 0, NULL);
    yac_cdef_comps_instance(yac_id, NULL, 0, NULL);

  } else {

    // ranks of type 1 or higher define a component
    char comp_name[16];
    int comp_id, comp_id_instance;
    sprintf(comp_name, "comp_%d", global_rank);
    yac_cdef_comp(comp_name, &comp_id);
    yac_cdef_comp_instance(yac_id, comp_name, &comp_id_instance);

    // test C to Python interfaces
    {
      int size, rank;

      yac_cget_comp_size_c2py(comp_id, &size);
      if (size != 1) PUT_ERR("error in yac_cget_comp_size_c2py");

      yac_cget_comp_rank_c2py(comp_id, &rank);
      if (rank != 0) PUT_ERR("erro in yac_cget_comp_rank_c2py");
    }

    // ranks of type 2 or higher define a grid and attach some
    // meta data to their component
    if (rank_type > 1) {

      char comp_meta_data[64];
      sprintf(comp_meta_data, "%s METADATA_C", comp_name);
      yac_cdef_component_metadata(comp_name, comp_meta_data);
      sprintf(comp_meta_data, "%s METADATA_C_instance", comp_name);
      yac_cdef_component_metadata_instance(yac_id, comp_name, comp_meta_data);

      char grid_name[16];
      int grid_id;
      sprintf(grid_name, "grid_%d_0", global_rank);
      yac_cdef_grid_reg2d(
        grid_name, (int[]){2, 2}, (int[]){0,0},
        (double[]){0.,180.}, (double[]){-45.,45.}, &grid_id);

      int point_id;
      yac_cdef_points_reg2d(
        grid_id, (int[]){1,1}, YAC_LOCATION_CELL,
        (double[]){90.}, (double[]){0.}, &point_id);

      // ranks of type 3 or higher define a second grid, attach some
      // meta data to it
      if (rank_type > 2) {

        sprintf(grid_name, "grid_%d_1", global_rank);
        yac_cdef_grid_reg2d(
          grid_name, (int[]){2, 2}, (int[]){0,0},
          (double[]){0.,180.}, (double[]){-45.,45.}, &grid_id);

        yac_cdef_points_reg2d(
          grid_id, (int[]){1,1}, YAC_LOCATION_CELL,
          (double[]){90.}, (double[]){0.}, &point_id);

        char grid_meta_data[64];
        sprintf(grid_meta_data, "%s METADATA_G", grid_name);
        yac_cdef_grid_metadata(grid_name, grid_meta_data);
        sprintf(grid_meta_data, "%s METADATA_G_instance", grid_name);
        yac_cdef_grid_metadata_instance(yac_id, grid_name, grid_meta_data);

        // ranks of type 4 or higher define a third grid, attach some
        // meta data to it and define two fields (one with meta data
        // and one without)
        if (rank_type > 3) {

          sprintf(grid_name, "grid_%d_2", global_rank);
          yac_cdef_grid_reg2d(
            grid_name, (int[]){2, 2}, (int[]){0,0},
            (double[]){0.,180.}, (double[]){-45.,45.}, &grid_id);

          yac_cdef_points_reg2d(
            grid_id, (int[]){1,1}, YAC_LOCATION_CELL,
            (double[]){90.}, (double[]){0.}, &point_id);

          char grid_meta_data[64];
          sprintf(grid_meta_data, "%s METADATA_G", grid_name);
          yac_cdef_grid_metadata(grid_name, grid_meta_data);
          sprintf(grid_meta_data, "%s METADATA_G_instance", grid_name);
          yac_cdef_grid_metadata_instance(yac_id, grid_name, grid_meta_data);

          char field_name[16];
          int field_id, field_id_instance;
          int point_ids[1];
          point_ids[0] = point_id;
          sprintf(field_name, "field_%d_0", global_rank);
          yac_cdef_field(
            field_name, comp_id, point_ids, 1, 1, "PT5M",
            YAC_TIME_UNIT_ISO_FORMAT, &field_id);
          yac_cdef_field(
            field_name, comp_id_instance, point_ids, 1, 1, "PT5M",
            YAC_TIME_UNIT_ISO_FORMAT, &field_id_instance);

          sprintf(field_name, "field_%d_1", global_rank);
          yac_cdef_field(
            field_name, comp_id, point_ids, 1, 1, "PT5M",
            YAC_TIME_UNIT_ISO_FORMAT, &field_id);
          yac_cdef_field(
            field_name, comp_id_instance, point_ids, 1, 1, "PT5M",
            YAC_TIME_UNIT_ISO_FORMAT, &field_id_instance);

          char field_meta_data[64];
          sprintf(field_meta_data, "%s METADATA_F", field_name);
          yac_cdef_field_metadata(
            comp_name, grid_name, field_name, field_meta_data);
          sprintf(field_meta_data, "%s METADATA_F_instance", field_name);
          yac_cdef_field_metadata_instance(
            yac_id, comp_name, grid_name, field_name, field_meta_data);

          //----------------------------------------
          // test some field_id based query routines
          //----------------------------------------
          if (field_id != yac_cget_field_id(comp_name, grid_name, field_name))
            PUT_ERR("ERROR in yac_cget_field_id");
          if (field_id_instance !=
              yac_cget_field_id_instance(
                yac_id, comp_name, grid_name, field_name))
            PUT_ERR("ERROR in yac_cget_field_id_instance");
          if (strcmp(comp_name, yac_cget_component_name_from_field_id(field_id)))
            PUT_ERR("ERROR in yac_cget_component_name_from_field_id");
          if (strcmp(comp_name,
                     yac_cget_component_name_from_field_id(field_id_instance)))
            PUT_ERR("ERROR in yac_cget_component_name_from_field_id");
          if (strcmp(grid_name, yac_cget_grid_name_from_field_id(field_id)))
            PUT_ERR("ERROR in yac_cget_grid_name_from_field_id");
          if (strcmp(grid_name,
                     yac_cget_grid_name_from_field_id(field_id_instance)))
            PUT_ERR("ERROR in yac_cget_grid_name_from_field_id");
          if (strcmp(field_name, yac_cget_field_name_from_field_id(field_id)))
            PUT_ERR("ERROR in yac_cget_field_name_from_field_id");
          if (strcmp(field_name,
                     yac_cget_field_name_from_field_id(field_id_instance)))
            PUT_ERR("ERROR in yac_cget_field_name_from_field_id");
          if (strcmp("PT5M", yac_cget_timestep_from_field_id(field_id)))
            PUT_ERR("ERROR in yac_cget_timestep_from_field_id");
          if (strcmp("PT5M",
                     yac_cget_timestep_from_field_id(field_id_instance)))
            PUT_ERR("ERROR in yac_cget_timestep_from_field_id");
          if (1 != yac_cget_collection_size_from_field_id(field_id))
            PUT_ERR("ERROR in yac_cget_collection_size_from_field_id");
          if (1 != yac_cget_collection_size_from_field_id(field_id_instance))
            PUT_ERR("ERROR in yac_cget_collection_size_from_field_id");
          if (YAC_EXCHANGE_TYPE_NONE !=
              yac_cget_role_from_field_id(field_id))
            PUT_ERR("ERROR in yac_cget_role_from_field_id");
          if (YAC_EXCHANGE_TYPE_NONE !=
              yac_cget_role_from_field_id(field_id_instance))
            PUT_ERR("ERROR in yac_cget_role_from_field_id");
        }
      }
    }
  }

  // end definition phase and synchronise definitions across all processes
  yac_csync_def();
  yac_csync_def_instance(yac_id);

  //--------------------
  // test query routines
  //--------------------

  int const ref_nbr_comps = 4 * ((global_size - 1) / 5) + ((global_size - 1) % 5);

  if (ref_nbr_comps != yac_cget_nbr_comps())
    PUT_ERR("ERROR in yac_cget_nbr_comps");
  if (ref_nbr_comps != yac_cget_nbr_comps_instance(yac_id))
    PUT_ERR("ERROR in yac_cget_nbr_comps_instance");

  char const ** comp_names =
    malloc((size_t)ref_nbr_comps * sizeof(*comp_names));
  char const ** comp_names_instance =
    malloc((size_t)ref_nbr_comps * sizeof(*comp_names_instance));
  yac_cget_comp_names(ref_nbr_comps, comp_names);
  yac_cget_comp_names_instance(yac_id, ref_nbr_comps, comp_names_instance);

  int const ref_nbr_grids =
    3 * ((global_size - 1) / 5) +
    ((global_size - 1) % 5 > 2) +
    2 * ((global_size - 1) % 5 > 3);

  if (ref_nbr_grids != yac_cget_nbr_grids())
    PUT_ERR("ERROR in yac_cget_nbr_grids");
  if (ref_nbr_grids != yac_cget_nbr_grids_instance(yac_id))
    PUT_ERR("ERROR in yac_cget_nbr_grids_instance");

  char const ** grid_names =
    malloc((size_t)ref_nbr_grids * sizeof(*grid_names));
  char const ** grid_names_instance =
    malloc((size_t)ref_nbr_grids * sizeof(*grid_names_instance));
  yac_cget_grid_names(ref_nbr_grids, grid_names);
  yac_cget_grid_names_instance(yac_id, ref_nbr_grids, grid_names_instance);

  for (int rank = 0; rank < global_size; ++rank) {

    int idx;
    int rank_type = rank % 5;

    char ref_comp_name[16];
    char ref_grid_names[3][32];
    char ref_field_names[2][32];
    sprintf(ref_comp_name, "comp_%d", rank);
    for (int i = 0; i < 3; ++i)
      sprintf(ref_grid_names[i], "grid_%d_%d", rank, i);
    for (int i = 0; i < 2; ++i)
      sprintf(ref_field_names[i], "field_%d_%d", rank, i);

    if (rank_type == 0) {

      // no component or grid is defined on this rank
      for (idx = 0; idx < ref_nbr_comps; ++idx)
        if (!strcmp(ref_comp_name, comp_names[idx])) break;
      if (idx != ref_nbr_comps) PUT_ERR("ERROR in yac_cget_comp_names");
      for (idx = 0; idx < ref_nbr_comps; ++idx)
        if (!strcmp(ref_comp_name, comp_names_instance[idx])) break;
      if (idx != ref_nbr_comps) PUT_ERR("ERROR in yac_cget_comp_names_instance");
      for (int i = 0; i < 3; ++i) {
        for (idx = 0; idx < ref_nbr_grids; ++idx)
          if (!strcmp(ref_grid_names[i], grid_names[idx])) break;
        if (idx != ref_nbr_grids) PUT_ERR("ERROR in yac_cget_grid_names");
        for (idx = 0; idx < ref_nbr_grids; ++idx)
          if (!strcmp(ref_grid_names[i], grid_names_instance[idx])) break;
        if (idx != ref_nbr_grids) PUT_ERR("ERROR in yac_cget_grid_names_instance");
      }

    } else {

      // a component is defined on this process
      for (idx = 0; idx < ref_nbr_comps; ++idx)
        if (!strcmp(ref_comp_name, comp_names[idx])) break;
      if (idx == ref_nbr_comps) PUT_ERR("ERROR in yac_cget_comp_names");
      for (idx = 0; idx < ref_nbr_comps; ++idx)
        if (!strcmp(ref_comp_name, comp_names_instance[idx])) break;
      if (idx == ref_nbr_comps) PUT_ERR("ERROR in yac_cget_comp_names_instance");

      if (rank_type == 1) {

        // the component of this rank has no meta data
        if (yac_cget_component_metadata(ref_comp_name))
          PUT_ERR("ERROR in yac_cget_component_metadata");
        if (yac_cget_component_metadata_instance(yac_id, ref_comp_name))
          PUT_ERR("ERROR in yac_cget_component_metadata_instance");

        // no grid is registered in the coupling configuration on this rank
        for (int i = 0; i < 3; ++i) {
          for (idx = 0; idx < ref_nbr_grids; ++idx)
            if (!strcmp(ref_grid_names[i], grid_names[idx])) break;
          if (idx != ref_nbr_grids) PUT_ERR("ERROR in yac_cget_grid_names");
          for (idx = 0; idx < ref_nbr_grids; ++idx)
            if (!strcmp(ref_grid_names[i], grid_names_instance[idx])) break;
          if (idx != ref_nbr_grids) PUT_ERR("ERROR in yac_cget_grid_names_instance");
        }

        if (yac_cget_comp_nbr_grids(ref_comp_name) != 0)
          PUT_ERR("ERROR in yac_cget_comp_nbr_grids");
        if (yac_cget_comp_nbr_grids_instance(yac_id, ref_comp_name) != 0)
          PUT_ERR("ERROR in yac_cget_comp_nbr_grids_instance");

      } else {

        char ref_comp_meta_data[64], ref_comp_meta_data_instance[64];
        sprintf(ref_comp_meta_data, "%s METADATA_C", ref_comp_name);
        sprintf(ref_comp_meta_data_instance, "%s METADATA_C_instance",
                ref_comp_name);

        // the component has meta data on this rank
        if (strcmp(
              ref_comp_meta_data,
              yac_cget_component_metadata(ref_comp_name)))
          PUT_ERR("ERROR in yac_cget_component_metadata");
        if (strcmp(
              ref_comp_meta_data_instance,
              yac_cget_component_metadata_instance(yac_id, ref_comp_name)))
          PUT_ERR("ERROR in yac_cget_component_metadata_instance");

        if (rank_type == 2) {

          // no grid is registered in the coupling configuration on this rank
          for (int i = 0; i < 3; ++i) {
            for (idx = 0; idx < ref_nbr_grids; ++idx)
              if (!strcmp(ref_grid_names[i], grid_names[idx])) break;
            if (idx != ref_nbr_grids) PUT_ERR("ERROR in yac_cget_grid_names");
            for (idx = 0; idx < ref_nbr_grids; ++idx)
              if (!strcmp(ref_grid_names[i], grid_names_instance[idx])) break;
            if (idx != ref_nbr_grids)
              PUT_ERR("ERROR in yac_cget_grid_names_instance");
          }

          if (yac_cget_comp_nbr_grids(ref_comp_name) != 0)
            PUT_ERR("ERROR in yac_cget_comp_nbr_grids");
          if (yac_cget_comp_nbr_grids_instance(yac_id, ref_comp_name) != 0)
            PUT_ERR("ERROR in yac_cget_comp_nbr_grids_instance");

        } else {

          // grid_*_0 is not registered in the coupling configuration
          for (idx = 0; idx < ref_nbr_grids; ++idx)
            if (!strcmp(ref_grid_names[0], grid_names[idx])) break;
          if (idx != ref_nbr_grids) PUT_ERR("ERROR in yac_cget_grid_names");
          for (idx = 0; idx < ref_nbr_grids; ++idx)
            if (!strcmp(ref_grid_names[0], grid_names_instance[idx])) break;
          if (idx != ref_nbr_grids)
            PUT_ERR("ERROR in yac_cget_grid_names_instance");

          char ref_grid_meta_data[2][64];
          char ref_grid_meta_data_instance[2][64];
          for (int i = 0; i < 2; ++i) {
            sprintf(ref_grid_meta_data[i], "%s METADATA_G",
                    ref_grid_names[i+1]);
            sprintf(ref_grid_meta_data_instance[i], "%s METADATA_G_instance",
                    ref_grid_names[i+1]);
          }

          // grid_*_1 is registered in the coupling configuration, has meta
          // data, but no fields
          for (idx = 0; idx < ref_nbr_grids; ++idx)
            if (!strcmp(ref_grid_names[1], grid_names[idx])) break;
          if (idx == ref_nbr_grids) PUT_ERR("ERROR in yac_cget_grid_names");
          for (idx = 0; idx < ref_nbr_grids; ++idx)
            if (!strcmp(ref_grid_names[1], grid_names_instance[idx])) break;
          if (idx == ref_nbr_grids)
            PUT_ERR("ERROR in yac_cget_grid_names_instance");
          if (strcmp(
                ref_grid_meta_data[0],
                yac_cget_grid_metadata(ref_grid_names[1])))
            PUT_ERR("ERROR in yac_cget_grid_metadata");
          if (strcmp(
                ref_grid_meta_data_instance[0],
                yac_cget_grid_metadata_instance(yac_id, ref_grid_names[1])))
            PUT_ERR("ERROR in yac_cget_grid_metadata_instance");
          if (yac_cget_nbr_fields(ref_comp_name, ref_grid_names[1]))
            PUT_ERR("ERROR in yac_cget_nbr_fields");
          if (yac_cget_nbr_fields_instance(
                yac_id, ref_comp_name, ref_grid_names[1]))
            PUT_ERR("ERROR in yac_cget_nbr_fields_instance");

          if (rank_type == 3) {

            // grid_*_2 is not defined on this rank
            for (idx = 0; idx < ref_nbr_grids; ++idx)
              if (!strcmp(ref_grid_names[2], grid_names[idx])) break;
            if (idx != ref_nbr_grids) PUT_ERR("ERROR in yac_cget_grid_names");
            for (idx = 0; idx < ref_nbr_grids; ++idx)
              if (!strcmp(ref_grid_names[2], grid_names_instance[idx])) break;
            if (idx != ref_nbr_grids)
              PUT_ERR("ERROR in yac_cget_grid_names_instance");

            if (yac_cget_comp_nbr_grids(ref_comp_name) != 0)
              PUT_ERR("ERROR in yac_cget_comp_nbr_grids");
            if (yac_cget_comp_nbr_grids_instance(yac_id, ref_comp_name) != 0)
              PUT_ERR("ERROR in yac_cget_comp_nbr_grids_instance");

          } else {

            // grid_*_2 is registered in the coupling configuration, has
            // meta data, and two fields
            for (idx = 0; idx < ref_nbr_grids; ++idx)
              if (!strcmp(ref_grid_names[2], grid_names[idx])) break;
            if (idx == ref_nbr_grids) PUT_ERR("ERROR in yac_cget_grid_names");
            for (idx = 0; idx < ref_nbr_grids; ++idx)
              if (!strcmp(ref_grid_names[2], grid_names_instance[idx])) break;
            if (idx == ref_nbr_grids)
              PUT_ERR("ERROR in yac_cget_grid_names_instance");
            if (strcmp(
                  ref_grid_meta_data[1],
                  yac_cget_grid_metadata(ref_grid_names[2])))
              PUT_ERR("ERROR in yac_cget_grid_metadata");
            if (strcmp(
                  ref_grid_meta_data_instance[1],
                  yac_cget_grid_metadata_instance(yac_id, ref_grid_names[2])))
              PUT_ERR("ERROR in yac_cget_grid_metadata_instance");

            int const ref_nbr_fields = 2;
            char const * field_names[2];
            char const * field_names_instance[2];
            char ref_field_meta_data[64];
            char ref_field_meta_data_instance[64];
            sprintf(ref_field_meta_data, "%s METADATA_F", ref_field_names[1]);
            sprintf(ref_field_meta_data_instance, "%s METADATA_F_instance",
                    ref_field_names[1]);

            // field_*_0 and field_*_1 are defined and the second one has
            // meta data
            if (yac_cget_nbr_fields(
                  ref_comp_name, ref_grid_names[2]) != ref_nbr_fields)
              PUT_ERR("ERROR in yac_cget_nbr_fields");
            if (yac_cget_nbr_fields_instance(
                  yac_id, ref_comp_name, ref_grid_names[2]) != ref_nbr_fields)
              PUT_ERR("ERROR in yac_cget_nbr_fields_instance");
            yac_cget_field_names(
              ref_comp_name, ref_grid_names[2], ref_nbr_fields, field_names);
            yac_cget_field_names_instance(
              yac_id, ref_comp_name, ref_grid_names[2], ref_nbr_fields,
              field_names_instance);
            for (int i = 0; i < ref_nbr_fields; ++i) {
              for (idx = 0; idx < ref_nbr_fields  ; ++idx)
                if (!strcmp(ref_field_names[i], field_names[idx])) break;
              if (idx == ref_nbr_fields)
                PUT_ERR("ERROR in yac_cget_field_names");
              for (idx = 0; idx < ref_nbr_fields; ++idx)
                if (!strcmp(
                       ref_field_names[i], field_names_instance[idx])) break;
              if (idx == ref_nbr_fields)
                PUT_ERR("ERROR in yac_cget_field_names_instance");
            }
            if (yac_cget_field_metadata(
                  ref_comp_name, ref_grid_names[2], ref_field_names[0]))
              PUT_ERR("ERROR in yac_cget_field_metadata");
            if (yac_cget_field_metadata_instance(
                  yac_id, ref_comp_name, ref_grid_names[2],
                  ref_field_names[0]))
              PUT_ERR("ERROR in yac_cget_field_metadata_instance");
            if (strcmp(
                  ref_field_meta_data,
                  yac_cget_field_metadata(
                    ref_comp_name, ref_grid_names[2], ref_field_names[1])))
              PUT_ERR("ERROR in yac_cget_field_metadata");
            if (strcmp(
                  ref_field_meta_data_instance,
                  yac_cget_field_metadata_instance(
                    yac_id, ref_comp_name, ref_grid_names[2],
                    ref_field_names[1])))
              PUT_ERR("ERROR in yac_cget_field_metadata_instance");

            if (yac_cget_comp_nbr_grids(ref_comp_name) != 1)
              PUT_ERR("ERROR in yac_cget_comp_nbr_grids");
            if (yac_cget_comp_nbr_grids_instance(yac_id, ref_comp_name) != 1)
              PUT_ERR("ERROR in yac_cget_comp_nbr_grids_instance");
            char const * comp_grid_name;
            yac_cget_comp_grid_names(ref_comp_name, 1, &comp_grid_name);
            if (strcmp(comp_grid_name, ref_grid_names[2]))
              PUT_ERR("ERROR in yac_cget_comp_grid_names");
            yac_cget_comp_grid_names_instance(yac_id, ref_comp_name, 1, &comp_grid_name);
            if (strcmp(comp_grid_name, ref_grid_names[2]))
              PUT_ERR("ERROR in yac_cget_comp_grid_names_instance");

            if (yac_cget_field_role(
                  ref_comp_name, ref_grid_names[2],
                  ref_field_names[0]) != YAC_EXCHANGE_TYPE_NONE)
              PUT_ERR("ERROR in yac_cget_field_role");
            if (yac_cget_field_role_instance(
                  yac_id, ref_comp_name, ref_grid_names[2],
                  ref_field_names[0]) != YAC_EXCHANGE_TYPE_NONE)
              PUT_ERR("ERROR in yac_cget_field_role_instance");
            if (yac_cget_field_collection_size(
                  ref_comp_name, ref_grid_names[2],
                  ref_field_names[0]) != 1)
              PUT_ERR("ERROR in yac_cget_field_collection_size");
            if (yac_cget_field_collection_size_instance(
                  yac_id, ref_comp_name, ref_grid_names[2],
                  ref_field_names[0]) != 1)
              PUT_ERR("ERROR in yac_cget_field_collection_size_instance");
            if (strcmp(
                  yac_cget_field_timestep(
                    ref_comp_name, ref_grid_names[2],
                    ref_field_names[0]), "PT5M"))
              PUT_ERR("ERROR in yac_cget_field_timestep");
            if (strcmp(
                  yac_cget_field_timestep_instance(
                    yac_id, ref_comp_name, ref_grid_names[2],
                    ref_field_names[0]), "PT5M"))
              PUT_ERR("ERROR in yac_cget_field_timestep");
          }
        }
      }
    }
  }

  free(grid_names_instance);
  free(grid_names);
  free(comp_names_instance);
  free(comp_names);

  yac_cfinalize_instance(yac_id);
  yac_cfinalize();
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}
