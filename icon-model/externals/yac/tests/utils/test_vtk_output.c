// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <unistd.h>
#include <limits.h>

#include "vtk_output.h"
#include "generate_cubed_sphere.h"
#include "test_function.h"
#include "geometry.h"

int main (void) {

  {
    // initialise and open vtk file
    char const filename[] = "test_vtk_output_1.vtk";
    YAC_VTK_FILE * vtk_file = yac_vtk_open(filename, "test");

    enum {
      NUM_POINTS = 13,
      LON = 0,
      LAT = 1,
    };
    double point_data[2][NUM_POINTS] = {{0,
                                         45,165,285,
                                         0,90,180,270,
                                         180,252,324,396,468},
                                        {90,
                                         85,85,85,
                                         80,80,80,80,
                                         75,75,75,75,75}};
    for (size_t i = 0; i < NUM_POINTS; ++i) {
      point_data[LON][i] *= YAC_RAD;
      point_data[LAT][i] *= YAC_RAD;
    }
    unsigned cell_corners[NUM_POINTS] = {0,
                                         1,2,3,
                                         4,5,6,7,
                                         8,9,10,11,12};
    unsigned num_points_per_cell[] = {1,3,4,5};
    unsigned num_cells = 4;

    // writes point data to file
    yac_vtk_write_point_data_ll(
      vtk_file, point_data[LON], point_data[LAT], NUM_POINTS);
    yac_vtk_write_cell_data(
      vtk_file, cell_corners, num_points_per_cell, num_cells);

    unsigned cell_scalars_uint[] = {0,1,2};
    int cell_scalars_int[] = {0,1,2};
    float cell_scalars_float[] = {0,1,2};
    double cell_scalars_double[] = {0,1,2};
    yac_vtk_write_cell_scalars_uint(
      vtk_file, cell_scalars_uint, num_cells, "cell_scalars_uint");
    yac_vtk_write_cell_scalars_int(
      vtk_file, cell_scalars_int, num_cells, "cell_scalars_int");
    yac_vtk_write_cell_scalars_float(
      vtk_file, cell_scalars_float, num_cells, "cell_scalars_float");
    yac_vtk_write_cell_scalars_double(
      vtk_file, cell_scalars_double, num_cells, "cell_scalars_double");

    unsigned point_scalars_uint[NUM_POINTS] = {UINT_MAX, 0,1,2,3,4,5,6,7,8,9,10,11};
    int point_scalars_int[NUM_POINTS] = {-1, 0,1,2,3,4,5,6,7,8,9,10,11};
    float point_scalars_float[NUM_POINTS] = {-1, 0,1,2,3,4,5,6,7,8,9,10,11};
    double point_scalars_double[NUM_POINTS];
    for (size_t i = 0; i < NUM_POINTS; ++i)
      point_scalars_double[i] =
        yac_test_func_deg(point_data[LON][i], point_data[LAT][i]);
    yac_vtk_write_point_scalars_uint(
      vtk_file, point_scalars_uint, NUM_POINTS, "point_scalars_uint");
    yac_vtk_write_point_scalars_int(
      vtk_file, point_scalars_int, NUM_POINTS, "point_scalars_int");
    yac_vtk_write_point_scalars_float(
      vtk_file, point_scalars_float, NUM_POINTS, "point_scalars_float");
    yac_vtk_write_point_scalars_double(
      vtk_file, point_scalars_double, NUM_POINTS, "point_scalars_double");

    // closes vtk file
    yac_vtk_close(vtk_file);

    unlink(filename);
  }

  {
    unsigned n = 30;
    unsigned num_cells;
    unsigned num_vertices;
    double * x_vertices;
    double * y_vertices;
    double * z_vertices;
    unsigned * vertices_of_cell;
    unsigned * face_id;

    yac_generate_cubed_sphere_grid_information(
      n, &num_cells, &num_vertices, &x_vertices, &y_vertices, &z_vertices,
      &vertices_of_cell, &face_id);

    struct {
      char * func_name;
      double (*func)(double, double);
    } test_funcs[] = {
        {.func_name = "yac_test_func", .func = yac_test_func},
        {.func_name = "yac_test_ana_fcos", .func = yac_test_ana_fcos},
        {.func_name = "yac_test_ana_fcossin", .func = yac_test_ana_fcossin},
        {.func_name = "yac_test_one", .func = yac_test_one},
        {.func_name = "yac_test_gulfstream", .func = yac_test_gulfstream},
        {.func_name = "yac_test_harmonic", .func = yac_test_harmonic},
        {.func_name = "yac_test_vortex", .func = yac_test_vortex}
      };

    // initialise and open vtk file
    char const filename[] = "test_vtk_output_2.vtk";
    YAC_VTK_FILE * vtk_file = yac_vtk_open(filename, "test2");

    double * point_data = xmalloc(3 * num_vertices * sizeof(*point_data));
    for (unsigned i = 0; i < num_vertices; ++i) {
      point_data[3 * i + 0] = x_vertices[i];
      point_data[3 * i + 1] = y_vertices[i];
      point_data[3 * i + 2] = z_vertices[i];
    }
    unsigned * num_points_per_cell =
      xmalloc(num_cells * sizeof(*num_points_per_cell));
    for (unsigned i = 0; i < num_cells; ++i) num_points_per_cell[i] = 4;

    // writes point data to file
    yac_vtk_write_point_data(vtk_file, point_data, num_vertices);
    yac_vtk_write_cell_data(
      vtk_file, vertices_of_cell, num_points_per_cell, num_cells);

    double * field_data = xmalloc(num_vertices * sizeof(*field_data));
    for (size_t i = 0; i < sizeof(test_funcs) / sizeof(test_funcs[0]); ++i) {

      for (unsigned j = 0; j < num_vertices; ++j) {

        double lon, lat;
        XYZtoLL(point_data + 3 * j, &lon, &lat);

        field_data[j] = test_funcs[i].func(lon, lat);
      }

      yac_vtk_write_point_scalars_double(
        vtk_file, field_data, num_vertices, test_funcs[i].func_name);
    }

    free(field_data);
    free(num_points_per_cell);
    free(point_data);

    // closes vtk file
    yac_vtk_close(vtk_file);

    unlink(filename);

    free(x_vertices);
    free(y_vertices);
    free(z_vertices);
    free(vertices_of_cell);
    free(face_id);
  }

  return EXIT_SUCCESS;
}
