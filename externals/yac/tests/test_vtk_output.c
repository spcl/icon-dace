/**
 * @file test_vtk_output.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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
#include <unistd.h>

#include "vtk_output.h"
#include "generate_cubed_sphere.h"
#include "test_function.h"
#include "geometry.h"

int main (void) {

  {
    // initialise and open vtk file
    char const filename[] = "test_vtk_output.vtk";
    VTK_FILE * vtk_file = vtk_open(filename, "test");

    size_t num_points = 3+4+5;
    double point_data_lon[] = {45,165,285,
                               0,90,180,270,
                               180,252,324,396,468};
    double point_data_lat[] = {85,85,85,
                               80,80,80,80,
                               75,75,75,75,75};
    double point_data[3+4+5][3];
    for (size_t i = 0; i < num_points; ++i)
      LLtoXYZ_deg(point_data_lon[i], point_data_lat[i], point_data[i]);
    unsigned cell_corners[] = {0,1,2,
                               3,4,5,6,
                               7,8,9,10,11};
    unsigned num_points_per_cell[] = {3,4,5};
    unsigned num_cells = 3;

    // writes point data to file
    vtk_write_point_data(vtk_file, &point_data[0][0], num_points);
    vtk_write_cell_data(vtk_file, cell_corners, num_points_per_cell, num_cells);

    unsigned cell_scalars_uint[] = {0,1,2};
    int cell_scalars_int[] = {0,1,2};
    float cell_scalars_float[] = {0,1,2};
    double cell_scalars_double[] = {0,1,2};
    vtk_write_cell_scalars_uint(
      vtk_file, cell_scalars_uint, num_cells, "cell_scalars_uint");
    vtk_write_cell_scalars_int(
      vtk_file, cell_scalars_int, num_cells, "cell_scalars_int");
    vtk_write_cell_scalars_float(
      vtk_file, cell_scalars_float, num_cells, "cell_scalars_float");
    vtk_write_cell_scalars_double(
      vtk_file, cell_scalars_double, num_cells, "cell_scalars_double");

    unsigned point_scalars_uint[] = {0,1,2,3,4,5,6,7,8,9,10,11};
    int point_scalars_int[] = {0,1,2,3,4,5,6,7,8,9,10,11};
    float point_scalars_float[] = {0,1,2,3,4,5,6,7,8,9,10,11};
    double point_scalars_double[12];
    for (size_t i = 0; i < num_points; ++i)
      point_scalars_double[i] =
        test_func_deg(point_data_lon[i], point_data_lat[i]);
    vtk_write_point_scalars_uint(
      vtk_file, point_scalars_uint, num_points, "point_scalars_uint");
    vtk_write_point_scalars_int(
      vtk_file, point_scalars_int, num_points, "point_scalars_int");
    vtk_write_point_scalars_float(
      vtk_file, point_scalars_float, num_points, "point_scalars_float");
    vtk_write_point_scalars_double(
      vtk_file, point_scalars_double, num_points, "point_scalars_double");

    // closes vtk file
    vtk_close(vtk_file);

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

    generate_cubed_sphere_grid_information(
      n, &num_cells, &num_vertices, &x_vertices, &y_vertices, &z_vertices,
      &vertices_of_cell, &face_id);

    struct {
      char * func_name;
      double (*func)(double, double);
    } test_funcs[] = {
        {.func_name = "test_func", .func = test_func},
        {.func_name = "test_ana_fcos", .func = test_ana_fcos},
        {.func_name = "test_ana_fcossin", .func = test_ana_fcossin},
        {.func_name = "test_one", .func = test_one},
        {.func_name = "test_gulfstream", .func = test_gulfstream},
        {.func_name = "test_harmonic", .func = test_harmonic},
        {.func_name = "test_vortex", .func = test_vortex}
      };

    // initialise and open vtk file
    char const filename[] = "test_vtk_output_2.vtk";
    VTK_FILE * vtk_file = vtk_open(filename, "test2");

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
    vtk_write_point_data(vtk_file, point_data, num_vertices);
    vtk_write_cell_data(
      vtk_file, vertices_of_cell, num_points_per_cell, num_cells);

    double * field_data = xmalloc(num_vertices * sizeof(*field_data));
    for (size_t i = 0; i < sizeof(test_funcs) / sizeof(test_funcs[0]); ++i) {

      for (unsigned j = 0; j < num_vertices; ++j) {

        double lon, lat;
        XYZtoLL(point_data + 3 * j, &lon, &lat);

        field_data[j] = test_funcs[i].func(lon, lat);
      }

      vtk_write_point_scalars_double(
        vtk_file, field_data, num_vertices, test_funcs[i].func_name);
    }

    free(field_data);
    free(num_points_per_cell);
    free(point_data);

    // closes vtk file
    vtk_close(vtk_file);

    unlink(filename);

    free(x_vertices);
    free(y_vertices);
    free(z_vertices);
    free(vertices_of_cell);
    free(face_id);
  }

  return EXIT_SUCCESS;
}
