// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "generate_cubed_sphere.h"
#include "generate_reg2d.h"
#include "geometry.h"
#include "ppm/ppm_xfuncs.h"

#ifdef YAC_NETCDF_ENABLED
#include "io_utils.h"
#include <netcdf.h>
#endif

// This routine is based on Mathlab code provided by Mike Hobson and written by
// Mike Rezny (both from MetOffice)
void yac_generate_cubed_sphere_grid_information(
  unsigned n, unsigned * num_cells, unsigned * num_vertices,
  double ** x_vertices, double ** y_vertices, double ** z_vertices,
  unsigned ** cell_to_vertex, unsigned ** face_id) {

  YAC_ASSERT((n >= 1) && (n <= 40000), "invalid number of linear subdivisions")

  *num_cells = n * n * 6;
  *num_vertices = *num_cells + 2;

  // allocation of output variables
  *x_vertices = xmalloc(*num_vertices * sizeof(**x_vertices));
  *y_vertices = xmalloc(*num_vertices * sizeof(**y_vertices));
  *z_vertices = xmalloc(*num_vertices * sizeof(**z_vertices));

  *cell_to_vertex = xmalloc(*num_cells * 4 * sizeof(**cell_to_vertex));

  *face_id = xmalloc(*num_cells * sizeof(**face_id));

  // allocation of temporary coordinate variables
  unsigned num_edge_coords = n - 1;
  double * cube_edge_x_vertices = xmalloc(12 * num_edge_coords * sizeof(*cube_edge_x_vertices));
  double * cube_edge_y_vertices = xmalloc(12 * num_edge_coords * sizeof(*cube_edge_y_vertices));
  double * cube_edge_z_vertices = xmalloc(12 * num_edge_coords * sizeof(*cube_edge_z_vertices));

  unsigned num_inner_coords = num_edge_coords * num_edge_coords;
  double * cube_inner_x_vertices = xmalloc(num_inner_coords * sizeof(*cube_inner_x_vertices));
  double * cube_inner_y_vertices = xmalloc(num_inner_coords * sizeof(*cube_inner_y_vertices));
  double * cube_inner_z_vertices = xmalloc(num_inner_coords * sizeof(*cube_inner_z_vertices));

  double * temp_x_vertices = xmalloc((n + 1) * (n + 1) * sizeof(*temp_x_vertices));
  double * temp_y_vertices = xmalloc((n + 1) * (n + 1) * sizeof(*temp_y_vertices));
  double * temp_z_vertices = xmalloc((n + 1) * (n + 1) * sizeof(*temp_z_vertices));

  {
    double * theta = xmalloc((n + 1) * sizeof(*theta));

    {
      double d = (M_PI_2 / (double)n);
      for (unsigned i = 0; i < n + 1; ++i)
        theta[i] = tan(- M_PI_4 + d * (double)i);
    }

    for (unsigned i = 0; i < n + 1; ++i) {
      for (unsigned j = 0; j < n + 1; ++j) {

        double scale =
          sqrt(1.0 / (theta[i] * theta[i] + theta[j] * theta[j] + 1.0));

        temp_x_vertices[i * (n + 1) + j] = theta[j] * scale;
        temp_y_vertices[i * (n + 1) + j] = theta[i] * scale;
        temp_z_vertices[i * (n + 1) + j] = - scale;
      }
    }
    free(theta);
  }

  // Store the coordinates for the 4 vertices on face 1
  (*x_vertices)[0] = temp_x_vertices[ 0 * (n + 1) + 0];
  (*y_vertices)[0] = temp_y_vertices[ 0 * (n + 1) + 0];
  (*z_vertices)[0] = temp_z_vertices[ 0 * (n + 1) + 0];
  (*x_vertices)[1] = temp_x_vertices[ n * (n + 1) + 0];
  (*y_vertices)[1] = temp_y_vertices[ n * (n + 1) + 0];
  (*z_vertices)[1] = temp_z_vertices[ n * (n + 1) + 0];
  (*x_vertices)[2] = temp_x_vertices[ 0 * (n + 1) + n];
  (*y_vertices)[2] = temp_y_vertices[ 0 * (n + 1) + n];
  (*z_vertices)[2] = temp_z_vertices[ 0 * (n + 1) + n];
  (*x_vertices)[3] = temp_x_vertices[ n * (n + 1) + n];
  (*y_vertices)[3] = temp_y_vertices[ n * (n + 1) + n];
  (*z_vertices)[3] = temp_z_vertices[ n * (n + 1) + n];
  // Store the coordinates for the 4 vertices on face 2
  (*x_vertices)[4] =  temp_x_vertices[ 0 * (n + 1) + 0];
  (*y_vertices)[4] =  temp_y_vertices[ 0 * (n + 1) + 0];
  (*z_vertices)[4] = -temp_z_vertices[ 0 * (n + 1) + 0];
  (*x_vertices)[5] =  temp_x_vertices[ n * (n + 1) + 0];
  (*y_vertices)[5] =  temp_y_vertices[ n * (n + 1) + 0];
  (*z_vertices)[5] = -temp_z_vertices[ n * (n + 1) + 0];
  (*x_vertices)[6] =  temp_x_vertices[ 0 * (n + 1) + n];
  (*y_vertices)[6] =  temp_y_vertices[ 0 * (n + 1) + n];
  (*z_vertices)[6] = -temp_z_vertices[ 0 * (n + 1) + n];
  (*x_vertices)[7] =  temp_x_vertices[ n * (n + 1) + n];
  (*y_vertices)[7] =  temp_y_vertices[ n * (n + 1) + n];
  (*z_vertices)[7] = -temp_z_vertices[ n * (n + 1) + n];

  // Store the coordinates for the edges
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    cube_edge_x_vertices[ 0 * num_edge_coords + i] = temp_x_vertices[(1 + i) * (n + 1) + 0];
    cube_edge_y_vertices[ 0 * num_edge_coords + i] = temp_y_vertices[(1 + i) * (n + 1) + 0];
    cube_edge_z_vertices[ 0 * num_edge_coords + i] = temp_z_vertices[(1 + i) * (n + 1) + 0];
    cube_edge_x_vertices[ 1 * num_edge_coords + i] = temp_x_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_y_vertices[ 1 * num_edge_coords + i] = temp_y_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_z_vertices[ 1 * num_edge_coords + i] = temp_z_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_x_vertices[ 2 * num_edge_coords + i] = temp_x_vertices[n * (n + 1) + (1 + i)];
    cube_edge_y_vertices[ 2 * num_edge_coords + i] = temp_y_vertices[n * (n + 1) + (1 + i)];
    cube_edge_z_vertices[ 2 * num_edge_coords + i] = temp_z_vertices[n * (n + 1) + (1 + i)];
    cube_edge_x_vertices[ 3 * num_edge_coords + i] = temp_x_vertices[(1 + i) * (n + 1) + n];
    cube_edge_y_vertices[ 3 * num_edge_coords + i] = temp_y_vertices[(1 + i) * (n + 1) + n];
    cube_edge_z_vertices[ 3 * num_edge_coords + i] = temp_z_vertices[(1 + i) * (n + 1) + n];
    cube_edge_x_vertices[ 4 * num_edge_coords + i] =  temp_x_vertices[(1 + i) * (n + 1) + 0];
    cube_edge_y_vertices[ 4 * num_edge_coords + i] =  temp_y_vertices[(1 + i) * (n + 1) + 0];
    cube_edge_z_vertices[ 4 * num_edge_coords + i] = -temp_z_vertices[(1 + i) * (n + 1) + 0];
    cube_edge_x_vertices[ 5 * num_edge_coords + i] =  temp_x_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_y_vertices[ 5 * num_edge_coords + i] =  temp_y_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_z_vertices[ 5 * num_edge_coords + i] = -temp_z_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_x_vertices[ 6 * num_edge_coords + i] =  temp_x_vertices[n * (n + 1) + (1 + i)];
    cube_edge_y_vertices[ 6 * num_edge_coords + i] =  temp_y_vertices[n * (n + 1) + (1 + i)];
    cube_edge_z_vertices[ 6 * num_edge_coords + i] = -temp_z_vertices[n * (n + 1) + (1 + i)];
    cube_edge_x_vertices[ 7 * num_edge_coords + i] =  temp_x_vertices[(1 + i) * (n + 1) + n];
    cube_edge_y_vertices[ 7 * num_edge_coords + i] =  temp_y_vertices[(1 + i) * (n + 1) + n];
    cube_edge_z_vertices[ 7 * num_edge_coords + i] = -temp_z_vertices[(1 + i) * (n + 1) + n];
    cube_edge_x_vertices[ 8 * num_edge_coords + i] =  temp_z_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_y_vertices[ 8 * num_edge_coords + i] =  temp_y_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_z_vertices[ 8 * num_edge_coords + i] =  temp_x_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_x_vertices[ 9 * num_edge_coords + i] =  temp_z_vertices[n * (n + 1) + (1 + i)];
    cube_edge_y_vertices[ 9 * num_edge_coords + i] =  temp_y_vertices[n * (n + 1) + (1 + i)];
    cube_edge_z_vertices[ 9 * num_edge_coords + i] =  temp_x_vertices[n * (n + 1) + (1 + i)];
    cube_edge_x_vertices[10 * num_edge_coords + i] = -temp_z_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_y_vertices[10 * num_edge_coords + i] =  temp_y_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_z_vertices[10 * num_edge_coords + i] =  temp_x_vertices[0 * (n + 1) + (1 + i)];
    cube_edge_x_vertices[11 * num_edge_coords + i] = -temp_z_vertices[n * (n + 1) + (1 + i)];
    cube_edge_y_vertices[11 * num_edge_coords + i] =  temp_y_vertices[n * (n + 1) + (1 + i)];
    cube_edge_z_vertices[11 * num_edge_coords + i] =  temp_x_vertices[n * (n + 1) + (1 + i)];
  }
  
  // Move the 12 edges to the final vertices array
  unsigned Estart = 9 - 1;
  unsigned Eend = Estart + 12 * num_edge_coords;
  for (unsigned i = 0; i < 12 * num_edge_coords; ++i) {
    (*x_vertices)[i + Estart] = cube_edge_x_vertices[i];
    (*y_vertices)[i + Estart] = cube_edge_y_vertices[i];
    (*z_vertices)[i + Estart] = cube_edge_z_vertices[i];
  }

  free(cube_edge_x_vertices);
  free(cube_edge_y_vertices);
  free(cube_edge_z_vertices);

  // store the internal vertices for face 1
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    for (unsigned j = 0; j < num_edge_coords; ++j) {

      cube_inner_x_vertices[i * num_edge_coords + j] = temp_x_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_y_vertices[i * num_edge_coords + j] = temp_y_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_z_vertices[i * num_edge_coords + j] = temp_z_vertices[(1 + i) * (n + 1) + (1 + j)];
    }
  }

  // Move face 1 to final Vertices array
  unsigned Fstart  = Eend;
  for (unsigned i = 0; i < num_inner_coords; ++i) {
    (*x_vertices)[i + Fstart] = cube_inner_x_vertices[i];
    (*y_vertices)[i + Fstart] = cube_inner_y_vertices[i];
    (*z_vertices)[i + Fstart] = cube_inner_z_vertices[i];
  }

  // store the internal vertices for face 2
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    for (unsigned j = 0; j < num_edge_coords; ++j) {

      cube_inner_x_vertices[i * num_edge_coords + j] =  temp_x_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_y_vertices[i * num_edge_coords + j] =  temp_y_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_z_vertices[i * num_edge_coords + j] = -temp_z_vertices[(1 + i) * (n + 1) + (1 + j)];
    }
  }

  // Move face 2 to final Vertices array
  Fstart += num_inner_coords;
  for (unsigned i = 0; i < num_inner_coords; ++i) {
    (*x_vertices)[i + Fstart] = cube_inner_x_vertices[i];
    (*y_vertices)[i + Fstart] = cube_inner_y_vertices[i];
    (*z_vertices)[i + Fstart] = cube_inner_z_vertices[i];
  }

  // store the internal vertices for face 3
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    for (unsigned j = 0; j < num_edge_coords; ++j) {

      cube_inner_x_vertices[i * num_edge_coords + j] = temp_z_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_y_vertices[i * num_edge_coords + j] = temp_y_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_z_vertices[i * num_edge_coords + j] = temp_x_vertices[(1 + i) * (n + 1) + (1 + j)];
    }
  }

  // Move face 3 to final Vertices array
  Fstart += num_inner_coords;
  for (unsigned i = 0; i < num_inner_coords; ++i) {
    (*x_vertices)[i + Fstart] = cube_inner_x_vertices[i];
    (*y_vertices)[i + Fstart] = cube_inner_y_vertices[i];
    (*z_vertices)[i + Fstart] = cube_inner_z_vertices[i];
  }

  // store the internal vertices for face 4
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    for (unsigned j = 0; j < num_edge_coords; ++j) {

      cube_inner_x_vertices[i * num_edge_coords + j] = -temp_z_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_y_vertices[i * num_edge_coords + j] =  temp_y_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_z_vertices[i * num_edge_coords + j] =  temp_x_vertices[(1 + i) * (n + 1) + (1 + j)];
    }
  }

  // Move face 4 to final Vertices array
  Fstart += num_inner_coords;
  for (unsigned i = 0; i < num_inner_coords; ++i) {
    (*x_vertices)[i + Fstart] = cube_inner_x_vertices[i];
    (*y_vertices)[i + Fstart] = cube_inner_y_vertices[i];
    (*z_vertices)[i + Fstart] = cube_inner_z_vertices[i];
  }

  // store the internal vertices for face 5
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    for (unsigned j = 0; j < num_edge_coords; ++j) {

      cube_inner_x_vertices[i * num_edge_coords + j] = temp_x_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_y_vertices[i * num_edge_coords + j] = temp_z_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_z_vertices[i * num_edge_coords + j] = temp_y_vertices[(1 + i) * (n + 1) + (1 + j)];
    }
  }

  // Move face 5 to final Vertices array
  Fstart += num_inner_coords;
  for (unsigned i = 0; i < num_inner_coords; ++i) {
    (*x_vertices)[i + Fstart] = cube_inner_x_vertices[i];
    (*y_vertices)[i + Fstart] = cube_inner_y_vertices[i];
    (*z_vertices)[i + Fstart] = cube_inner_z_vertices[i];
  }

  // store the internal vertices for face 6
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    for (unsigned j = 0; j < num_edge_coords; ++j) {

      cube_inner_x_vertices[i * num_edge_coords + j] =  temp_x_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_y_vertices[i * num_edge_coords + j] = -temp_z_vertices[(1 + i) * (n + 1) + (1 + j)];
      cube_inner_z_vertices[i * num_edge_coords + j] =  temp_y_vertices[(1 + i) * (n + 1) + (1 + j)];
    }
  }

  // Move face 6 to final Vertices array
  Fstart += num_inner_coords;
  for (unsigned i = 0; i < num_inner_coords; ++i) {
    (*x_vertices)[i + Fstart] = cube_inner_x_vertices[i];
    (*y_vertices)[i + Fstart] = cube_inner_y_vertices[i];
    (*z_vertices)[i + Fstart] = cube_inner_z_vertices[i];
  }

  free(temp_x_vertices);
  free(temp_y_vertices);
  free(temp_z_vertices);
  free(cube_inner_x_vertices);
  free(cube_inner_y_vertices);
  free(cube_inner_z_vertices);

  for (unsigned i = 0; i < 6; ++i)
    for (unsigned j = 0; j < n * n; ++j)
      (*face_id)[i * n * n + j] = i + 1;

  unsigned * edge_vertices = xmalloc(num_edge_coords * 12 * sizeof(*edge_vertices));

  Fstart = Estart + 12 * num_edge_coords;

  for (unsigned i = 0; i < num_edge_coords * 12; ++i) edge_vertices[i] = i + Estart;

  unsigned * F = xmalloc((n + 1) * (n + 1) * sizeof(*F));
  unsigned cell_to_vertex_offset = 0;

// Calculate the indices for all the vertices on face 1
//        v2-----e3-----v4
//        |              |
//        e1     f1     e4
//        |              |
//        v1-----e2-----v3

  F[0] = 0; F[n] = 2; F[n * (n + 1) + 0] = 1; F[n * (n + 1) + n] = 3;
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    F[i + 1] =                 edge_vertices[1 * num_edge_coords + i];
    F[(i + 1) * (n + 1)] =     edge_vertices[0 * num_edge_coords + i];
    F[(i + 1) * (n + 1) + n] = edge_vertices[3 * num_edge_coords + i];
    F[n * (n + 1) + i + 1] =   edge_vertices[2 * num_edge_coords + i];
  }

  for (unsigned i = 0; i < n-1; ++i)
    for (unsigned j = 0; j < n-1; ++j)
      F[(i + 1) * (n + 1) + (j + 1)] = Fstart + i * (n - 1) + j;

  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 1)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 1)];
    }
  }

// Calculate the indices for all the vertices on face 2
//        v6-----e7-----v8
//        |              |
//        e5     f2     e8
//        |              |
//        v5-----e6-----v7

  F[0] = 4; F[n] = 6; F[n * (n + 1) + 0] = 5; F[n * (n + 1) + n] = 7;
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    F[i + 1] =                 edge_vertices[5 * num_edge_coords + i];
    F[(i + 1) * (n + 1)] =     edge_vertices[4 * num_edge_coords + i];
    F[(i + 1) * (n + 1) + n] = edge_vertices[7 * num_edge_coords + i];
    F[n * (n + 1) + i + 1] =   edge_vertices[6 * num_edge_coords + i];
  }

  for (unsigned i = 0; i < n-1; ++i)
    for (unsigned j = 0; j < n-1; ++j)
      F[(i + 1) * (n + 1) + (j + 1)] += (n - 1) * (n - 1);

  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 1)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 1)];
    }
  }

// Calculate the indices for all the vertices on face 3
//        v2-----e10----v6
//        |              |
//        e1     f3     e5
//        |              |
//        v1-----e9-----v5

  F[0] = 0; F[n] = 4; F[n * (n + 1) + 0] = 1; F[n * (n + 1) + n] = 5;
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    F[i + 1] =                 edge_vertices[8 * num_edge_coords + i];
    F[(i + 1) * (n + 1)] =     edge_vertices[0 * num_edge_coords + i];
    F[(i + 1) * (n + 1) + n] = edge_vertices[4 * num_edge_coords + i];
    F[n * (n + 1) + i + 1] =   edge_vertices[9 * num_edge_coords + i];
  }

  for (unsigned i = 0; i < n-1; ++i)
    for (unsigned j = 0; j < n-1; ++j)
      F[(i + 1) * (n + 1) + (j + 1)] += (n - 1) * (n - 1);

  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 1)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 1)];
    }
  }

// Calculate the indices for all the vertices on face 4
//        v4-----e12----v8
//        |              |
//        e4     f4     e8
//        |              |
//        v3-----e11----v7

  F[0] = 2; F[n] = 6; F[n * (n + 1) + 0] = 3; F[n * (n + 1) + n] = 7;
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    F[i + 1] =                 edge_vertices[10 * num_edge_coords + i];
    F[(i + 1) * (n + 1)] =     edge_vertices[ 3 * num_edge_coords + i];
    F[(i + 1) * (n + 1) + n] = edge_vertices[ 7 * num_edge_coords + i];
    F[n * (n + 1) + i + 1] =   edge_vertices[11 * num_edge_coords + i];
  }

  for (unsigned i = 0; i < n-1; ++i)
    for (unsigned j = 0; j < n-1; ++j)
      F[(i + 1) * (n + 1) + (j + 1)] += (n - 1) * (n - 1);

  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 1)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 1)];
    }
  }

// Calculate the indices for all the vertices on face 5
//        v5-----e6-----v7
//        |              |
//        e9     f5     e11
//        |              |
//        v1-----e2-----v3

  F[0] = 0; F[n] = 2; F[n * (n + 1) + 0] = 4; F[n * (n + 1) + n] = 6;
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    F[i + 1] =                 edge_vertices[ 1 * num_edge_coords + i];
    F[(i + 1) * (n + 1)] =     edge_vertices[ 8 * num_edge_coords + i];
    F[(i + 1) * (n + 1) + n] = edge_vertices[10 * num_edge_coords + i];
    F[n * (n + 1) + i + 1] =   edge_vertices[ 5 * num_edge_coords + i];
  }

  for (unsigned i = 0; i < n-1; ++i)
    for (unsigned j = 0; j < n-1; ++j)
      F[(i + 1) * (n + 1) + (j + 1)] += (n - 1) * (n - 1);

  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 1)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 1)];
    }
  }

// Calculate the indices for all the vertices on face 6
//        v6-----e7-----v8
//        |              |
//       e10     f6     e12
//        |              |
//        v2-----e3-----v4

  F[0] = 1; F[n] = 3; F[n * (n + 1) + 0] = 5; F[n * (n + 1) + n] = 7;
  for (unsigned i = 0; i < num_edge_coords; ++i) {
    F[i + 1] =                 edge_vertices[ 2 * num_edge_coords + i];
    F[(i + 1) * (n + 1)] =     edge_vertices[ 9 * num_edge_coords + i];
    F[(i + 1) * (n + 1) + n] = edge_vertices[11 * num_edge_coords + i];
    F[n * (n + 1) + i + 1] =   edge_vertices[ 6 * num_edge_coords + i];
  }

  for (unsigned i = 0; i < n-1; ++i)
    for (unsigned j = 0; j < n-1; ++j)
      F[(i + 1) * (n + 1) + (j + 1)] += (n - 1) * (n - 1);

  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 0)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 1) * (n + 1) + (j + 1)];
      (*cell_to_vertex)[cell_to_vertex_offset++] = F[(i + 0) * (n + 1) + (j + 1)];
    }
  }

  free(edge_vertices);
  free(F);
}

struct yac_basic_grid_data yac_generate_cubed_sphere_grid(unsigned n) {

  unsigned num_cells, num_vertices;
  double * x_vertices, * y_vertices, * z_vertices;
  unsigned * cell_to_vertex;
  unsigned * face_id;

  yac_generate_cubed_sphere_grid_information(
    n, &num_cells, &num_vertices, &x_vertices, &y_vertices, &z_vertices,
    &cell_to_vertex, &face_id);

  free(face_id);

  int * num_vertices_per_cell =
    xmalloc(num_cells * sizeof(*num_vertices_per_cell));
  for (unsigned i = 0; i < num_cells; ++i) num_vertices_per_cell[i] = 4;

  struct yac_basic_grid_data grid_data =
    yac_generate_basic_grid_data_unstruct(
      (size_t)num_vertices, (size_t)num_cells, num_vertices_per_cell,
      x_vertices, y_vertices, (int *)cell_to_vertex);
  grid_data.cell_ids =
    xmalloc((size_t)num_cells * sizeof(*(grid_data.cell_ids)));
  grid_data.vertex_ids =
    xmalloc((size_t)num_vertices * sizeof(*(grid_data.vertex_ids)));
  grid_data.edge_ids =
    xmalloc(grid_data.num_edges * sizeof(*(grid_data.edge_ids)));
  for (size_t i = 0; i < (size_t)num_cells; ++i)
    grid_data.cell_ids[i] = (yac_int)i;
  for (size_t i = 0; i < (size_t)num_vertices; ++i) {
    grid_data.vertex_ids[i] = (yac_int)i;
    grid_data.vertex_coordinates[i][0] = x_vertices[i];
    grid_data.vertex_coordinates[i][1] = y_vertices[i];
    grid_data.vertex_coordinates[i][2] = z_vertices[i];
  }
  for (size_t i = 0; i < grid_data.num_edges; ++i)
    grid_data.edge_ids[i] = (yac_int)i;
  free(num_vertices_per_cell);
  free(x_vertices);
  free(y_vertices);
  free(z_vertices);
  free(cell_to_vertex);

  return grid_data;
}

struct yac_basic_grid * yac_generate_cubed_sphere_basic_grid(
  char const * name, size_t n) {

  return
    yac_basic_grid_new(
      name,
      yac_generate_cubed_sphere_grid((unsigned)n));
}

static void decompose_domain_simple(unsigned n, int size, int * cell_owner) {

  unsigned nbr_cells = n * n * 6;

  for (unsigned i = 0; i < nbr_cells; ++i)
    cell_owner[i] = (i * size + size - 1) / nbr_cells;
}

static void decompose_domain_2d(unsigned n, int size, int * cell_owner) {

  // distribute processes among all six sides of the cubed sphere
  int proc_start[6 + 1];
  for (int i = 0; i < 6 + 1; ++i)
    proc_start[i] = (size * i) / 6;

  // for each side of the cube
  for (unsigned i = 0, offset = 0; i < 6; ++i) {

    int num_procs[2];
    int curr_num_procs = proc_start[i+1] - proc_start[i];

    yac_generate_reg2d_decomp((int[2]){(int)n,(int)n}, curr_num_procs, num_procs);

    int local_start_x[num_procs[0] + 1];
    int local_start_y[num_procs[1] + 1];

    for (int j = 0; j < num_procs[0] + 1; ++j)
      local_start_x[j] = (n * j) / num_procs[0];
    for (int j = 0; j < num_procs[1] + 1; ++j)
      local_start_y[j] = (n * j) / num_procs[1];

    for (int p_y = 0; p_y < num_procs[1]; ++p_y) {
      for (int row = local_start_y[p_y]; row < local_start_y[p_y+1]; ++row) {
        for (int p_x = 0; p_x < num_procs[0]; ++p_x) {
          for (int column = local_start_x[p_x]; column < local_start_x[p_x+1];
               ++column) {

            cell_owner[offset++] = proc_start[i] + p_x + num_procs[0] * p_y;
          }
        }
      }
    }
  }
}

static void decompose_domain(unsigned n, int size, int * cell_owner) {

  if (size <= 6) {
    decompose_domain_simple(n, size, cell_owner);
  } else {
    decompose_domain_2d(n, size, cell_owner);
  }
}

void yac_generate_part_cube_grid_information(
  unsigned n, unsigned * nbr_vertices, unsigned * nbr_cells,
  unsigned ** num_vertices_per_cell, unsigned ** cell_to_vertex,
  double ** x_vertices, double ** y_vertices, double ** x_cells,
  double ** y_cells, int ** global_cell_id, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size) {

  double * x_vertices_3d, * y_vertices_3d, * z_vertices_3d;
  unsigned * dummy;

  // generate global grid
  yac_generate_cubed_sphere_grid_information(n, nbr_cells, nbr_vertices,
                                             &x_vertices_3d, &y_vertices_3d,
                                             &z_vertices_3d, cell_to_vertex,
                                             &dummy);
  *num_vertices_per_cell = xmalloc(*nbr_cells * sizeof(**num_vertices_per_cell));
  for (unsigned i = 0; i < *nbr_cells; ++i) (*num_vertices_per_cell)[i] = 4;
  free(dummy);

  double * x_cell_3d = xmalloc(*nbr_cells * sizeof(*x_cell_3d));
  double * y_cell_3d = xmalloc(*nbr_cells * sizeof(*y_cell_3d));
  double * z_cell_3d = xmalloc(*nbr_cells * sizeof(*z_cell_3d));

  for (unsigned i = 0; i < *nbr_cells; ++i) {

    x_cell_3d[i] = y_cell_3d[i] = z_cell_3d[i] = 0;

    for (unsigned j = 0; j < 4; ++j) {

      x_cell_3d[i] += x_vertices_3d[(*cell_to_vertex)[4 * i + j]];
      y_cell_3d[i] += y_vertices_3d[(*cell_to_vertex)[4 * i + j]];
      z_cell_3d[i] += z_vertices_3d[(*cell_to_vertex)[4 * i + j]];
    }

    double scale = 1.0 / sqrt(x_cell_3d[i] * x_cell_3d[i] + 
                              y_cell_3d[i] * y_cell_3d[i] + 
                              z_cell_3d[i] * z_cell_3d[i]);

    x_cell_3d[i] *= scale;
    y_cell_3d[i] *= scale;
    z_cell_3d[i] *= scale;
  }

  int * cell_is_on_rank = xmalloc(*nbr_cells * sizeof(*cell_is_on_rank));

  decompose_domain(n, size, cell_is_on_rank);

  // mask for required vertices and cells
  int * required_vertices = xcalloc(*nbr_vertices, sizeof(*required_vertices));
  int * required_cells = xcalloc(*nbr_cells, sizeof(*required_cells));

  // mark all local cells and their vertices as required
  for (unsigned i = 0, offset = 0; i < *nbr_cells;
       offset += (*num_vertices_per_cell)[i++]) {

    if (cell_is_on_rank[i] != rank) continue;

    cell_is_on_rank[i] = -1;

    required_cells[i] = 1;
    for (unsigned j = 0; j < (*num_vertices_per_cell)[i]; ++j)
      required_vertices[(*cell_to_vertex)[offset+j]] = 1;
  }

  // mark all halo cells as required and generate cell_is_on_rank
  for (unsigned i = 0, offset = 0; i < *nbr_cells;
       offset += (*num_vertices_per_cell)[i++]) {

    if (!required_cells[i]) {

      for (unsigned j = 0; j < (*num_vertices_per_cell)[i]; ++j) {

        if (required_vertices[(*cell_to_vertex)[offset+j]]) {

          required_cells[i] = 1;
          break;
        }
      }

    }
  }

  // mask for halo vertices
  int * vertex_is_on_rank = xmalloc(*nbr_vertices * sizeof(*vertex_is_on_rank));

  // mark all halo vertices as required and generate vertex_is_on_rank
  for (unsigned i = 0; i < *nbr_vertices; ++i)
    if (required_vertices[i])
      vertex_is_on_rank[i] = -1;
  for (unsigned i = 0, offset = 0; i < *nbr_cells;
       offset += (*num_vertices_per_cell)[i++]) {

    if (required_cells[i] && cell_is_on_rank[i] != -1) {

      for (unsigned j = 0; j < (*num_vertices_per_cell)[i]; ++j) {

        if (!required_vertices[(*cell_to_vertex)[offset+j]]) {

          required_vertices[(*cell_to_vertex)[offset+j]] = 1;
          vertex_is_on_rank[(*cell_to_vertex)[offset+j]] = cell_is_on_rank[i];
        }
      }
    }
  }

  // count the number of cells and vertices
  int part_num_vertices = 0;
  int part_num_cells = 0;
  for (unsigned i = 0; i < *nbr_vertices; ++i)
    if (required_vertices[i])
      part_num_vertices++;
  for (unsigned i = 0; i < *nbr_cells; ++i)
    if(required_cells[i])
      part_num_cells++;

  *global_cell_id = xmalloc(part_num_cells * sizeof(**global_cell_id));
  *cell_core_mask = xmalloc(part_num_cells * sizeof(**cell_core_mask));
  *global_corner_id = xmalloc(part_num_vertices * sizeof(**global_corner_id));
  *corner_core_mask = xmalloc(part_num_vertices * sizeof(**corner_core_mask));

  *x_cells = xmalloc(part_num_cells * sizeof(**x_cells));
  *y_cells = xmalloc(part_num_cells * sizeof(**y_cells));
  *x_vertices = xmalloc(part_num_vertices * sizeof(**x_vertices));
  *y_vertices = xmalloc(part_num_vertices * sizeof(**y_vertices));

  // generate final vertex data
  part_num_vertices = 0;
  int * global_to_local_vertex = xmalloc(*nbr_vertices * sizeof(*global_to_local_vertex));
  for (unsigned i = 0; i < *nbr_vertices; ++i) {

    if (required_vertices[i]) {

      (*global_corner_id)[part_num_vertices] = i;
      (*corner_core_mask)[part_num_vertices] = vertex_is_on_rank[i] == -1;
      double p[3] = {x_vertices_3d[i], y_vertices_3d[i], z_vertices_3d[i]};
      XYZtoLL(p, (*x_vertices) + part_num_vertices, (*y_vertices) + part_num_vertices);
      global_to_local_vertex[i] = part_num_vertices;
      part_num_vertices++;
    }
  }

  free(vertex_is_on_rank);
  *nbr_vertices = part_num_vertices;
  free(required_vertices);

  // generate final cell data
  int num_cell_vertex_dependencies = 0;
  part_num_cells = 0;
  for (unsigned i = 0, offset = 0; i < *nbr_cells;
       offset += (*num_vertices_per_cell)[i++]) {

    if (required_cells[i]) {

      (*global_cell_id)[part_num_cells] = i;
      (*cell_core_mask)[part_num_cells] = (cell_is_on_rank[i] == -1);
      double middle_point[3] = {x_cell_3d[i],y_cell_3d[i],z_cell_3d[i]};

      for (unsigned j = 0; j < (*num_vertices_per_cell)[i]; ++j) {
        (*cell_to_vertex)[num_cell_vertex_dependencies+j] =
          global_to_local_vertex[(*cell_to_vertex)[offset+j]];
      }

      XYZtoLL(middle_point, (*x_cells)+part_num_cells, (*y_cells)+part_num_cells);

      num_cell_vertex_dependencies += (*num_vertices_per_cell)[i];

      (*num_vertices_per_cell)[part_num_cells] = (*num_vertices_per_cell)[i];

      part_num_cells++;
    }
  }

  free(x_cell_3d);
  free(y_cell_3d);
  free(z_cell_3d);
  free(x_vertices_3d);
  free(y_vertices_3d);
  free(z_vertices_3d);

  *num_vertices_per_cell = xrealloc(*num_vertices_per_cell, part_num_cells *
                                    sizeof(**num_vertices_per_cell));
  *cell_to_vertex = xrealloc(*cell_to_vertex, num_cell_vertex_dependencies *
                             sizeof(**cell_to_vertex));
  free(cell_is_on_rank);
  *nbr_cells = part_num_cells;
  free(required_cells);
  free(global_to_local_vertex);
}

#ifdef YAC_NETCDF_ENABLED
void yac_write_cubed_sphere_grid(unsigned n, char const * filename) {

  // generate basic grid information
  unsigned num_cells, num_vertices;
  double * x_vertices, * y_vertices, * z_vertices;
  unsigned * vertices_of_cell, * dummy;
  yac_generate_cubed_sphere_grid_information(
    n, &num_cells, &num_vertices, &x_vertices, &y_vertices, &z_vertices,
    &vertices_of_cell, &dummy);
  free(dummy);

  // create file and define dimensions and variables
  size_t nv = 4, ne = 4;
  int ncid;
  int dim_cell_id, dim_vertex_id, dim_nv_id, dim_ne_id;
  int var_vlon_id, var_vlat_id, var_clon_id, var_clat_id, var_mask_id,
      var_v2c_id, var_c2v_id;
  yac_nc_create(filename, NC_CLOBBER, &ncid);
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "cell", (size_t)num_cells, &dim_cell_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "vertex", (size_t)num_vertices, &dim_vertex_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "nv", 4, &dim_nv_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "ne", 4, &dim_ne_id));
  YAC_HANDLE_ERROR(
    nc_def_var(ncid, "vlon", NC_DOUBLE, 1, &dim_vertex_id, &var_vlon_id));
  YAC_HANDLE_ERROR(
    nc_def_var(ncid, "vlat", NC_DOUBLE, 1, &dim_vertex_id, &var_vlat_id));
  YAC_HANDLE_ERROR(
    nc_def_var(ncid, "clon", NC_DOUBLE, 1, &dim_cell_id, &var_clon_id));
  YAC_HANDLE_ERROR(
    nc_def_var(ncid, "clat", NC_DOUBLE, 1, &dim_cell_id, &var_clat_id));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "cell_sea_land_mask", NC_INT, 1, &dim_cell_id, &var_mask_id));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "vertex_of_cell", NC_INT, 2, (int[]){dim_nv_id, dim_cell_id},
      &var_c2v_id));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "cells_of_vertex", NC_INT, 2, (int[]){dim_ne_id, dim_vertex_id},
      &var_v2c_id));
  YAC_HANDLE_ERROR(nc_enddef(ncid));

    // generate and write grid data

  double * lon_buffer =
    xmalloc(MAX(num_cells, num_vertices) * sizeof(*lon_buffer));
  double * lat_buffer =
    xmalloc(MAX(num_cells, num_vertices) * sizeof(*lat_buffer));

  for (unsigned i = 0; i < num_vertices; ++i) {
    double coord_3d[3] = {x_vertices[i], y_vertices[i], z_vertices[i]};
    XYZtoLL(coord_3d, lon_buffer + i, lat_buffer + i);
  }
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_vlon_id, lon_buffer));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_vlat_id, lat_buffer));

  for (unsigned i = 0; i < num_cells; ++i) {
    double coord_3d[3] = {0.0, 0.0, 0.0};
    for (unsigned j = 0; j < 4; ++j) {
      unsigned vertex_idx = vertices_of_cell[4 * i + j];
      coord_3d[0] += x_vertices[vertex_idx];
      coord_3d[1] += y_vertices[vertex_idx];
      coord_3d[2] += z_vertices[vertex_idx];
    }
    normalise_vector(coord_3d);
    XYZtoLL(coord_3d, lon_buffer + i, lat_buffer + i);
  }
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_clon_id, lon_buffer));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_clat_id, lat_buffer));

  free(lat_buffer);
  free(lon_buffer);
  free(z_vertices);
  free(y_vertices);
  free(x_vertices);

  int * int_buffer =
    xmalloc(MAX(nv, ne) * MAX(num_cells, num_vertices) * sizeof(*int_buffer));

  for (unsigned i = 0; i < num_cells; ++i) int_buffer[i] = 1;
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_mask_id, int_buffer));

  for (unsigned i = 0; i < ne * num_vertices; ++i)
    int_buffer[i] = INT_MAX;
  for (unsigned i = 0; i < num_cells; ++i) {
    int cell_id = (int)(i + 1);
    for (unsigned j = 0; j < ne; ++j) {
      unsigned vertex_idx = vertices_of_cell[ne * i + j];
      for (unsigned k = 0; k < ne; ++k) {
        if (int_buffer[k * num_vertices + vertex_idx] == INT_MAX) {
          int_buffer[k * num_vertices + vertex_idx] = cell_id;
          break;
        }
      }
    }
  }
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_v2c_id, int_buffer));

  for (unsigned i = 0; i < num_cells; ++i)
    for (unsigned j = 0; j < nv; ++j)
      int_buffer[j * num_cells + i] =
        vertices_of_cell[nv * i + j] + 1;
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_c2v_id, int_buffer));

  free(int_buffer);
  free(vertices_of_cell);

  YAC_HANDLE_ERROR(nc_close(ncid));
}
#else
void yac_write_cubed_sphere_grid(unsigned n, char const * filename) {

   UNUSED(n);
   UNUSED(filename);
   die(
     "ERROR(yac_write_cubed_sphere_grid): "
     "YAC is built without the NetCDF support");
}
#endif // YAC_NETCDF_ENABLED
