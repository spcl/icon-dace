/**
 * @file clipping.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
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

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "geometry.h"
#include "clipping.h"
#include "area.h"
#include "ensure_array_size.h"
#include "utils.h"

//#define YAC_VERBOSE_CLIPPING

static double const tol = 1.0e-12;

struct point_list_element {

  double vec_coords[3];
  struct yac_circle * edge_circle;
  int to_be_removed;
  struct point_list_element * next;
};

struct point_list {

  struct point_list_element * start;
  struct point_list_element * free_elements;
};

/* internal helper routines for working with linked lists of points */

static void init_point_list(struct point_list * list);

static void reset_point_list(struct point_list * list);

static size_t generate_point_list(
  struct point_list * list, struct grid_cell cell, int cell_ordering,
  struct yac_circle * circle_buffer);

static struct point_list_element *
get_free_point_list_element(struct point_list * list);

static size_t remove_points(struct point_list * list);
static size_t remove_zero_length_edges(struct point_list * list);

static void free_point_list(struct point_list * list);

static int get_cell_points_ordering(struct grid_cell cell);

static void generate_cell(struct point_list * list,
                          struct grid_cell * cell);

static enum yac_cell_type get_cell_type(struct grid_cell target_cell);

/* ------------------------- */

static struct grid_cell * overlap_cell_buffer = NULL;
static size_t overlap_cell_buffer_size = 0;

static inline struct grid_cell * get_overlap_cell_buffer(size_t N) {

  // ensure that there are enough buffer cells

  if (overlap_cell_buffer_size < N) {

    size_t old_overlap_cell_buffer_size = overlap_cell_buffer_size;

    ENSURE_ARRAY_SIZE(overlap_cell_buffer, overlap_cell_buffer_size, N);

    for (; old_overlap_cell_buffer_size < overlap_cell_buffer_size;
         ++old_overlap_cell_buffer_size)
      yac_init_grid_cell(overlap_cell_buffer + old_overlap_cell_buffer_size);
  }

  return overlap_cell_buffer;
}

/* ------------------------- */

static double get_edge_direction(
  double * ref_corner, double * corner_a, double * corner_b) {

  double edge_norm[3];
  crossproduct_kahan(corner_a, corner_b, edge_norm);
  normalise_vector(edge_norm);

  // sine of the angle between the edge and the reference corner
  double angle =
    edge_norm[0] * ref_corner[0] +
    edge_norm[1] * ref_corner[1] +
    edge_norm[2] * ref_corner[2];

  // if the reference corner is directly on the edge
  // (for small angles sin(x)==x)
  if (fabs(angle) < yac_angle_tol) return 0.0;

  return copysign(1.0, angle);
}

void yac_compute_overlap_info (size_t N,
                               struct grid_cell * source_cell,
                               struct grid_cell target_cell,
                               double * overlap_areas,
                               double (*overlap_barycenters)[3]) {

  YAC_ASSERT(
    target_cell.num_corners > 2,
    "ERROR(yac_compute_overlap_info): "
    "target cell has too few corners")

  struct grid_cell * overlap_buffer = get_overlap_cell_buffer(N);
  enum yac_cell_type target_cell_type;

  // initialise barycenter coordinates if available
  if (overlap_barycenters != NULL)
    for (size_t i = 0; i < N; ++i)
      for (int j = 0; j < 3; ++j)
        overlap_barycenters[i][j] = 0.0;

  // if the target is not a triangle
  if ( target_cell.num_corners > 3 )
    target_cell_type = get_cell_type (target_cell);

  // special case: target is a triangle or a lon-lat cell -> no triangulation
  // of the target is required
  if ( target_cell.num_corners < 4 || target_cell_type == LON_LAT_CELL ) {
    yac_cell_clipping ( N, source_cell, target_cell, overlap_buffer);
    for (size_t i = 0; i < N; ++i) {
      if (overlap_buffer[i].num_corners > 1) {
        if (overlap_barycenters == NULL)
          overlap_areas[i] = yac_huiliers_area (overlap_buffer[i]);
        else {
          overlap_areas[i] =
            yac_huiliers_area_info(
              overlap_buffer[i], overlap_barycenters[i], 1.0);
          if ((overlap_barycenters[i][0] != 0.0) ||
              (overlap_barycenters[i][1] != 0.0) ||
              (overlap_barycenters[i][2] != 0.0))
            normalise_vector(overlap_barycenters[i]);
          if (overlap_areas[i] < 0.0) {
            overlap_areas[i] = -overlap_areas[i];
            overlap_barycenters[i][0] = -overlap_barycenters[i][0];
            overlap_barycenters[i][1] = -overlap_barycenters[i][1];
            overlap_barycenters[i][2] = -overlap_barycenters[i][2];
          }
        }
      } else {
        overlap_areas[i] = 0.0;
      }
    }
    return;
  }

  // the triangulation algorithm only works for cells that only have great
  // circle edges
  // in order to also support other edge types, additional work would be
  // required
  YAC_ASSERT(
    target_cell_type == GREAT_CIRCLE_CELL,
    "ERROR(yac_compute_overlap_info): invalid target cell type")

  // data structure to hold the triangles of the target cell
  struct grid_cell target_partial_cell =
    {.coordinates_xyz = (double[3][3]){{-1}},
     .edge_type       = (enum yac_edge_type[3]){GREAT_CIRCLE_EDGE,
                                                GREAT_CIRCLE_EDGE,
                                                GREAT_CIRCLE_EDGE},
     .num_corners     = 3};
  // common node point to all partial target cells
  double * base_corner = target_cell.coordinates_xyz[0];
  target_partial_cell.coordinates_xyz[0][0] = base_corner[0];
  target_partial_cell.coordinates_xyz[0][1] = base_corner[1];
  target_partial_cell.coordinates_xyz[0][2] = base_corner[2];

  // initialise overlap areas
  for ( size_t n = 0; n < N; n++) overlap_areas[n] = 0.0;

  // for all triangles of the target cell
  // (triangles a formed by first corner of the target cells and each edge of
  //  the cell; the first and last edge of the cell already, contain the
  //  first corner, therefore we can skip them)
  for ( size_t corner_idx = 1;
        corner_idx < target_cell.num_corners - 1; ++corner_idx ) {

    double * corner_a = target_cell.coordinates_xyz[corner_idx];
    double * corner_b = target_cell.coordinates_xyz[corner_idx+1];

    // if the current edge has a length of zero
    if (points_are_identically(corner_a, corner_b)) continue;

    double edge_direction =
      get_edge_direction(base_corner, corner_a, corner_b);

    target_partial_cell.coordinates_xyz[1][0] = corner_a[0];
    target_partial_cell.coordinates_xyz[1][1] = corner_a[1];
    target_partial_cell.coordinates_xyz[1][2] = corner_a[2];
    target_partial_cell.coordinates_xyz[2][0] = corner_b[0];
    target_partial_cell.coordinates_xyz[2][1] = corner_b[1];
    target_partial_cell.coordinates_xyz[2][2] = corner_b[2];

    // clip the current target cell triangle with all source cells
    yac_cell_clipping(N, source_cell, target_partial_cell, overlap_buffer);

    // for all source cells
    for (size_t n = 0; n < N; n++) {

      if (overlap_buffer[n].num_corners == 0) continue;

      if (overlap_barycenters == NULL)
        overlap_areas[n] +=
          yac_huiliers_area(overlap_buffer[n]) * edge_direction;
      else
        overlap_areas[n] +=
          yac_huiliers_area_info(
            overlap_buffer[n], overlap_barycenters[n], edge_direction);
    }
  }

  for (size_t n = 0; n < N; n++) {

    if (overlap_areas[n] < 0.0) {

      overlap_areas[n] = -overlap_areas[n];

      if (overlap_barycenters != NULL) {
        overlap_barycenters[n][0] = -overlap_barycenters[n][0];
        overlap_barycenters[n][1] = -overlap_barycenters[n][1];
        overlap_barycenters[n][2] = -overlap_barycenters[n][2];
      }
    }
  }

  if (overlap_barycenters != NULL)
    for (size_t n = 0; n < N; n++)
      if ((overlap_areas[n] > 0.0) &&
          ((overlap_barycenters[n][0] != 0.0) ||
           (overlap_barycenters[n][1] != 0.0) ||
           (overlap_barycenters[n][2] != 0.0)))
        normalise_vector(overlap_barycenters[n]);

#ifdef YAC_VERBOSE_CLIPPING
  for (size_t n = 0; n < N; n++)
    printf("overlap area %zu: %lf \n", n, overlap_areas[n]);
#endif
}

static inline double dotproduct(double a[], double b[]) {

  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* ------------------------- */

void yac_compute_overlap_areas (size_t N,
                                struct grid_cell * source_cell,
                                struct grid_cell target_cell,
                                double * partial_areas) {

  yac_compute_overlap_info(
    N, source_cell, target_cell, partial_areas, NULL);
}

/* ------------------------- */

static enum yac_cell_type get_cell_type(struct grid_cell target_cell) {

  YAC_ASSERT(
    target_cell.num_corners > 0, "ERROR(get_cell_type): CELL without edge")

  enum yac_cell_type cell_type = MIXED_CELL;

  // if the cell is a typical lon-lat cell
  if ((target_cell.num_corners == 4) &&
      ((target_cell.edge_type[0] == LAT_CIRCLE_EDGE &&
        target_cell.edge_type[1] == LON_CIRCLE_EDGE &&
        target_cell.edge_type[2] == LAT_CIRCLE_EDGE &&
        target_cell.edge_type[3] == LON_CIRCLE_EDGE) ||
       (target_cell.edge_type[0] == LON_CIRCLE_EDGE &&
        target_cell.edge_type[1] == LAT_CIRCLE_EDGE &&
        target_cell.edge_type[2] == LON_CIRCLE_EDGE &&
        target_cell.edge_type[3] == LAT_CIRCLE_EDGE))) {

    cell_type = LON_LAT_CELL;

  } else {

    size_t count_lat_edges = 0, count_great_circle_edges = 0;

    // count the number of edges for each type
    // (lon edges are counted as gc ones)
    for (size_t i = 0; i < target_cell.num_corners; ++i)
      if (target_cell.edge_type[i] == LON_CIRCLE_EDGE ||
          target_cell.edge_type[i] == GREAT_CIRCLE_EDGE)
        count_great_circle_edges++;
      else
        count_lat_edges++;

    // if the cell is a lon lat cell with one lat edge having a length of zero
    // due to being at a pole
    if ((count_lat_edges == 1) && (count_great_circle_edges == 2)) {

      size_t i;
      for (i = 0; i < 3; ++i) if (target_cell.edge_type[i] == LAT_CIRCLE_EDGE) break;
      size_t pol_index = (3 + i - 1) % 3;
      if (fabs(fabs(target_cell.coordinates_xyz[pol_index][2])-1.0) <
          yac_angle_low_tol) cell_type = LON_LAT_CELL;

    // if the cell only consists of great circle edges
    } else if (count_lat_edges == 0) cell_type = GREAT_CIRCLE_CELL;

    // if the cell only consists of lat circle edges
    else if (count_great_circle_edges == 0) cell_type = LAT_CELL;
  }

  YAC_ASSERT(
    cell_type != MIXED_CELL,
    "invalid cell type (cell contains edges consisting "
    "of great circles and circles of latitude)")

  return cell_type;
}

int yac_circle_compare(void const * a, void const * b) {

  struct yac_circle const * circle_a = *(struct yac_circle const **)a;
  struct yac_circle const * circle_b = *(struct yac_circle const **)b;

  YAC_ASSERT(
    ((circle_a->type == GREAT_CIRCLE) ||
     (circle_a->type == LON_CIRCLE) ||
     (circle_a->type == LAT_CIRCLE)/* ||
     (circle_a->type == POINT)*/) &&
    ((circle_b->type == GREAT_CIRCLE) ||
     (circle_b->type == LON_CIRCLE) ||
     (circle_b->type == LAT_CIRCLE)/* ||
     (circle_b->type == POINT)*/),
    "ERROR(yac_circle_compare): unsupported circle type")

  // put lat circles to the end
  int ret = (circle_a->type == LAT_CIRCLE) - (circle_b->type == LAT_CIRCLE);
  if (!ret) ret = (int)(circle_a->type) - (int)(circle_b->type);
  if (!ret) {
    switch (circle_a->type) {
      default:
      case (GREAT_CIRCLE):
        for (int i = 0; !ret && (i < 3); ++i) {
          ret =
            (circle_a->data.gc.norm_vector[i] >
             circle_b->data.gc.norm_vector[i]) -
            (circle_a->data.gc.norm_vector[i] <
             circle_b->data.gc.norm_vector[i]);
        }
        break;
      case (LON_CIRCLE):
        for (int i = 0; !ret && (i < 3); ++i) {
          ret =
            (circle_a->data.lon.norm_vector[i] >
             circle_b->data.lon.norm_vector[i]) -
            (circle_a->data.lon.norm_vector[i] <
             circle_b->data.lon.norm_vector[i]);
        }
        break;
      case (LAT_CIRCLE):
        ret = circle_a->data.lat.north_is_out - circle_b->data.lat.north_is_out;
        if (!ret)
          ret = (circle_a->data.lat.z > circle_b->data.lat.z) -
                (circle_a->data.lat.z < circle_b->data.lat.z);
        break;
      // case (POINT):
        // for (int i = 0; !ret && (i < 3); ++i)
          // ret = (circle_a->data.p.vec[i] > circle_b->data.p.vec[i]) -
                // (circle_a->data.p.vec[i] < circle_b->data.p.vec[i]);
        // break;
    }
  }

  return ret;
}


static void compute_norm_vector_kahan(
  double const a[3], double const b[3], double norm_vector[3]) {

  crossproduct_kahan(a, b, norm_vector);

  double scale = 1.0 / sqrt(norm_vector[0] * norm_vector[0] +
                            norm_vector[1] * norm_vector[1] +
                            norm_vector[2] * norm_vector[2]);
  norm_vector[0] *= scale;
  norm_vector[1] *= scale;
  norm_vector[2] *= scale;
}

void yac_circle_generate(
  double const * a, double const * b, enum yac_edge_type type,
  int edge_ordering, struct yac_circle * circle) {

  if (points_are_identically(a, b)) {
    circle->type = POINT;
    circle->data.p.vec[0] = a[0];
    circle->data.p.vec[1] = a[1];
    circle->data.p.vec[2] = a[2];
  } else {

    // switch edges to ensure that "inside/outside" of the circle is computed
    // correctly
    if (edge_ordering < 0) {
      double const * temp = a;
      a = b;
      b = temp;
    }

    YAC_ASSERT(
      (type == GREAT_CIRCLE_EDGE) ||
      (type == LON_CIRCLE_EDGE) ||
      (type == LAT_CIRCLE_EDGE),
      "ERROR(yac_circle_generate): invalid edge type")

    switch (type) {
      default:
      case(GREAT_CIRCLE_EDGE):
        circle->type = GREAT_CIRCLE;
        compute_norm_vector_kahan(a, b, circle->data.gc.norm_vector);
        break;
      case(LON_CIRCLE_EDGE):
        circle->type = LON_CIRCLE;
        compute_norm_vector_kahan(a, b, circle->data.lon.norm_vector);
        break;
      case(LAT_CIRCLE_EDGE):
        circle->type = LAT_CIRCLE;
        circle->data.lat.north_is_out = (a[0] * b[1] - a[1] * b[0]) < 0.0;
        if (circle->data.lat.north_is_out == 1337) exit(EXIT_FAILURE);
        circle->data.lat.z = a[2];
        break;
    }
  }
}

static inline struct yac_circle
  generate_lat_circle(double z, int north_is_out) {

  struct yac_circle circle = {
    .type = LAT_CIRCLE,
    .data.lat.north_is_out = north_is_out,
    .data.lat.z = z
  };
  return circle;
}

int yac_circle_point_is_inside(
  double const point[3], struct yac_circle * circle) {

  YAC_ASSERT(
    (circle->type == GREAT_CIRCLE) ||
    (circle->type == LON_CIRCLE) ||
    (circle->type == LAT_CIRCLE)/* ||
    (circle->type == POINT)*/,
    "ERROR(yac_circle_point_is_inside): unsupported circle type")

  switch (circle->type) {
    default:
    case (GREAT_CIRCLE): {
      double dot = point[0] * circle->data.gc.norm_vector[0] +
                   point[1] * circle->data.gc.norm_vector[1] +
                   point[2] * circle->data.gc.norm_vector[2];
      if      (dot < (- yac_angle_tol * 1e-3)) return 0;
      else if (dot > (+ yac_angle_tol * 1e-3)) return 1;
      else return 2;
    }
    case (LON_CIRCLE): {
      double dot = point[0] * circle->data.lon.norm_vector[0] +
                   point[1] * circle->data.lon.norm_vector[1] +
                   point[2] * circle->data.lon.norm_vector[2];
      if      (dot < (- yac_angle_tol * 1e-3)) return 0;
      else if (dot > (+ yac_angle_tol * 1e-3)) return 1;
      else return 2;
    }
    case (LAT_CIRCLE): {
      double diff_z = circle->data.lat.z - point[2];
      if (fabs(diff_z) < yac_angle_tol * 1e-3) return 2;
      else return (diff_z < 0.0) ^ circle->data.lat.north_is_out;
    }
    // case (POINT):
      // return points_are_identically(point, circle->data.p.vec)?2:0;
  }
}

int yac_circle_compare_distances(
  double const a[3], double const b[3], struct yac_circle * circle) {

  YAC_ASSERT(
    (circle->type == GREAT_CIRCLE) ||
    (circle->type == LON_CIRCLE) ||
    (circle->type == LAT_CIRCLE) ||
    (circle->type == POINT),
    "ERROR(yac_circle_compare_distances): invalid circle type")

  double dist_a, dist_b;
  switch (circle->type) {
    default:
    case(GREAT_CIRCLE): {
      double norm_vector[3] = {circle->data.gc.norm_vector[0],
                               circle->data.gc.norm_vector[1],
                               circle->data.gc.norm_vector[2]};
      dist_a = fabs(a[0] * norm_vector[0] +
                    a[1] * norm_vector[1] +
                    a[2] * norm_vector[2]);
      dist_b = fabs(b[0] * norm_vector[0] +
                    b[1] * norm_vector[1] +
                    b[2] * norm_vector[2]);
      break;
    }
    case(LON_CIRCLE): {
      double norm_vector[3] = {circle->data.lon.norm_vector[0],
                               circle->data.lon.norm_vector[1],
                               circle->data.lon.norm_vector[2]};
      dist_a = fabs(a[0] * norm_vector[0] +
                    a[1] * norm_vector[1] +
                    a[2] * norm_vector[2]);
      dist_b = fabs(b[0] * norm_vector[0] +
                    b[1] * norm_vector[1] +
                    b[2] * norm_vector[2]);
      break;
    }
    case(LAT_CIRCLE): {
      double circle_lat = acos(circle->data.lat.z);
      dist_a = fabs(circle_lat - acos(a[2]));
      dist_b = fabs(circle_lat - acos(b[2]));
      break;
    }
    case(POINT): {
      dist_a = get_vector_angle(circle->data.p.vec, a);
      dist_b = get_vector_angle(circle->data.p.vec, b);
      break;
    }
  }

  return (dist_a > dist_b + yac_angle_tol) -
         (dist_a + yac_angle_tol < dist_b);
}

int yac_circle_contains_north_pole(struct yac_circle * circle) {

  YAC_ASSERT(
    (circle->type == GREAT_CIRCLE) ||
    (circle->type == LAT_CIRCLE) ||
    (circle->type == LON_CIRCLE),
    "ERROR(yac_circle_contains_north_pole): circle type has to be either "
    "GREAT_CIRCLE or LAT_CIRCLE")

  switch (circle->type) {
    default:
    case (GREAT_CIRCLE):
      return circle->data.gc.norm_vector[2] > 0.0;
    case (LAT_CIRCLE):
      return !circle->data.lat.north_is_out;
    case (LON_CIRCLE):
      return 0;
  }
}

/**
 * cell clipping using Sutherland-Hodgman algorithm;
 */
static void circle_clipping(
  struct point_list * cell, size_t num_cell_edges,
  struct yac_circle ** clipping_circles, size_t num_clipping_circles) {

  // to avoid some problems that can occur close the the pole, we process the
  // target lat-circle edges at the end
  qsort(clipping_circles, num_clipping_circles, sizeof(*clipping_circles),
        yac_circle_compare);

  // for all clipping circles
  for (size_t i = 0; (i < num_clipping_circles) && (num_cell_edges > 1); ++i) {

    struct point_list_element * cell_edge_start = cell->start;
    struct point_list_element * cell_edge_end = cell_edge_start->next;

    struct yac_circle * clipping_circle = clipping_circles[i];

    int start_is_inside, first_start_is_inside;

    start_is_inside =
      yac_circle_point_is_inside(
        cell_edge_start->vec_coords, clipping_circle);
    first_start_is_inside = start_is_inside;

    // for all edges of the cell
    for (size_t cell_edge_idx = 0; cell_edge_idx < num_cell_edges;
         ++cell_edge_idx) {

      int end_is_inside =
        (cell_edge_idx != num_cell_edges - 1)?
          (yac_circle_point_is_inside(
             cell_edge_end->vec_coords, clipping_circle)):first_start_is_inside;

      double * cell_edge_coords[2] =
        {cell_edge_start->vec_coords, cell_edge_end->vec_coords};
      struct yac_circle * cell_edge_circle = cell_edge_start->edge_circle;

      double intersection[2][3];
      int num_edge_intersections = -1;

      // One intersection between the clipping circle and the cell edge is
      // possible, if either the start or the end vertex of the current cell
      // edge is inside and the respective other one is outside.
      // Two intersections are possible, if either the clipping circle or the
      // cell edge is a circle of latitude, while this other is not.
      // Additionally, this is not possible if one of the two cell edge
      // vertices is inside while the other is not or if both cell edge
      // vertices are on the plane of the clipping edge.
      int one_intersect_expected = (end_is_inside + start_is_inside == 1);
      int two_intersect_possible =
        ((cell_edge_circle->type == LAT_CIRCLE) ^
         (clipping_circle->type == LAT_CIRCLE)) &&
        (end_is_inside + start_is_inside != 1) &&
        (end_is_inside + start_is_inside != 4);

      // if we need to check for intersections
      if (one_intersect_expected || two_intersect_possible) {

        // get intersection points between the clipping circle and the circle
        // of the cell edge
        int num_circle_intersections =
          yac_circle_intersect(
            *clipping_circle, *cell_edge_circle,
            intersection[0], intersection[1]);

        YAC_ASSERT(
          (num_circle_intersections >= -1) && (num_circle_intersections <= 2),
          "ERROR(circle_clipping): invalid number of intersections")

        // determine the intersections of the clipping circle with the cell
        // edge based on the circle intersections
        switch (num_circle_intersections) {

          // special case:
          // both circles are on the same plane
          case (-1): {

            // MoHa: this part of the code should not be reachable...but I
            // leave it just in case...
            start_is_inside = 2;
            end_is_inside = 2;
            num_edge_intersections = 0;

            one_intersect_expected = 0;

            // if this is the first iteration
            if (cell_edge_idx == 0) first_start_is_inside = 2;
            break;
          }
          // special case:
          // no intersections between the two circles
          // (can occure if one of the two circles is a latitude circle while
          //  the other is a great circle)
          case (0): {
            num_edge_intersections = 0;
            break;
          }
          // standart case:
          // two intersections between the two circles
          // (one intersection between two circles is a special case, but is
          //  being handled here anyway)
          default:
          case (1):
          case (2): {

            // check the relation between the intersection points and the
            // cell edge
            int is_on_edge[2] = {
              yac_point_on_edge(
                intersection[0], cell_edge_coords[0], cell_edge_coords[1],
                cell_edge_circle->type),
              (num_circle_intersections == 2)?
                yac_point_on_edge(
                  intersection[1], cell_edge_coords[0], cell_edge_coords[1],
                  cell_edge_circle->type):0};

            // if both intersection points are on the edge
            if (is_on_edge[0] && is_on_edge[1]) {

              // if only one intersection was expected
              YAC_ASSERT_F(
                !one_intersect_expected,
                "ERROR: two intersections found, even "
                "though no circle of latitude involed and "
                "both edge vertices on different sides of "
                "the cutting plane.\n"
                "cell edge (%lf %lf %lf) (%lf %lf %lf) (edge type %d)\n"
                "circle (gc: norm_vec %lf %lf %lf\n"
                "        lon: norm_vec %lf %lf %lf\n"
                "        lat: z %lf north_is_out %d\n"
                "        point: vec %lf %lf %lf) (circle type %d)\n"
                "intersections points (%lf %lf %lf) (%lf %lf %lf)\n",
                cell_edge_coords[0][0],
                cell_edge_coords[0][1],
                cell_edge_coords[0][2],
                cell_edge_coords[1][0],
                cell_edge_coords[1][1],
                cell_edge_coords[1][2], (int)cell_edge_circle->type,
                  clipping_circle->data.gc.norm_vector[0],
                  clipping_circle->data.gc.norm_vector[1],
                  clipping_circle->data.gc.norm_vector[2],
                  clipping_circle->data.lon.norm_vector[0],
                  clipping_circle->data.lon.norm_vector[1],
                  clipping_circle->data.lon.norm_vector[2],
                  clipping_circle->data.lat.z,
                  clipping_circle->data.lat.north_is_out,
                  clipping_circle->data.p.vec[0],
                  clipping_circle->data.p.vec[1],
                  clipping_circle->data.p.vec[2], (int)(clipping_circle->type),
                  intersection[0][0], intersection[0][1], intersection[0][2],
                  intersection[1][0], intersection[1][1], intersection[1][2])

              // if both cell edge vertices are on the same side of the
              // clipping edge
              if (end_is_inside == start_is_inside) {

                // if the two intersection points are basically the same point
                if (points_are_identically(intersection[0], intersection[1])) {

                  num_edge_intersections = 1;

                } else {

                  // if the first intersection point is further away from the
                  // cell start edge vertex than the second one ->
                  // switch them
                  if (sq_len_diff_vec(cell_edge_coords[0], intersection[0]) >
                      sq_len_diff_vec(cell_edge_coords[0], intersection[1])) {

                    double temp_intersection[3];
                    temp_intersection[0] = intersection[1][0];
                    temp_intersection[1] = intersection[1][1];
                    temp_intersection[2] = intersection[1][2];
                    intersection[1][0] = intersection[0][0];
                    intersection[1][1] = intersection[0][1];
                    intersection[1][2] = intersection[0][2];
                    intersection[0][0] = temp_intersection[0];
                    intersection[0][1] = temp_intersection[1];
                    intersection[0][2] = temp_intersection[2];
                  }

                  num_edge_intersections = 2;
                }

              // one of the two cell edge vertices is on the clipping circle
              // => one of the two intersection points must be more or less
              //    identical to a vertex of the cell edge
              } else {

                double * cell_edge_coord = cell_edge_coords[end_is_inside == 2];
                double distances[2] = {
                  sq_len_diff_vec(cell_edge_coord, intersection[0]),
                  sq_len_diff_vec(cell_edge_coord, intersection[1])};

                // if both intersection points are nearly identical to the cell
                // edge vertex
                if ((distances[0] < yac_sq_angle_tol) &&
                    (distances[1] < yac_sq_angle_tol)) {

                  num_edge_intersections = 0;

                } else {

                  num_edge_intersections = 1;

                  // if the fist intersection points is closer than the second
                  if (distances[0] < distances[1]) {
                    intersection[0][0] = intersection[1][0];
                    intersection[0][1] = intersection[1][1];
                    intersection[0][2] = intersection[1][2];
                  }
                }
              }
            // if only one intersection point is on the cell edge
            } else if (is_on_edge[0] || is_on_edge[1]) {

              // if one of the two cell edge vertices is on the clipping edge
              if ((end_is_inside == 2) || (start_is_inside == 2)) {

                // we assume that the intersection point is the respective
                // cell edge vertex, which is on the clipping edge -> we do not
                // need this intersection point
                num_edge_intersections = 0;

              } else {

                // if the second intersection point is on the cell edge
                if (is_on_edge[1]) {
                  intersection[0][0] = intersection[1][0];
                  intersection[0][1] = intersection[1][1];
                  intersection[0][2] = intersection[1][2];
                }

                num_edge_intersections = 1;
              }

            // if none of the two intersection points is on the cell edge
            } else {
              num_edge_intersections = 0;
            }
            break;
          } // case one or two circle intersections
        } // switch (num_circle_intersections)

        // If an intersection was expected but none was found, the clipping
        // circle was most propably to close to an edge vertex. We can now
        // assume that the respective vertex is directly on the clipping circle.
        if ((one_intersect_expected)  && (num_edge_intersections == 0)) {

          if (yac_circle_compare_distances(
                cell_edge_coords[0], cell_edge_coords[1],
                clipping_circle) <= 0)
            start_is_inside = 2;
          else
            end_is_inside = 2;

          one_intersect_expected = 0;

          // if this is the first iteration
          if (cell_edge_idx == 0) first_start_is_inside = start_is_inside;
        }
      // else -> not edge intersection was excepted
      } else {
        num_edge_intersections = 0;
      }

      YAC_ASSERT(
        num_edge_intersections != -1,
        "ERROR(circle_clipping): internal error");

      // here we know the number of intersections and their location and we
      // know the relation of the cell edge vertices to the clipping edge
      // (start_is_inside and end_is_inside)

      // if the start cell edge vertex is outside -> dump it after clipping
      // is finished
      cell_edge_start->to_be_removed = start_is_inside == 0;

      // the easiest case is that we expected one intersection and got one
      if (one_intersect_expected) {

        // insert an intersection point in the cell point list in the
        // current edge
        struct point_list_element * intersect_point =
          get_free_point_list_element(cell);
        cell_edge_start->next = intersect_point;
        intersect_point->next = cell_edge_end;

        intersect_point->vec_coords[0] = intersection[0][0];
        intersect_point->vec_coords[1] = intersection[0][1];
        intersect_point->vec_coords[2] = intersection[0][2];

        intersect_point->edge_circle =
          (start_is_inside)?clipping_circle:cell_edge_circle;

      // if the clipping edge goes through both of the two cell edge vertices
      } else if ((start_is_inside == 2) && (end_is_inside == 2)) {

        // if one of the two edges is a circle of latitude while the other is
        // not, then we may have to replace the cell edge circle with the
        // clipping circle
        if ((cell_edge_circle->type == LAT_CIRCLE) ^
            (clipping_circle->type == LAT_CIRCLE)) {

          int clipping_circle_contains_north =
            yac_circle_contains_north_pole(clipping_circle);
          int same_inside_direction =
            clipping_circle_contains_north ==
            yac_circle_contains_north_pole(cell_edge_circle);
          int cell_edge_is_on_south_hemisphere =
            (cell_edge_end->vec_coords[2] < 0.0);
          int clipping_circle_is_lat = clipping_circle->type == LAT_CIRCLE;

          // after generating a truth table with these four logical values
          // I came up with the following formula to determine whether
          // the cell edge type should switch to the one of the clipping circle
          // or not
          if (same_inside_direction &&
              (cell_edge_is_on_south_hemisphere ^
               clipping_circle_contains_north ^
               clipping_circle_is_lat))
            cell_edge_start->edge_circle = clipping_circle;

        }

      // if we have no intersection, but the cell edge start vertex
      // is on the clipping edge
      } else if ((num_edge_intersections == 0) && (start_is_inside == 2)) {

        // if the end cell edge vertex is outside
        if (end_is_inside == 0)
          cell_edge_start->edge_circle = clipping_circle;

      // if we have two intersections (only happens if one of the two edges is
      // a circle of latitude while the other is not)
      } else if (num_edge_intersections == 2) {

        struct point_list_element * intersect_points[2] =
          {get_free_point_list_element(cell),
           get_free_point_list_element(cell)};

        // add two points between the current source edge vertices
        cell_edge_start->next = intersect_points[0];
        intersect_points[0]->next = intersect_points[1];
        intersect_points[1]->next = cell_edge_end;

        intersect_points[0]->vec_coords[0] = intersection[0][0];
        intersect_points[0]->vec_coords[1] = intersection[0][1];
        intersect_points[0]->vec_coords[2] = intersection[0][2];
        intersect_points[1]->vec_coords[0] = intersection[1][0];
        intersect_points[1]->vec_coords[1] = intersection[1][1];
        intersect_points[1]->vec_coords[2] = intersection[1][2];

        YAC_ASSERT(
          ((start_is_inside == 0) && (end_is_inside == 0)) ||
          ((start_is_inside == 1) && (end_is_inside == 1)),
          "ERROR: one cell edge vertex is on the clipping edge, therefore we "
          "should not have two intersections.")

        // if a and b are outside
        if ((start_is_inside == 0) && (end_is_inside == 0)) {
          intersect_points[0]->edge_circle = cell_edge_circle;
          intersect_points[1]->edge_circle = clipping_circle;

        // if a and b are inside
        } else /*if ((start_is_inside == 1) && (end_is_inside == 1))*/ {
          intersect_points[0]->edge_circle = clipping_circle;
          intersect_points[1]->edge_circle = cell_edge_circle;
        }

      // if we have one intersection even though the two cell edge vertices
      // are not on opposite sides of the clipping edge
      } else if (two_intersect_possible && (num_edge_intersections == 1)) {

        // ensure that both cell edge vertices are on the same side of the
        // clipping edge or that at least one cell vertex is directly on
        // the clipping edge
        YAC_ASSERT(
          (start_is_inside == end_is_inside) ||
          ((start_is_inside == 2) || (end_is_inside == 2)),
          "ERROR: unhandled intersection case")

        switch (MAX(start_is_inside, end_is_inside)) {

          // if both cell edge vertices are outside -> circle of latitude and
          // greate circle touch at a single vertex
          default:
          case (0): {
            // insert an intersection point in the cell point list in the
            // current edge
            struct point_list_element * intersect_point =
              get_free_point_list_element(cell);
            cell_edge_start->next = intersect_point;
            intersect_point->next = cell_edge_end;

            intersect_point->vec_coords[0] = intersection[0][0];
            intersect_point->vec_coords[1] = intersection[0][1];
            intersect_point->vec_coords[2] = intersection[0][2];

            intersect_point->edge_circle = clipping_circle;
            break;
          }

          // if both cell edge vertices are inside -> nothing to be done
          case (1): break;

          // if one of the two cell edge vertices is on the clipping edge
          // while the other is either inside or outside
          case (2): {
            // insert an intersection point in the cell point list in the
            // current edge
            struct point_list_element * intersect_point =
              get_free_point_list_element(cell);
            cell_edge_start->next = intersect_point;
            intersect_point->next = cell_edge_end;

            intersect_point->vec_coords[0] = intersection[0][0];
            intersect_point->vec_coords[1] = intersection[0][1];
            intersect_point->vec_coords[2] = intersection[0][2];

            // if the start cell edge vertex is on the clipping edge
            if (start_is_inside == 2) {
              // if the end cell edge vertex in on the outside
              if (end_is_inside == 0) {
                // cell_edge_start->edge_circle = cell_edge_circle;
                intersect_point->edge_circle = clipping_circle;
              } else {
                cell_edge_start->edge_circle = clipping_circle;
                intersect_point->edge_circle = cell_edge_circle;
              }
            // if the end cell edge vertex is on the clipping edge
            } else {
              // if the start cell edge vertex in on the outside
              if (start_is_inside == 0) {
                cell_edge_start->edge_circle = clipping_circle;
                intersect_point->edge_circle = cell_edge_circle;
              } else {
                // cell_edge_start->edge_circle = cell_edge_circle;
                intersect_point->edge_circle = clipping_circle;
              }
            }
          }
        }
      }

      cell_edge_start = cell_edge_end;
      cell_edge_end = cell_edge_end->next;
      start_is_inside = end_is_inside;

    } // for all cell edges

    // remove all points that are to be deleted
    num_cell_edges = remove_points(cell);
  }
}

static void point_list_clipping (
  struct point_list * source_list, struct point_list target_list, size_t nct) {

  struct yac_circle * clipping_circles[nct];

  struct point_list_element * curr_tgt_point = target_list.start;

  for (size_t i = 0; i < nct; ++i, curr_tgt_point = curr_tgt_point->next)
    clipping_circles[i] = curr_tgt_point->edge_circle;

  YAC_ASSERT(
    source_list->start != NULL,
    "ERROR(point_list_clipping): source cell without corners")

  // count the number of edges in the source cell
  size_t ncs = 1;
  for (struct point_list_element * curr = source_list->start->next,
                                 * end = source_list->start;
       curr != end; curr = curr->next, ++ncs);

  circle_clipping(source_list, ncs, clipping_circles, nct);
}

static void copy_point_list(struct point_list in, struct point_list * out) {

  reset_point_list(out);

  struct point_list_element * curr = in.start;

  if (curr == NULL) return;

  struct point_list_element * new_point_list = get_free_point_list_element(out);
  out->start = new_point_list;
  *new_point_list = *curr;
  curr = curr->next;

  do {

    new_point_list->next = get_free_point_list_element(out);
    new_point_list = new_point_list->next;
    *new_point_list = *curr;
    curr = curr->next;

  } while (curr != in.start);

  new_point_list->next = out->start;
}

void yac_cell_clipping (size_t N,
                        struct grid_cell * source_cell,
                        struct grid_cell target_cell,
                        struct grid_cell * overlap_buffer) {

  if (target_cell.num_corners < 2) {
    for (size_t n = 0; n < N; n++ ) overlap_buffer[n].num_corners = 0;
    return;
  }

  enum yac_cell_type tgt_cell_type = get_cell_type(target_cell);

  YAC_ASSERT(
   tgt_cell_type != MIXED_CELL,
   "invalid target cell type (cell contains edges consisting "
   "of great circles and circles of latitude)")

  // determine ordering of target cell corners
  int target_ordering = get_cell_points_ordering(target_cell);
  // if all corners of the target cell are on the same great circle
  if (!target_ordering) {
    for (size_t n = 0; n < N; n++ ) overlap_buffer[n].num_corners = 0;
    return;
  }

  size_t max_num_src_cell_corners = 0;
  for (size_t n = 0; n < N; ++n)
    if (source_cell[n].num_corners > max_num_src_cell_corners)
      max_num_src_cell_corners = source_cell[n].num_corners;
  struct yac_circle * circle_buffer =
    xmalloc(
      (target_cell.num_corners + max_num_src_cell_corners) *
       sizeof(*circle_buffer));
  struct yac_circle * src_circle_buffer =
    circle_buffer + target_cell.num_corners;

  // generate point list for target cell (clip cell)
  struct point_list target_list;
  init_point_list(&target_list);
  size_t nct =
    generate_point_list(
      &target_list, target_cell, target_ordering, circle_buffer);

  struct point_list source_list, temp_list;
  init_point_list(&temp_list);
  init_point_list(&source_list);

  // for all source cells
  for (size_t n = 0; n < N; n++ ) {

    overlap_buffer[n].num_corners = 0;

    enum yac_cell_type src_cell_type = get_cell_type(source_cell[n]);

    YAC_ASSERT(
      src_cell_type != MIXED_CELL,
      "invalid source cell type (cell contains edges consisting "
      "of great circles and circles of latitude)")

    if (source_cell[n].num_corners < 2) continue;

    // determine ordering of source cell corners
    int source_ordering = get_cell_points_ordering(source_cell[n]);

    // if all corners of the source cell are on the same great circle
    if (!source_ordering) continue;

    // generate point list for current source list
    size_t ncs =
      generate_point_list(
        &source_list, source_cell[n], source_ordering, src_circle_buffer);

    struct point_list * overlap;
    double fabs_tgt_coordinate_z = fabs(target_cell.coordinates_xyz[0][2]);
    double fabs_src_coordinate_z = fabs(source_cell[n].coordinates_xyz[0][2]);

    // in the case that source and target cell are both LAT_CELL's, than the
    // bigger one has to be the target cell
    // a similar problem occurs when the target cell is a LAT_CELL and the
    // source is a GREAT_CIRCLE_CELL which overlaps with the pole that is also
    // include in the target cell
    if (((tgt_cell_type == LAT_CELL) && (src_cell_type == GREAT_CIRCLE_CELL)) ||
        ((tgt_cell_type == LAT_CELL) && (src_cell_type == LAT_CELL) &&
         (fabs_tgt_coordinate_z > fabs_src_coordinate_z))) {

      copy_point_list(target_list, &temp_list);

      point_list_clipping(&temp_list, source_list, ncs);

      overlap = &temp_list;

    } else {

      point_list_clipping(&source_list, target_list, nct);

      overlap = &source_list;
    }

    generate_cell(overlap, overlap_buffer + n);

  }
  free(circle_buffer);
  free_point_list(&source_list);
  free_point_list(&target_list);
  free_point_list(&temp_list);
}

static void yac_lon_lat_cell_lat_clipping(
  struct grid_cell * cell, double z_upper_bound, double z_lower_bound,
  struct grid_cell * overlap_buffer) {

  double cell_upper_bound = cell->coordinates_xyz[0][2];
  double cell_lower_bound = cell->coordinates_xyz[2][2];

  int upper_idx[2], lower_idx[2];

  if (cell_upper_bound < cell_lower_bound) {
    double temp = cell_upper_bound;
    cell_upper_bound = cell_lower_bound;
    cell_lower_bound = temp;
    if (cell->edge_type[0] == LAT_CIRCLE_EDGE) {
      upper_idx[0] = 2;
      upper_idx[1] = 3;
      lower_idx[0] = 1;
      lower_idx[1] = 0;
    } else {
      upper_idx[0] = 2;
      upper_idx[1] = 1;
      lower_idx[0] = 3;
      lower_idx[1] = 0;
    }
  } else {
    if (cell->edge_type[0] == LAT_CIRCLE_EDGE) {
      upper_idx[0] = 0;
      upper_idx[1] = 1;
      lower_idx[0] = 3;
      lower_idx[1] = 2;
    } else {
      upper_idx[0] = 0;
      upper_idx[1] = 3;
      lower_idx[0] = 1;
      lower_idx[1] = 2;
    }
  }

  // if z_upper_bound and z_lower_bound are identical or
  // if cell does not overlap with latitude band
  if ((z_upper_bound == z_lower_bound) ||
      (cell_upper_bound <= z_lower_bound) ||
      (cell_lower_bound >= z_upper_bound)) {
    overlap_buffer->num_corners = 0;
    return;
  }

  if (overlap_buffer->array_size < 4) {
    overlap_buffer->coordinates_xyz =
      xmalloc(4 * sizeof(*(overlap_buffer->coordinates_xyz)));
    overlap_buffer->edge_type =
      xmalloc(4 * sizeof(*(overlap_buffer->edge_type)));
    overlap_buffer->array_size = 4;
  }

  memcpy(overlap_buffer->coordinates_xyz, cell->coordinates_xyz,
         4 * sizeof(*(cell->coordinates_xyz)));
  memcpy(overlap_buffer->edge_type, cell->edge_type,
         4 * sizeof(*(cell->edge_type)));
  overlap_buffer->num_corners = 4;


  double tmp_scale;
  double * p[2];
  if (fabs(cell_lower_bound) < fabs(cell_upper_bound)) {
    p[0] = cell->coordinates_xyz[lower_idx[0]];
    p[1] = cell->coordinates_xyz[lower_idx[1]];
    tmp_scale = cell->coordinates_xyz[lower_idx[0]][0] *
                cell->coordinates_xyz[lower_idx[0]][0] +
                cell->coordinates_xyz[lower_idx[0]][1] *
                cell->coordinates_xyz[lower_idx[0]][1];
  } else {
    p[0] = cell->coordinates_xyz[upper_idx[0]];
    p[1] = cell->coordinates_xyz[upper_idx[1]];
    tmp_scale = cell->coordinates_xyz[upper_idx[0]][0] *
                cell->coordinates_xyz[upper_idx[0]][0] +
                cell->coordinates_xyz[upper_idx[0]][1] *
                cell->coordinates_xyz[upper_idx[0]][1];
  }

  // if upper bound overlaps with cell
  if ((z_upper_bound < cell_upper_bound) &&
      (z_upper_bound > cell_lower_bound)) {

    double scale = sqrt((1.0 - z_upper_bound * z_upper_bound) / tmp_scale);

    overlap_buffer->coordinates_xyz[upper_idx[0]][0] = p[0][0] * scale;
    overlap_buffer->coordinates_xyz[upper_idx[0]][1] = p[0][1] * scale;
    overlap_buffer->coordinates_xyz[upper_idx[0]][2] = z_upper_bound;
    overlap_buffer->coordinates_xyz[upper_idx[1]][0] = p[1][0] * scale;
    overlap_buffer->coordinates_xyz[upper_idx[1]][1] = p[1][1] * scale;
    overlap_buffer->coordinates_xyz[upper_idx[1]][2] = z_upper_bound;
  }

  // if lower bound overlaps with cell
  if ((z_lower_bound < cell_upper_bound) &&
      (z_lower_bound > cell_lower_bound)) {

    double scale = sqrt((1.0 - z_lower_bound * z_lower_bound) / tmp_scale);

    overlap_buffer->coordinates_xyz[lower_idx[0]][0] = p[0][0] * scale;
    overlap_buffer->coordinates_xyz[lower_idx[0]][1] = p[0][1] * scale;
    overlap_buffer->coordinates_xyz[lower_idx[0]][2] = z_lower_bound;
    overlap_buffer->coordinates_xyz[lower_idx[1]][0] = p[1][0] * scale;
    overlap_buffer->coordinates_xyz[lower_idx[1]][1] = p[1][1] * scale;
    overlap_buffer->coordinates_xyz[lower_idx[1]][2] = z_lower_bound;
  }
}

static double get_closest_pole(struct point_list * cell_list) {

  struct point_list_element * curr = cell_list->start;
  struct point_list_element * start = cell_list->start;

  if (curr == NULL) return 1.0;

  double max_z = 0.0;

  do {

    double curr_z = curr->vec_coords[2];
    if (fabs(curr_z) > fabs(max_z)) max_z = curr_z;

    curr = curr->next;
  } while (curr != start);

  return (max_z > 0.0)?1.0:-1.0;
}

/*this routine is potentially being used in the CDO*/
void yac_cell_lat_clipping (size_t N,
                            struct grid_cell * cells,
                            double lat_bounds[2], // lat in rad
                            struct grid_cell * overlap_buffer) {

  double z_bounds[2] = {sin(lat_bounds[0]), sin(lat_bounds[1])};
  int upper_bound_idx = lat_bounds[0] < lat_bounds[1];
  int lower_bound_idx = upper_bound_idx ^ 1;
  int is_pole[2] = {fabs(fabs(lat_bounds[0]) - M_PI_2) < yac_angle_tol,
                    fabs(fabs(lat_bounds[1]) - M_PI_2) < yac_angle_tol};
  int upper_is_north_pole =
    is_pole[upper_bound_idx] && (lat_bounds[upper_bound_idx] > 0.0);
  int lower_is_south_pole =
    is_pole[lower_bound_idx] && (lat_bounds[lower_bound_idx] < 0.0);

  // if both bounds are nearly identical,
  // or if the upper bound is the south pole,
  // or if the lower bound is the north pole
  if ((fabs(lat_bounds[0] - lat_bounds[1]) < yac_angle_tol) ||
      (is_pole[upper_bound_idx] ^ upper_is_north_pole) ||
      (is_pole[lower_bound_idx] ^ lower_is_south_pole)) {

    for (size_t n = 0; n < N; ++n)
      overlap_buffer[n].num_corners = 0;
    return;
  }

  struct yac_circle lat_circle_buffer[2];
  struct yac_circle * lat_circles[2] =
    {&(lat_circle_buffer[0]), &(lat_circle_buffer[1])};
  size_t num_lat_circles = 0;
  if (!lower_is_south_pole)
    lat_circle_buffer[num_lat_circles++] =
      generate_lat_circle(z_bounds[lower_bound_idx], 0);
  if (!upper_is_north_pole)
    lat_circle_buffer[num_lat_circles++] =
      generate_lat_circle(z_bounds[upper_bound_idx], 1);

  size_t max_num_cell_corners = 0;
  for (size_t n = 0; n < N; ++n)
    if (cells[n].num_corners > max_num_cell_corners)
      max_num_cell_corners = cells[n].num_corners;
  struct yac_circle * circle_buffer =
    xmalloc(max_num_cell_corners * sizeof(*circle_buffer));

  struct point_list cell_list;
  init_point_list(&cell_list);

  // for all source cells
  for (size_t n = 0; n < N; n++) {

    if (cells[n].num_corners < 2) continue;

    overlap_buffer[n].num_corners = 0;

    enum yac_cell_type src_cell_type = get_cell_type(cells[n]);

    YAC_ASSERT(
      src_cell_type != MIXED_CELL,
      "invalid source cell type (cell contains edges consisting "
      "of great circles and circles of latitude)\n")

    if (src_cell_type == LON_LAT_CELL) {

      yac_lon_lat_cell_lat_clipping(
        cells + n, z_bounds[upper_bound_idx], z_bounds[lower_bound_idx],
        overlap_buffer + n);

    } else {

      int cell_ordering = get_cell_points_ordering(cells[n]);

      // generate point list for current source list
      size_t num_corners =
        generate_point_list(
          &cell_list, cells[n], cell_ordering, circle_buffer);

      // if the cell contains a pole, we need to add this pole
      double pole = get_closest_pole(&cell_list);
      if (yac_point_in_cell(
            (double[3]){0.0, 0.0, pole}, cells[n])) {

        int flag = 0;

        // if the cell contains the north pole
        if (pole > 0.0) {
          // use lower bound if upper bound is north pole
          double z_bound = z_bounds[upper_bound_idx ^ upper_is_north_pole];

          for (size_t i = 0; i < num_corners; ++i) {
            flag |= (z_bound < cells[n].coordinates_xyz[i][2]);
          }
        } else {
          // use upper bound if lower bound is south pole
          double z_bound = z_bounds[lower_bound_idx ^ lower_is_south_pole];

          for (size_t i = 0; i < num_corners; ++i) {
            flag |= (z_bound > cells[n].coordinates_xyz[i][2]);
          }
        }

        YAC_ASSERT(
          flag,
          "ERROR(yac_cell_lat_clipping): Latitude bounds are within a cell "
          "covering a pole, this is not supported. Increased grid resolution "
          "or widen lat bounds may help.")
      }

      circle_clipping(&cell_list, num_corners, lat_circles, num_lat_circles);

      generate_cell(&cell_list, overlap_buffer + n);
    }
  }

  free(circle_buffer);
  free_point_list(&cell_list);
}

/* ---------------------------------------------------- */

void yac_correct_weights ( size_t nSourceCells, double * weight ) {

  // number of iterations to get better accuracy of the weights
  enum {maxIter = 10};

  for ( size_t iter = 1; iter < maxIter; iter++ ) {

    double weight_sum = 0.0;

    for ( size_t n = 0; n < nSourceCells; n++ ) weight_sum += weight[n];

    if ( fabs(1.0 - weight_sum) < tol ) break;

    double scale = 1.0 / weight_sum;

    for ( size_t n = 0; n < nSourceCells; n++ ) weight[n] *= scale;
  }
}

/* ---------------------------------------------------- */

static int get_cell_points_ordering(struct grid_cell cell) {

  YAC_ASSERT(
    cell.num_corners > 2, "ERROR(get_cell_points_ordering): invalid cell")

  double edge_norm_vectors[cell.num_corners][3];
  double vertices[cell.num_corners][3];
  size_t num_corners = 0;

  for (size_t i = 0, j = cell.num_corners - 1; i < cell.num_corners; i++) {
    if (!points_are_identically(
           cell.coordinates_xyz[j], cell.coordinates_xyz[i])) {
      crossproduct_d(
        cell.coordinates_xyz[j], cell.coordinates_xyz[i],
        edge_norm_vectors[num_corners]);
      normalise_vector(edge_norm_vectors[num_corners]);
      vertices[num_corners][0] = cell.coordinates_xyz[i][0];
      vertices[num_corners][1] = cell.coordinates_xyz[i][1];
      vertices[num_corners][2] = cell.coordinates_xyz[i][2];
      ++num_corners;
      j = i;
    }
  }

  struct sin_cos_angle angle_sum = SIN_COS_ZERO;
  int cell_direction = 0;
  int valid_angles_count = 0;

  for (size_t i = 0, j = num_corners - 1; i < num_corners; j = i++) {

    struct sin_cos_angle edge_angle =
      get_vector_angle_2(edge_norm_vectors[j], edge_norm_vectors[i]);

    // if both edges are on the same plane
    if ((edge_angle.cos > 0.0) && (edge_angle.sin < yac_angle_tol))
      continue;

    valid_angles_count++;

    double edge_direction =
      edge_norm_vectors[j][0] * vertices[i][0] +
      edge_norm_vectors[j][1] * vertices[i][1] +
      edge_norm_vectors[j][2] * vertices[i][2];

    struct sin_cos_angle temp = angle_sum;
    if (edge_direction < 0.0)
      cell_direction -= sub_angles(temp, edge_angle, &angle_sum);
    else
      cell_direction += sum_angles(temp, edge_angle, &angle_sum);
  }

  if (valid_angles_count == 0) return 0;
  return (cell_direction >= 0)?1:-1;
}

static void init_point_list(struct point_list * list) {

  list->start = NULL;
  list->free_elements = NULL;
}

static void reset_point_list(struct point_list * list) {

  if (list->start != NULL) {

    struct point_list_element * temp = list->start->next;
    list->start->next = list->free_elements;
    list->free_elements = temp;

    list->start = NULL;
  }
}

static size_t remove_points(struct point_list * list) {

  struct point_list_element * curr = list->start;
  struct point_list_element * start = list->start;
  struct point_list_element * prev = NULL;

  if (curr == NULL) return 0;

  // find the first point that is not to be removed
  while(curr->to_be_removed) {
    prev = curr;
    curr = curr->next;
    if (curr == start) break;
  };

  // there is no point to remain
  if (curr->to_be_removed) {
    reset_point_list(list);
    return 0;
  }

  list->start = curr;
  start = curr;
  size_t num_remaining_points = 1;

  prev = curr;
  curr = curr->next;

  while (curr != start) {

    if (curr->to_be_removed) {
      prev->next = curr->next;
      curr->next = list->free_elements;
      list->free_elements = curr;
      curr = prev;
    } else {
      num_remaining_points++;
    }

    prev = curr;
    curr = curr->next;

  };

  if (list->start == list->start->next) {

    list->start->next = list->free_elements;
    list->free_elements = list->start;
    list->start = NULL;
    num_remaining_points = 0;
  }

  return num_remaining_points;
}

//! returns number of edges/corners
static size_t remove_zero_length_edges(struct point_list * list) {

  struct point_list_element * curr = list->start;
  struct point_list_element * start = list->start;

  do {

    if (!curr->to_be_removed) {

      // if both points are nearly identical (angle between them is very small)
      if (points_are_identically(curr->vec_coords, curr->next->vec_coords)) {
        curr->to_be_removed = 1;

      // check whether we can merge two lat circle edges
      } else if (curr->edge_circle->type == LAT_CIRCLE &&
                 curr->next->edge_circle->type == LAT_CIRCLE) {
        double temp_a = atan2(curr->vec_coords[1], curr->vec_coords[0]);
        double temp_b = atan2(curr->next->next->vec_coords[1],
                              curr->next->next->vec_coords[0]);
        if (fabs(get_angle(temp_a, temp_b)) < M_PI_2)
          curr->next->to_be_removed = 1;
      }
    }

    curr = curr->next;
  } while (curr != start);

  return remove_points(list);
}

static size_t generate_point_list(
  struct point_list * list, struct grid_cell cell, int cell_ordering,
  struct yac_circle * circle_buffer) {

  reset_point_list(list);

  if (cell.num_corners == 0) return 0;

  struct point_list_element * curr = get_free_point_list_element(list);

  list->start = curr;

  for (size_t i = 0; i < cell.num_corners; ++i, curr = curr->next) {

    curr->vec_coords[0] = cell.coordinates_xyz[i][0];
    curr->vec_coords[1] = cell.coordinates_xyz[i][1];
    curr->vec_coords[2] = cell.coordinates_xyz[i][2];
    yac_circle_generate(
      cell.coordinates_xyz[i], cell.coordinates_xyz[(i+1)%(cell.num_corners)],
      cell.edge_type[i], cell_ordering,
      ((curr->edge_circle = circle_buffer + i)));
    curr->next =
      ((i + 1) == cell.num_corners)?
        (list->start):get_free_point_list_element(list);
  }

  return remove_zero_length_edges(list);
}

static struct point_list_element *
get_free_point_list_element(struct point_list * list) {

  struct point_list_element * element;

  if (list->free_elements == NULL) {

    for (int i = 0; i < 7; ++i) {

      element = (struct point_list_element *)xmalloc(1 * sizeof(*element));

      element->next = list->free_elements;
      list->free_elements = element;
    }

    element = (struct point_list_element *)xmalloc(1 * sizeof(*element));

  } else {

    element = list->free_elements;
    list->free_elements = list->free_elements->next;
  }

  element->next = NULL;
  element->to_be_removed = 0;

  return element;
}

static void free_point_list(struct point_list * list) {

  struct point_list_element * element;

  if (list->start != NULL) {

    struct point_list_element * temp = list->free_elements;
    list->free_elements = list->start->next;
    list->start->next = temp;
  }

  while (list->free_elements != NULL) {

    element = list->free_elements;
    list->free_elements = element->next;
    free(element);
  }

  list->start = NULL;
  list->free_elements = NULL;
}

static void compute_norm_vector(double a[], double b[], double norm[]) {

  crossproduct_kahan(a, b, norm);

  YAC_ASSERT(
    (fabs(norm[0]) >= tol) || (fabs(norm[1]) >= tol) || (fabs(norm[2]) >= tol),
    "ERROR: a and b are identical -> no norm vector")

  normalise_vector(norm);
}

static int is_empty_gc_cell(struct point_list * list, size_t num_edges) {

  if (num_edges < 2) return 0;

  // if the polygon only has two vertices
  if (num_edges == 2)
    // the cell is empty if both edges have the same type
    // (lon and gc are counted as the same type)
    return !((list->start->edge_circle->type == LAT_CIRCLE) ^
             (list->start->next->edge_circle->type == LAT_CIRCLE));

  struct point_list_element * curr = list->start;

  // if there is at least one lat circle edge, the cell is not empty
  for (size_t i = 0; i < num_edges; ++i, curr = curr->next)
    if (curr->edge_circle->type == LAT_CIRCLE) return 0;

  // compute the norm vector for the first edge
  double norm_vec[3];
  compute_norm_vector(curr->vec_coords, curr->next->vec_coords, norm_vec);
  curr = curr->next->next;

  // check whether at least one vertex of the polygon is on a diffent plane
  // than the one defined by the first edge of the polygon
  for (size_t i = 0; i < num_edges - 2; ++i, curr = curr->next)
    if (fabs(dotproduct(norm_vec, curr->vec_coords)) > yac_angle_tol)
      return 0;

  return 1;
}

static inline enum yac_edge_type circle2edge_type(enum yac_circle_type type) {
  switch (type) {
    default:
    case (GREAT_CIRCLE): return GREAT_CIRCLE_EDGE;
    case (LON_CIRCLE): return LON_CIRCLE_EDGE;
    case (LAT_CIRCLE): return LAT_CIRCLE_EDGE;
  }
}

static void generate_cell(struct point_list * list,
                          struct grid_cell * cell) {

  if (list->start == NULL) {

    reset_point_list(list);
    cell->num_corners = 0;
    return;
  }

  size_t num_edges = remove_zero_length_edges(list);

  if (is_empty_gc_cell(list, num_edges)) {

    reset_point_list(list);
    cell->num_corners = 0;
    return;
  }

  if (num_edges > cell->array_size) {
    free(cell->coordinates_xyz);
    free(cell->edge_type);
    cell->coordinates_xyz = (double(*)[3])xmalloc(num_edges * sizeof(*cell->coordinates_xyz));
    cell->edge_type = (enum yac_edge_type *)xmalloc(num_edges * sizeof(*cell->edge_type));
    cell->array_size = num_edges;
  }
  cell->num_corners = num_edges;

  struct point_list_element * curr = list->start;

  for (size_t i = 0; i < num_edges; ++i) {

    cell->coordinates_xyz[i][0] = curr->vec_coords[0];
    cell->coordinates_xyz[i][1] = curr->vec_coords[1];
    cell->coordinates_xyz[i][2] = curr->vec_coords[2];
    cell->edge_type[i] = circle2edge_type(curr->edge_circle->type);

    curr = curr->next;
  }
}
