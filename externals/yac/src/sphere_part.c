/**
 * @file sphere_part.c
 *
 * @copyright Copyright  (C)  2014 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
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

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <float.h>

#include "sphere_part.h"
#include "geometry.h"
#include "grid.h"
#include "interval_tree.h"
#include "utils.h"
#include "ensure_array_size.h"

union I_list {
  struct {
    struct interval_node * head_node;
    size_t num_nodes;
  } ivt;
  size_t *list;
};

enum {
   I_list_tree_min_size = 2, //!< make I list into tree when list is
                             //!<  larger than this
};
enum {
  U_FLAG = 1,
  T_FLAG = 2,
};

enum yac_node_flags {
   U_IS_LEAF = 1,
   T_IS_LEAF = 2,
   I_IS_INTERVAL_TREE = 4,
};

struct sphere_part_node {

   int flags;

   union I_list I;
   void * U, * T;

   size_t I_size, U_size, T_size;

   struct sin_cos_angle I_angle;

   double gc_norm_vector[3];
};

struct bnd_sphere_part_search {
  struct sphere_part_node base_node;
  size_t * ids;
};

enum node_type {
  I_NODE_FULL = 0,
  I_NODE = 1,
  U_NODE = 2,
  T_NODE = 3,
};

struct temp_partition_data {
  struct bounding_circle bnd_circle; // has to be first entry
  size_t local_id;
  int node_type;
};

struct point_id_xyz {
  double coordinates_xyz[3]; // has to be the first entry
  size_t idx; // index of the point in the coordinates array passed to
              // the constructor
};

struct point_id_xyz_angle {
  struct point_id_xyz point;
  double cos_angle;
};

struct point_sphere_part_node {

   int flags;

   void * U, * T;

   size_t U_size, T_size;

   double gc_norm_vector[3];
};

struct point_sphere_part_search {

   struct point_sphere_part_node base_node;
   struct point_id_xyz * points;
   size_t max_tree_depth;
};

struct point_id_xyz_int32 {
  size_t idx;
  int32_t coordinate_xyz[3];
};

static void init_sphere_part_node(struct sphere_part_node * node) {

   node->flags = 0;
   node->I_size = 0;
   node->U_size = 0;
   node->T_size = 0;
   node->I_angle = SIN_COS_ZERO;
   node->gc_norm_vector[0] = 0;
   node->gc_norm_vector[1] = 0;
   node->gc_norm_vector[2] = 1;
}

static struct sphere_part_node * get_sphere_part_node() {

   struct sphere_part_node * node = xmalloc(1 * sizeof(*node));

   init_sphere_part_node(node);

   return node;
}

// computes the norm vector for the plane the partitions the data in half
// (more or less)
static inline void compute_gc_norm_vector(
  void * coords_data, size_t coords_size, size_t coords_count,
  double prev_gc_norm_vector[], double gc_norm_vector[]) {


   double balance_point[3] = {0.0,0.0,0.0};

   // compute balance point
   for (size_t i = 0; i < coords_count; ++i) {

     double * coords =
      (double*)(((unsigned char*)coords_data) + i * coords_size);
     balance_point[0] += coords[0];
     balance_point[1] += coords[1];
     balance_point[2] += coords[2];
   }

   if ((fabs(balance_point[0]) > 1e-9) ||
       (fabs(balance_point[1]) > 1e-9) ||
       (fabs(balance_point[2]) > 1e-9)) {
      normalise_vector(balance_point);
   } else {
     balance_point[0] = prev_gc_norm_vector[2];
     balance_point[1] = prev_gc_norm_vector[0];
     balance_point[2] = prev_gc_norm_vector[1];
   }

   crossproduct_kahan(
      balance_point, prev_gc_norm_vector, gc_norm_vector);

  if ((fabs(gc_norm_vector[0]) > 1e-9) ||
      (fabs(gc_norm_vector[1]) > 1e-9) ||
      (fabs(gc_norm_vector[2]) > 1e-9)) {
     normalise_vector(gc_norm_vector);
  } else {
     gc_norm_vector[0] = prev_gc_norm_vector[2];
     gc_norm_vector[1] = prev_gc_norm_vector[0];
     gc_norm_vector[2] = prev_gc_norm_vector[1];
  }
}

static void swap_partition_data(struct temp_partition_data *a, struct temp_partition_data *b) 
{ 
  struct temp_partition_data temp = *a; 
  *a = *b; 
  *b = temp; 
}

static size_t swap_node_type(struct temp_partition_data *part_data, size_t i, int node_type, size_t begin, size_t end)
{
  for (size_t j = begin; j < end; ++j)
    {
      if (part_data[j].node_type == node_type)
        {
          swap_partition_data(&part_data[i], &part_data[j]);
          return j + 1;
        }
    }

  return end + 1;
}

static void sort_partition_data(struct temp_partition_data *part_data, size_t I_FULL_size, size_t I_size, size_t U_size, size_t T_size)
{
   size_t I_begin = I_FULL_size;
   size_t U_begin = I_FULL_size + I_size;
   size_t T_begin = I_FULL_size + I_size + U_size;

   size_t I_end = I_begin + I_size;
   size_t U_end = U_begin + U_size;
   size_t T_end = T_begin + T_size;

   while (I_end > I_begin && part_data[I_end-1].node_type == I_NODE) I_end--;
   while (U_end > U_begin && part_data[U_end-1].node_type == U_NODE) U_end--;
   while (T_end > T_begin && part_data[T_end-1].node_type == T_NODE) T_end--;

   for (size_t i = 0; i < I_FULL_size; ++i)
     {
       if (part_data[i].node_type != I_NODE_FULL)
         {
           I_begin = swap_node_type(part_data, i, I_NODE_FULL, I_begin, T_end);
         }
     }

   I_begin = I_FULL_size;
   for (size_t i = I_begin; i < I_end; ++i)
     {
       if (part_data[i].node_type != I_NODE)
         {
           U_begin = swap_node_type(part_data, i, I_NODE, U_begin, U_end);
           if (U_begin > U_end)
             {
               T_begin = swap_node_type(part_data, i, I_NODE, T_begin, T_end);
             }
         }
     }

   U_begin = I_FULL_size + I_size;
   T_begin = I_FULL_size + I_size + U_size;
   for (size_t i = U_begin; i < U_end; ++i)
     {
       if (part_data[i].node_type != U_NODE)
         {
           T_begin = swap_node_type(part_data, i, U_NODE, T_begin, T_end);
         }
     }
}

static void partition_data (
   size_t * local_cell_ids, struct temp_partition_data * part_data,
   size_t num_cell_ids, size_t threshold, struct sphere_part_node * parent_node,
   double prev_gc_norm_vector[]) {

   if (num_cell_ids == 0) {
     parent_node->flags = U_IS_LEAF + T_IS_LEAF;
     return;
   }

   compute_gc_norm_vector(
     part_data, sizeof(*part_data), num_cell_ids,
     prev_gc_norm_vector, parent_node->gc_norm_vector);

   // partition data into cells that overlap with the great circle and cells
   // that are on side of the circle

   size_t I_FULL_size = 0;
   size_t I_size = 0;
   size_t U_size = 0;
   size_t T_size = 0;

   struct sin_cos_angle max_inc_angle = SIN_COS_ZERO;

   for (size_t i = 0; i < num_cell_ids; ++i) {

      struct bounding_circle curr_bnd_circle = part_data[i].bnd_circle;

      // get angle between the norm vector of the great circle and the base
      // point of the bounding circle
      struct sin_cos_angle angle =
        get_vector_angle_2(
          curr_bnd_circle.base_vector, parent_node->gc_norm_vector);

      // get the angle between between the plane of the great circle and base
      // point of the bounding circle
      struct sin_cos_angle diff_angle_gc =
        sin_cos_angle_new(fabs(angle.cos), angle.sin);

      // if the point intersects with the great circle
      if (compare_angles(diff_angle_gc, curr_bnd_circle.inc_angle) <= 0) {

         // if gc_norm_vector or -gc_norm_vector is in the bounding circle
         if (angle.sin < curr_bnd_circle.inc_angle.sin) {

            // set node type for current cell
            part_data[i].node_type = I_NODE_FULL;
            I_FULL_size++;

            max_inc_angle = SIN_COS_M_PI_2;

         } else {

            // set node type for current cell
            part_data[i].node_type = I_NODE;
            I_size++;

            struct sin_cos_angle inc_angle =
               sum_angles_no_check(diff_angle_gc, curr_bnd_circle.inc_angle);

            if (compare_angles(max_inc_angle, inc_angle) < 0)
               max_inc_angle = inc_angle;
         }

      // angle > M_PI_2
      } else if (angle.cos < 0.0) {

         // set node type for current cell
         part_data[i].node_type = U_NODE;
         U_size++;

      } else {

         // set node type for current cell
         part_data[i].node_type = T_NODE;
         T_size++;
      }
   }

   sort_partition_data(part_data, I_FULL_size, I_size, U_size, T_size);

   I_size += I_FULL_size;
   parent_node->I_size = I_size;
   parent_node->U_size = U_size;
   parent_node->T_size = T_size;

   // if max_inc_angle > PI/2
   if (compare_angles(max_inc_angle, SIN_COS_M_PI_2) >= 0) {
      parent_node->I_angle = SIN_COS_M_PI_2;
   } else {
      parent_node->I_angle = max_inc_angle;
   }

   if (I_size > 0) {

      if (I_size > I_list_tree_min_size) {

         assert(sizeof(struct interval_node) > sizeof(size_t));
         struct interval_node * head_node = xmalloc(I_size * sizeof(*head_node));
         parent_node->I.ivt.head_node = head_node;
         parent_node->I.ivt.num_nodes = I_size;

         // for all bounding circles that include gc_norm_vector or
         // -gc_norm_vector
         for (size_t i = 0; i < I_FULL_size; ++i) {

            head_node[i].range.left = -M_PI;
            head_node[i].range.right = M_PI;
            head_node[i].value = part_data[i].local_id;
         }

         for (size_t i = I_FULL_size; i < I_size; ++i) {

            double GCp[3], bVp[3];
            struct bounding_circle curr_bnd_circle = part_data[i].bnd_circle;

            // project the base vector of the current bounding circle onto the
            // current partitioning great circle
            crossproduct_kahan(parent_node->gc_norm_vector,
                            curr_bnd_circle.base_vector, GCp);
            crossproduct_kahan(GCp, parent_node->gc_norm_vector, bVp);
            YAC_ASSERT(
              (fabs(bVp[0]) > 1e-9) || (fabs(bVp[1]) > 1e-9) ||
              (fabs(bVp[2]) > 1e-9),
              "ERROR(partition_data): "
              "projected vector is nearly identical to gc_norm_vector")
            normalise_vector(bVp);

            // the inc_angle of the bounding circle also needs to be projected
            // onto the great circle
            // (the formula for this is:
            //    new_inc_angle=acos(sin(inc_angle)/cos(base_angle))
            // To find this formula, you have to determine the longitude of
            // of longitude circle that has exactly one intersection with the
            // bounding circle.
            // (see https://www.dropbox.com/s/82bxzno1mjogbtf/function_2nd_try.pdf)
            double bnd_circle_lat_cos =
              get_vector_angle_2(bVp, curr_bnd_circle.base_vector).cos;
            double inc_angle =
              M_PI_2 - acos(curr_bnd_circle.inc_angle.sin/bnd_circle_lat_cos);

            // By definition the previous norm vector is on the current great
            // great circle. We use this as a reference point for I

            // compute "distance" of the projected base vector to the reference
            // point
            double base_angle = get_vector_angle(bVp, prev_gc_norm_vector);

            head_node[i].range.left = base_angle - inc_angle;
            head_node[i].range.right = base_angle + inc_angle;
            head_node[i].value = part_data[i].local_id;
         }

         yac_generate_interval_tree(head_node, I_size);
         parent_node->flags |= I_IS_INTERVAL_TREE;
      } else {
         for (size_t i = 0; i < I_size; ++i)
            local_cell_ids[i] = part_data[i].local_id;
         parent_node->I.list = (void*)local_cell_ids;
      }
   } else
      parent_node->I.list = NULL;

   part_data += I_size;
   local_cell_ids += I_size;

   // check whether the lists are small enough (if not -> partition again)
   if (U_size <= threshold) {

      for (size_t i = 0; i < U_size; ++i)
         local_cell_ids[i] = part_data[i].local_id;
      parent_node->U = (void*)local_cell_ids;
      parent_node->flags |= U_IS_LEAF;

   } else {

      parent_node->U = get_sphere_part_node();
      partition_data(local_cell_ids, part_data, U_size, threshold,
                     parent_node->U, parent_node->gc_norm_vector);
   }
   local_cell_ids += U_size;
   part_data += U_size;

   if (T_size <= threshold) {

      for (size_t i = 0; i < T_size; ++i)
         local_cell_ids[i] = part_data[i].local_id;
      parent_node->T = (void*)local_cell_ids;
      parent_node->flags |= T_IS_LEAF;
      local_cell_ids += T_size;

   } else {

      parent_node->T = get_sphere_part_node();
      partition_data(local_cell_ids, part_data, T_size, threshold,
                     parent_node->T, parent_node->gc_norm_vector);
   }
}

static int compare_point_idx_xyz(void const * a, void const * b) {
  return (((struct point_id_xyz *)a)->idx > ((struct point_id_xyz *)b)->idx) -
         (((struct point_id_xyz *)a)->idx < ((struct point_id_xyz *)b)->idx);
}

static struct point_sphere_part_node * partition_point_data (
  struct point_id_xyz * points, size_t num_points, size_t threshold,
  double prev_gc_norm_vector[], size_t curr_tree_depth,
  size_t * max_tree_depth) {

  if (curr_tree_depth > *max_tree_depth) *max_tree_depth = curr_tree_depth;

  struct point_sphere_part_node * node = xmalloc(1 * sizeof(*node));
  double * gc_norm_vector = &(node->gc_norm_vector[0]);

  compute_gc_norm_vector(
    points, sizeof(*points), num_points, prev_gc_norm_vector, gc_norm_vector);

  // angle between a point and the great circle plane
  // acos(dot(gc_norm_vector, point_xyz)) = angle(gc_norm_vector, point_xyz)
  // acos(dot(gc_norm_vector, point_xyz)) - PI/2 = angle(gc_plane, point_xyz)
  // dot <= 0.0    -> U list
  // dot >  0.0    -> T list

  struct point_id_xyz * left = points, * right = points + num_points - 1;

  // sort such that all points for the U list come first, followed by the
  // elements of the T list
  while (1) {
    // find element that does not belong into U-list
    while (left <= right) {
      double * curr_coordinates_xyz = &(left->coordinates_xyz[0]);
      double dot = curr_coordinates_xyz[0] * gc_norm_vector[0] +
                   curr_coordinates_xyz[1] * gc_norm_vector[1] +
                   curr_coordinates_xyz[2] * gc_norm_vector[2];

      // if (angle < M_PI_2)
      if (dot > 0.0) break;
      ++left;
    };

    // find element that does not belong into T-list
    while (left < right) {
      double * curr_coordinates_xyz = &(right->coordinates_xyz[0]);
      double dot = curr_coordinates_xyz[0] * gc_norm_vector[0] +
                   curr_coordinates_xyz[1] * gc_norm_vector[1] +
                   curr_coordinates_xyz[2] * gc_norm_vector[2];

      // if (angle >= M_PI_2)
      if (dot <= 0.0) break;
      --right;
    }

    if (left < right) {
      struct point_id_xyz tmp_point = *left;
      *left = *right;
      *right = tmp_point;
      ++left;
      --right;
    } else {
      break;
    }
  }

  size_t U_size = left - points;
  size_t T_size = num_points - U_size;

  node->U_size = U_size;
  node->T_size = T_size;
  node->flags = 0;

  // check whether the lists are small enough (if not -> partition again)
  if ((U_size <= threshold) || (U_size == num_points)) {

    node->U = points;
    node->flags |= U_IS_LEAF;
    qsort(points, U_size, sizeof(*points), compare_point_idx_xyz);

  } else {

    node->U = partition_point_data(points, U_size, threshold, gc_norm_vector,
                                   curr_tree_depth + 1, max_tree_depth);
  }

  if ((T_size <= threshold) || (T_size == num_points)) {

    node->T = points + U_size;
    node->flags |= T_IS_LEAF;
    qsort(points + U_size, T_size, sizeof(*points), compare_point_idx_xyz);

  } else {

    node->T =
      partition_point_data(points + U_size, T_size, threshold, gc_norm_vector,
                           curr_tree_depth + 1, max_tree_depth);
  }

  return node;
}

struct bnd_sphere_part_search * yac_bnd_sphere_part_search_new(
  struct bounding_circle * circles, size_t num_circles) {

  struct bnd_sphere_part_search * search = xmalloc(1 * sizeof(*search));

  double gc_norm_vector[3] = {0.0,0.0,1.0};

  init_sphere_part_node(&(search->base_node));

  size_t * ids = xmalloc(num_circles * sizeof(*ids));
  search->ids = ids;

  struct temp_partition_data * part_data =
    xmalloc(num_circles * sizeof(*part_data));
   for (size_t i = 0; i < num_circles; ++i) {
      part_data[i].bnd_circle = circles[i];
      part_data[i].local_id = i;
   }

   partition_data(ids, part_data, num_circles, I_list_tree_min_size,
                  &(search->base_node), gc_norm_vector);

   free(part_data);

   return search;
}

static inline int compare_points_int32_coord(
  const void * a,const void * b) {

  struct point_id_xyz_int32 * point_a = (struct point_id_xyz_int32 *)a;
  struct point_id_xyz_int32 * point_b = (struct point_id_xyz_int32 *)b;

  int ret;

  for (int i = 0; i < 3; ++i)
    if ((ret = (point_a->coordinate_xyz[i] > point_b->coordinate_xyz[i]) -
               (point_a->coordinate_xyz[i] < point_b->coordinate_xyz[i])))
      return ret;
  return 0;
}

static struct point_id_xyz *
  get_unique_points(
    size_t * num_points, coordinate_pointer coordinates_xyz,
    yac_int const * ids, int const * mask) {

  struct point_id_xyz_int32 * points_int32 =
    xmalloc(*num_points * sizeof(*points_int32));

  double const scale = (double)(2 << 21);

  size_t num_unmasked_points;

  if (mask == NULL) {
    num_unmasked_points = *num_points;
    for (size_t i = 0; i < num_unmasked_points; ++i) {

      points_int32[i].idx = i;
      for (size_t j = 0; j < 3; ++j)
        points_int32[i].coordinate_xyz[j] =
          (int32_t)round(coordinates_xyz[i][j] * scale);
    }
  } else {
    num_unmasked_points = 0;
    for (size_t i = 0; i < *num_points; ++i) {

      if (!mask[i]) continue;
      points_int32[num_unmasked_points].idx = i;
      for (size_t j = 0; j < 3; ++j)
        points_int32[num_unmasked_points].coordinate_xyz[j] =
          (int32_t)round(coordinates_xyz[i][j] * scale);
      num_unmasked_points++;
    }
  }

  // sort points
  qsort(points_int32, num_unmasked_points,
        sizeof(*points_int32), compare_points_int32_coord);

  struct point_id_xyz_int32 dummy;
  dummy.idx = SIZE_MAX;
  dummy.coordinate_xyz[0] = INT32_MAX;
  dummy.coordinate_xyz[1] = INT32_MAX;
  dummy.coordinate_xyz[2] = INT32_MAX;
  struct point_id_xyz_int32 * prev = &dummy, * curr = points_int32;
  yac_int prev_id = XT_INT_MAX;
  size_t new_num_points = 0;
  for (size_t i = 0; i < num_unmasked_points; ++i, ++curr) {

    size_t curr_idx = curr->idx;
    if (compare_points_int32_coord(prev, curr)) {
      prev = points_int32 + new_num_points++;
      if (prev != curr) *prev = *curr;
      prev_id = ids[curr_idx];
    } else {
      yac_int curr_id = ids[curr_idx];
      if (curr_id > prev_id) {
        prev_id = curr_id;
        prev->idx = curr_idx;
      }
    }
  }

  struct point_id_xyz * points = xmalloc(new_num_points * sizeof(*points));
  for (size_t i = 0; i < new_num_points; ++i) {
    size_t curr_idx = points_int32[i].idx;
    points[i].idx = curr_idx;
    memcpy(
      points[i].coordinates_xyz, coordinates_xyz[curr_idx], 3 * sizeof(double));
  }
  *num_points = new_num_points;
  free(points_int32);
  return points;
}

struct point_sphere_part_search * yac_point_sphere_part_search_new (
  size_t num_points, coordinate_pointer coordinates_xyz, yac_int const * ids) {

  if (num_points == 0) return NULL;

  struct point_sphere_part_search * search = xmalloc(1 * sizeof(*search));
  struct point_id_xyz * points =
    get_unique_points(&num_points, coordinates_xyz, ids, NULL);
  search->points = points;

  size_t max_tree_depth = 0;

  // emperical measurements have given a threshold for the leaf size of 2
  struct point_sphere_part_node * tmp_node =
    partition_point_data(
      points, num_points, I_list_tree_min_size, (double[3]){0.0,0.0,1.0},
      1, &max_tree_depth);

  search->base_node = *tmp_node;
  search->max_tree_depth = max_tree_depth;
  free(tmp_node);

  return search;
}

struct point_sphere_part_search * yac_point_sphere_part_search_mask_new (
  size_t num_points, coordinate_pointer coordinates_xyz,
  yac_int const * ids, int const * mask) {

  if (num_points == 0) return NULL;

  struct point_id_xyz * points =
    get_unique_points(&num_points, coordinates_xyz, ids, mask);

  struct point_sphere_part_search * search = xmalloc(1 * sizeof(*search));
  search->points = xrealloc(points, num_points * sizeof(*points));

  size_t max_tree_depth = 0;

  // emperical measurements have given a threshold for the leaf size of 2
  struct point_sphere_part_node * tmp_node =
    partition_point_data(
      search->points, num_points, I_list_tree_min_size,
      (double[3]){0.0,0.0,1.0}, 1, &max_tree_depth);

  search->base_node = *tmp_node;
  search->max_tree_depth = max_tree_depth;
  free(tmp_node);

  return search;
}

static void search_bnd_circle_I_node(
  struct sphere_part_node * node, struct bounding_circle bnd_circle,
  size_t ** restrict overlap_cells, size_t * overlap_cells_array_size,
  size_t * restrict num_overlap_cells,
  struct overlaps * search_interval_tree_buffer, double prev_gc_norm_vector[]) {

  if (node->flags & I_IS_INTERVAL_TREE) {

    struct sin_cos_angle angle =
      get_vector_angle_2(
        bnd_circle.base_vector, node->gc_norm_vector);
    angle.cos = fabs(angle.cos);

    // if gc_norm_vector or -gc_norm_vector is in the bounding circle
    if (compare_angles(angle, bnd_circle.inc_angle) <= 0) {
      ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                        *num_overlap_cells + node->I.ivt.num_nodes);
      for (size_t i = 0; i < node->I.ivt.num_nodes; ++i)
        (*overlap_cells)[(*num_overlap_cells)+i] =
          node->I.ivt.head_node[i].value;
      *num_overlap_cells += node->I.ivt.num_nodes;
      return;
    }

    double GCp[3], bVp[3];

    // project the base vector of the current bounding circle onto the
    // current great circle
    crossproduct_kahan(node->gc_norm_vector, bnd_circle.base_vector, GCp);
    crossproduct_kahan(GCp, node->gc_norm_vector, bVp);
    YAC_ASSERT(
      (fabs(bVp[0]) > 1e-9) || (fabs(bVp[1]) > 1e-9) || (fabs(bVp[2]) > 1e-9),
      "ERROR(search_bnd_circle_I_node): "
      "projected vector is nearly identical to gc_norm_vector")
    normalise_vector(bVp);

    // the inc_angle of the bounding circle also needs to be projected
    // onto the great circle (see routine partition data for a detailed
    // explanation)
    double bnd_circle_lat_cos =
      get_vector_angle_2(bVp, bnd_circle.base_vector).cos;
    double inc_angle =
      M_PI_2 - acos(bnd_circle.inc_angle.sin/bnd_circle_lat_cos);

    // compute "distance" of the projected base vector to the reference point
    double base_angle = get_vector_angle(bVp, prev_gc_norm_vector);

    search_interval_tree_buffer->num_overlaps = 0;

    yac_search_interval_tree(
      node->I.ivt.head_node, node->I.ivt.num_nodes,
      (struct interval){
        .left = base_angle - inc_angle,
        .right = base_angle + inc_angle},
      search_interval_tree_buffer);

    ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                      *num_overlap_cells +
                      search_interval_tree_buffer->num_overlaps);

    for (size_t i = 0; i < search_interval_tree_buffer->num_overlaps; ++i)
      (*overlap_cells)[(*num_overlap_cells)+i] =
        node->I.ivt.head_node[
          search_interval_tree_buffer->overlap_iv[i]].value;

    *num_overlap_cells += search_interval_tree_buffer->num_overlaps;

  } else {

    ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                      *num_overlap_cells + node->I_size);
    memcpy(*overlap_cells + *num_overlap_cells, node->I.list,
    node->I_size * sizeof(**overlap_cells));
    *num_overlap_cells += node->I_size;
  }
}

//! TODO change to iterative implementation and allocate overlap_cells first
static void search_big_bnd_circle(
  struct sphere_part_node * node, struct bounding_circle bnd_circle,
  size_t ** restrict overlap_cells, size_t * overlap_cells_array_size,
  size_t * restrict num_overlap_cells,
  struct overlaps * search_interval_tree_buffer, double prev_gc_norm_vector[]) {

  if (node->flags & T_IS_LEAF) {

    ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                      *num_overlap_cells + node->T_size);
    memcpy(*overlap_cells + *num_overlap_cells, node->T,
           node->T_size * sizeof(**overlap_cells));
    *num_overlap_cells += node->T_size;

  } else {
    search_big_bnd_circle(
       node->T, bnd_circle, overlap_cells, overlap_cells_array_size,
       num_overlap_cells, search_interval_tree_buffer, node->gc_norm_vector);
  }

  if (node->flags & U_IS_LEAF) {

    ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                      *num_overlap_cells + node->U_size);
    memcpy(*overlap_cells + *num_overlap_cells, node->U,
           node->U_size * sizeof(**overlap_cells));
    *num_overlap_cells += node->U_size;

  } else {
    search_big_bnd_circle(
       node->U, bnd_circle, overlap_cells, overlap_cells_array_size,
       num_overlap_cells, search_interval_tree_buffer, node->gc_norm_vector);
  }

  search_bnd_circle_I_node(
     node, bnd_circle, overlap_cells, overlap_cells_array_size,
     num_overlap_cells, search_interval_tree_buffer, prev_gc_norm_vector);
}

static void search_small_bnd_circle(
  struct sphere_part_node * node, struct bounding_circle bnd_circle,
  size_t ** restrict overlap_cells, size_t * overlap_cells_array_size,
  size_t * restrict num_overlap_cells,
  struct overlaps * search_interval_tree_buffer, double prev_gc_norm_vector[]) {

   double dot = bnd_circle.base_vector[0] * node->gc_norm_vector[0] +
                bnd_circle.base_vector[1] * node->gc_norm_vector[1] +
                bnd_circle.base_vector[2] * node->gc_norm_vector[2];

   // angle < M_PI_2 + bnd_circle.inc_angle
   if (dot > - bnd_circle.inc_angle.sin) {

      if (node->flags & T_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->T_size);
         memcpy(*overlap_cells + *num_overlap_cells, node->T,
                node->T_size * sizeof(**overlap_cells));
         *num_overlap_cells += node->T_size;

      } else {
         search_small_bnd_circle(
            node->T, bnd_circle, overlap_cells, overlap_cells_array_size,
            num_overlap_cells, search_interval_tree_buffer,
            node->gc_norm_vector);
      }
   }

   // angle > M_PI_2 - bnd_circle.inc_angle
   if (dot < bnd_circle.inc_angle.sin) {

      if (node->flags & U_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->U_size);
         memcpy(*overlap_cells + *num_overlap_cells, node->U,
                node->U_size * sizeof(**overlap_cells));
         *num_overlap_cells += node->U_size;

      } else {
         search_small_bnd_circle(
            node->U, bnd_circle, overlap_cells, overlap_cells_array_size,
            num_overlap_cells, search_interval_tree_buffer,
            node->gc_norm_vector);
      }
   }

   struct sin_cos_angle angle_sum =
      sum_angles_no_check(node->I_angle, bnd_circle.inc_angle);

   // if (I_angle + inc_angle > PI/2) ||
   //    (fabs(angle - M_PI_2) <= (I_angle + inc_angle))
   //
   // assumtion:
   //   I_angle >= 0 && I_angle <= PI/2
   //   inc_angle >= 0 && inc_angle <= PI/2
   //   angle >= 0 && angle <= PI
   //
   //   => I_angle + inc_angle >= 0 && I_angle + inc_angle <= PI
   //
   // I_angle + inc_angle >= PI/2
   //
   // fabs(angle - M_PI_2) <= (I_angle + inc_angle)
   // => sin(fabs(angle - M_PI_2)) <= sin(I_angle + inc_angle)
   //    this is wrong for (I_angle + inc_angle) > PI/2, however that case is
   //    already covered by the first condition
   // => fabs(cos(angle)) <= sin(I_angle + inc_angle)
   if (((angle_sum.sin < 0.0) || (angle_sum.cos <= 0.0)) ||
       (fabs(dot) <= angle_sum.sin)) {
      search_bnd_circle_I_node(
         node, bnd_circle, overlap_cells, overlap_cells_array_size,
         num_overlap_cells, search_interval_tree_buffer, prev_gc_norm_vector);
   }
}

static void search_bnd_circle(struct sphere_part_node * node,
                              struct bounding_circle bnd_circle,
                              size_t ** restrict overlap_cells,
                              size_t * overlap_cells_array_size,
                              size_t * restrict num_overlap_cells,
                              struct overlaps * search_interval_tree_buffer,
                              double prev_gc_norm_vector[]) {

  // if the bounding circle has an angle in the range of [0;PI/2[
  if (bnd_circle.inc_angle.cos > 0.0)
    search_small_bnd_circle(
      node, bnd_circle, overlap_cells, overlap_cells_array_size,
      num_overlap_cells, search_interval_tree_buffer, prev_gc_norm_vector);
  else
    search_big_bnd_circle(
      node, bnd_circle, overlap_cells, overlap_cells_array_size,
      num_overlap_cells, search_interval_tree_buffer, prev_gc_norm_vector);
}

static inline void check_leaf_NN(
  struct point_id_xyz * points, size_t num_points,
  double * point_coordinates_xyz, struct sin_cos_angle * best_angle,
  double (**result_coordinates_xyz)[3],
  size_t * result_coordinates_xyz_array_size, size_t ** local_point_ids,
  size_t * local_point_ids_array_size, size_t total_num_local_point_ids,
  size_t * num_local_point_ids) {

  size_t * local_point_ids_ = *local_point_ids;
  size_t local_point_ids_array_size_ = *local_point_ids_array_size;
  double (*result_coordinates_xyz_)[3];
  size_t result_coordinates_xyz_array_size_;
  size_t num_local_point_ids_ = *num_local_point_ids;

  if (result_coordinates_xyz != NULL) {
    result_coordinates_xyz_ = *result_coordinates_xyz;
    result_coordinates_xyz_array_size_ = *result_coordinates_xyz_array_size;
    ENSURE_ARRAY_SIZE(
      result_coordinates_xyz_, result_coordinates_xyz_array_size_,
      total_num_local_point_ids + num_local_point_ids_ + num_points);
    *result_coordinates_xyz = result_coordinates_xyz_;
    *result_coordinates_xyz_array_size = result_coordinates_xyz_array_size_;
    result_coordinates_xyz_ += total_num_local_point_ids;
  }
  ENSURE_ARRAY_SIZE(
    local_point_ids_, local_point_ids_array_size_,
    total_num_local_point_ids + num_local_point_ids_ + num_points);
  *local_point_ids = local_point_ids_;
  *local_point_ids_array_size = local_point_ids_array_size_;
  local_point_ids_ += total_num_local_point_ids;


  // check leaf for results
  for (size_t i = 0; i < num_points; ++i) {

    struct sin_cos_angle curr_angle =
      get_vector_angle_2(
        points[i].coordinates_xyz, point_coordinates_xyz);
    int compare = compare_angles(curr_angle, *best_angle);

    // if the point is worse than the currently best point
    if (compare > 0) continue;

    // if we found a better point
    if (compare < 0) {

      *best_angle = curr_angle;
      num_local_point_ids_ = 1;
      if (result_coordinates_xyz != NULL) {
        result_coordinates_xyz_[0][0] = points[i].coordinates_xyz[0];
        result_coordinates_xyz_[0][1] = points[i].coordinates_xyz[1];
        result_coordinates_xyz_[0][2] = points[i].coordinates_xyz[2];
      }
      local_point_ids_[0] = points[i].idx;

    } else {

      if (result_coordinates_xyz != NULL) {
        result_coordinates_xyz_[num_local_point_ids_][0] =
          points[i].coordinates_xyz[0];
        result_coordinates_xyz_[num_local_point_ids_][1] =
          points[i].coordinates_xyz[1];
        result_coordinates_xyz_[num_local_point_ids_][2] =
          points[i].coordinates_xyz[2];
      }
      local_point_ids_[num_local_point_ids_] = points[i].idx;
      num_local_point_ids_++;
    }
  }

  *num_local_point_ids = num_local_point_ids_;
}

static void point_search_NN(
  struct bounding_circle * bnd_circle, double (**result_coordinates_xyz)[3],
  size_t * result_coordinates_xyz_array_size, size_t ** local_point_ids,
  size_t * local_point_ids_array_size, size_t total_num_local_point_ids,
  size_t * num_local_point_ids, double * dot_stack,
  struct point_sphere_part_node ** node_stack,
  int * flags, size_t curr_tree_depth) {

  double * point_coordinates_xyz = bnd_circle->base_vector;
  struct sin_cos_angle best_angle = bnd_circle->inc_angle;

  double dot = dot_stack[curr_tree_depth];
  struct point_sphere_part_node * node = node_stack[curr_tree_depth];
  int skip_U = flags[curr_tree_depth] & U_FLAG;
  int skip_T = flags[curr_tree_depth] & T_FLAG;

  do {

    if (!skip_U) {

      flags[curr_tree_depth] |= U_FLAG;

      // angle + inc_angle >= M_PI_2
      if ((dot < best_angle.sin) | (best_angle.cos <= 0.0)) {

        if (node->flags & U_IS_LEAF) {

          check_leaf_NN(
            (struct point_id_xyz *)(node->U), node->U_size,
            point_coordinates_xyz, &best_angle, result_coordinates_xyz,
            result_coordinates_xyz_array_size, local_point_ids,
            local_point_ids_array_size, total_num_local_point_ids,
            num_local_point_ids);

        } else {

          // traverse down one level
          ++curr_tree_depth;
          node = (struct point_sphere_part_node *)(node->U);
          dot = node->gc_norm_vector[0] * point_coordinates_xyz[0] +
                node->gc_norm_vector[1] * point_coordinates_xyz[1] +
                node->gc_norm_vector[2] * point_coordinates_xyz[2];
          dot_stack[curr_tree_depth] = dot;
          node_stack[curr_tree_depth] = node;
          flags[curr_tree_depth] = 0;
          // skip_U = 0;
          skip_T = 0;
          continue;
        }
      }
    }

    if (!skip_T) {

      flags[curr_tree_depth] = U_FLAG + T_FLAG;

      // angle - inc_angle < M_PI_2
      if ((dot > - best_angle.sin) || (best_angle.cos <= 0.0)) {

        if (node->flags & T_IS_LEAF) {

           check_leaf_NN(
             (struct point_id_xyz *)(node->T), node->T_size,
             point_coordinates_xyz, &best_angle, result_coordinates_xyz,
             result_coordinates_xyz_array_size, local_point_ids,
             local_point_ids_array_size, total_num_local_point_ids,
             num_local_point_ids);

        } else {

           // traverse down one level
           ++curr_tree_depth;
           node = (struct point_sphere_part_node *)(node->T);
           dot = node->gc_norm_vector[0] * point_coordinates_xyz[0] +
                 node->gc_norm_vector[1] * point_coordinates_xyz[1] +
                 node->gc_norm_vector[2] * point_coordinates_xyz[2];
           dot_stack[curr_tree_depth] = dot;
           node_stack[curr_tree_depth] = node;
           flags[curr_tree_depth] = 0;
           skip_U = 0;
           // skip_T = 0;
           continue;
        }
      }
    }

    if (curr_tree_depth == 0) break;

    // go up one level in the tree

    curr_tree_depth--;
    dot = dot_stack[curr_tree_depth];
    node = node_stack[curr_tree_depth];
    skip_U = flags[curr_tree_depth] & U_FLAG;
    skip_T = flags[curr_tree_depth] & T_FLAG;

  } while (1);

  bnd_circle->inc_angle = best_angle;
}

static int point_check_bnd_circle(
   struct point_sphere_part_node * node, struct bounding_circle bnd_circle) {

   double dot = node->gc_norm_vector[0]*bnd_circle.base_vector[0] +
                node->gc_norm_vector[1]*bnd_circle.base_vector[1] +
                node->gc_norm_vector[2]*bnd_circle.base_vector[2];

   int ret = 0;

   // angle + inc_angle >= M_PI_2
   if (dot <= bnd_circle.inc_angle.sin) {

      if (node->flags & U_IS_LEAF) {

         struct point_id_xyz * U = (struct point_id_xyz *)(node->U);
         size_t U_size = node->U_size;
         for (size_t i = 0; i < U_size; ++i) {
            double cos_angle =
               U[i].coordinates_xyz[0] * bnd_circle.base_vector[0] +
               U[i].coordinates_xyz[1] * bnd_circle.base_vector[1] +
               U[i].coordinates_xyz[2] * bnd_circle.base_vector[2];
            if (cos_angle > bnd_circle.inc_angle.cos) return 1;
         }

      } else {
         ret = point_check_bnd_circle(node->U, bnd_circle);
      }
   }

   // angle - inc_angle < M_PI_2
   if ((!ret) && (dot > - bnd_circle.inc_angle.sin)) {

      if (node->flags & T_IS_LEAF) {

         struct point_id_xyz * T = (struct point_id_xyz *)(node->T);
         size_t T_size = node->T_size;
         for (size_t i = 0; i < T_size; ++i) {
            double cos_angle =
               T[i].coordinates_xyz[0] * bnd_circle.base_vector[0] +
               T[i].coordinates_xyz[1] * bnd_circle.base_vector[1] +
               T[i].coordinates_xyz[2] * bnd_circle.base_vector[2];
            if (cos_angle > bnd_circle.inc_angle.cos) return 1;
         }

      } else {
         ret = point_check_bnd_circle(node->T, bnd_circle);
      }
   }

   return ret;
}

static inline int leaf_contains_matching_point(
  struct point_id_xyz * points, size_t num_points, double coordinate_xyz[3],
  size_t ** local_point_ids, size_t * local_point_ids_array_size,
  double (**result_coordinates_xyz)[3],
  size_t * result_coordinates_xyz_array_size,
  size_t total_num_local_point_ids, size_t * num_local_point_ids) {

  for (size_t i = 0; i < num_points; ++i) {

    // if the points are nearly identical
    if (points_are_identically(points[i].coordinates_xyz, coordinate_xyz)) {

      ENSURE_ARRAY_SIZE(*local_point_ids, *local_point_ids_array_size,
                        total_num_local_point_ids + 1);
      (*local_point_ids)[total_num_local_point_ids] = points[i].idx;
      if (result_coordinates_xyz != NULL) {
        ENSURE_ARRAY_SIZE(*result_coordinates_xyz,
                          *result_coordinates_xyz_array_size,
                          total_num_local_point_ids + 1);
        memcpy((*result_coordinates_xyz) + total_num_local_point_ids,
               points[i].coordinates_xyz, 3 * sizeof(double));
      }
      *num_local_point_ids = 1;
      return 1;
    }
  }

  return 0;
}

void yac_point_sphere_part_search_NN(struct point_sphere_part_search * search,
                                     size_t num_points,
                                     double (*coordinates_xyz)[3],
                                     double * cos_angles,
                                     double (**result_coordinates_xyz)[3],
                                     size_t * result_coordinates_xyz_array_size,
                                     size_t ** local_point_ids,
                                     size_t * local_point_ids_array_size,
                                     size_t * num_local_point_ids) {

  memset(num_local_point_ids, 0, num_points * sizeof(*num_local_point_ids));

  if (search == NULL) return;

  struct point_sphere_part_node * base_node = &(search->base_node);

  size_t total_num_local_point_ids = 0;

  double * dot_stack = xmalloc(search->max_tree_depth * sizeof(*dot_stack));
  struct point_sphere_part_node ** node_stack =
    xmalloc(search->max_tree_depth * sizeof(*node_stack));
  int * flags = xmalloc(search->max_tree_depth * sizeof(*flags));

  for (size_t i = 0; i < num_points; ++i) {

    struct point_sphere_part_node * curr_node = base_node;

    double * curr_coordinates_xyz = coordinates_xyz[i];

    size_t curr_tree_depth = 0;
    struct point_id_xyz * points = NULL;
    size_t curr_num_points = 0;

    // get the matching leaf for the current point
    do {

      double dot = curr_node->gc_norm_vector[0]*curr_coordinates_xyz[0] +
                   curr_node->gc_norm_vector[1]*curr_coordinates_xyz[1] +
                   curr_node->gc_norm_vector[2]*curr_coordinates_xyz[2];

      dot_stack[curr_tree_depth] = dot;
      node_stack[curr_tree_depth] = curr_node;
      flags[curr_tree_depth] = 0;

      // angle > M_PI_2
      if (dot < yac_angle_tol) {

        flags[curr_tree_depth] = U_FLAG;

        if (curr_node->flags & U_IS_LEAF) {
          if (curr_node->U_size > 0) {
            points = (struct point_id_xyz*)(curr_node->U);
            curr_num_points = curr_node->U_size;
            break;
          } else {
            flags[curr_tree_depth] = T_FLAG;
            YAC_ASSERT(
              curr_node->flags & T_IS_LEAF,
              "ERROR(yac_point_sphere_part_search_NN): "
              "if one branch is empty, the other has to be a leaf");
            points = (struct point_id_xyz*)(curr_node->T);
            curr_num_points = curr_node->T_size;
            break;
          }
        } else curr_node = curr_node->U;

      // angle < M_PI_2
      } else if (dot > -yac_angle_tol) {

        flags[curr_tree_depth] = T_FLAG;

        if (curr_node->flags & T_IS_LEAF) {
          if (curr_node->T_size > 0) {
            points = (struct point_id_xyz*)(curr_node->T);
            curr_num_points = curr_node->T_size;
            break;
          } else {
            flags[curr_tree_depth] = U_FLAG;
            YAC_ASSERT(
              curr_node->flags & U_IS_LEAF,
              "ERROR(yac_point_sphere_part_search_NN): "
              "if one branch is empty, the other has to be a leaf");
            points = (struct point_id_xyz*)(curr_node->U);
            curr_num_points = curr_node->U_size;
            break;
          }
        } else curr_node = curr_node->T;
      }

      curr_tree_depth++;
    } while (1);

    // if we do not have to do a finer search
    if (leaf_contains_matching_point(
          points, curr_num_points, curr_coordinates_xyz, local_point_ids,
          local_point_ids_array_size, result_coordinates_xyz,
          result_coordinates_xyz_array_size, total_num_local_point_ids,
          num_local_point_ids + i)) {

      if (cos_angles != NULL) cos_angles[i] = 1.0;

    } else {

      struct bounding_circle bnd_circle;
      bnd_circle.base_vector[0] = curr_coordinates_xyz[0];
      bnd_circle.base_vector[1] = curr_coordinates_xyz[1];
      bnd_circle.base_vector[2] = curr_coordinates_xyz[2];
      bnd_circle.inc_angle = SIN_COS_M_PI;
      bnd_circle.sq_crd = DBL_MAX;

      check_leaf_NN(
        points, curr_num_points, curr_coordinates_xyz, &bnd_circle.inc_angle,
        result_coordinates_xyz, result_coordinates_xyz_array_size,
        local_point_ids, local_point_ids_array_size, total_num_local_point_ids,
        num_local_point_ids + i);

      // get best result points
      point_search_NN(
        &bnd_circle, result_coordinates_xyz, result_coordinates_xyz_array_size,
        local_point_ids, local_point_ids_array_size, total_num_local_point_ids,
        num_local_point_ids + i, dot_stack, node_stack, flags, curr_tree_depth);

      if (cos_angles != NULL) cos_angles[i] = bnd_circle.inc_angle.cos;
    }

    total_num_local_point_ids += num_local_point_ids[i];
  }

  free(flags);
  free(node_stack);
  free(dot_stack);
}

static inline int compare_point_id_xyz_angle(const void * a, const void * b) {

  const struct point_id_xyz_angle * p_a = (const struct point_id_xyz_angle *)a;
  const struct point_id_xyz_angle * p_b = (const struct point_id_xyz_angle *)b;

  int ret = (p_a->cos_angle < p_b->cos_angle) -
            (p_a->cos_angle > p_b->cos_angle);

  if (ret != 0) return ret;

  return (p_a->point.idx > p_b->point.idx) - (p_a->point.idx < p_b->point.idx);
}

static size_t initial_point_bnd_search_NNN(
  size_t n, struct point_id_xyz * points, size_t num_points,
  double * point_coordinates_xyz, struct point_id_xyz_angle ** results,
  size_t * results_array_size) {

  assert(num_points > 0);

  ENSURE_ARRAY_SIZE(*results, *results_array_size, num_points);
  struct point_id_xyz_angle * results_ = *results;

#ifdef __NEC__
// vectorization of the following loop leads to failues of
//
// test_couple_config_parallel1
// test_interp_method_nnn_parallel
// test_interp_method_hcsbb_parallel
//
// with NEC compiler when CFLAGS='-O2'
#pragma _NEC novector
#endif
  for (size_t i = 0; i < num_points; ++i) {

    results_[i].point = points[i];
    results_[i].cos_angle = clamp_abs_one(
      points[i].coordinates_xyz[0] * point_coordinates_xyz[0] +
      points[i].coordinates_xyz[1] * point_coordinates_xyz[1] +
      points[i].coordinates_xyz[2] * point_coordinates_xyz[2]);
  }

  qsort(results_, num_points, sizeof(*results_), compare_point_id_xyz_angle);

  if (num_points <= n) return num_points;

  size_t num_results;
  double min_cos_angle = results_[n - 1].cos_angle;

  for (num_results = n;
       (num_results < num_points) &&
       !(fabs(min_cos_angle - results_[num_results].cos_angle) > 0.0);
       ++num_results);

  return num_results;
}

static inline struct sin_cos_angle check_leaf_NNN(
  size_t n, double * point_coordinates_xyz,
  struct point_id_xyz * points, size_t num_points,
  struct point_id_xyz_angle ** results, size_t * results_array_size,
  size_t * num_results, struct sin_cos_angle curr_angle) {

  size_t num_results_ = *num_results;
  ENSURE_ARRAY_SIZE(*results, *results_array_size, num_results_ + num_points);
  struct point_id_xyz_angle * results_ = *results;

  int flag = 0;

  double min_cos_angle = results_[num_results_-1].cos_angle;

  // check leaf for results
  for (size_t i = 0; i < num_points; ++i) {

    double curr_cos_angle =
      points[i].coordinates_xyz[0] * point_coordinates_xyz[0] +
      points[i].coordinates_xyz[1] * point_coordinates_xyz[1] +
      points[i].coordinates_xyz[2] * point_coordinates_xyz[2];

    // if the point is worse than the currently best point
    if (curr_cos_angle < min_cos_angle) continue;

    struct point_id_xyz_angle point =
      {.point = points[i], .cos_angle = curr_cos_angle};

    // insert point
    size_t j;
    for (j = 0; j < num_results_; ++j) {

      if (compare_point_id_xyz_angle(
            &point, results_ + num_results_ - j - 1) < 0) {
        results_[num_results_ - j] = results_[num_results_ - j - 1];
      } else {
        break;
      }
    }
    results_[num_results_ - j] = point;

    ++num_results_;
    flag = 1;
  }

  if (flag) {

    if (num_results_ > n) {

      size_t new_num_results;
      min_cos_angle = results_[n - 1].cos_angle;

      for (new_num_results = n;
           (new_num_results < num_results_) &&
           !(fabs(min_cos_angle - results_[new_num_results].cos_angle) > 0.0);
           ++new_num_results);
      num_results_ = new_num_results;
    }
    *num_results = num_results_;

    return
      get_vector_angle_2(
        results_[num_results_-1].point.coordinates_xyz, point_coordinates_xyz);
  } else return curr_angle;
}

static void point_search_NNN(
  size_t n, double * point_coordinates_xyz,
  struct point_id_xyz_angle ** results, size_t * results_array_size,
  size_t * num_results, double * dot_stack,
  struct point_sphere_part_node ** node_stack, int * flags,
  size_t curr_tree_depth) {

  struct sin_cos_angle angle =
    get_vector_angle_2(
      (*results)[(*num_results)-1].point.coordinates_xyz,
                       point_coordinates_xyz);

  // if we have already found at least n exactly matching points
  if ((*num_results >= n) && (angle.sin <= yac_angle_tol)) return;

  double dot = dot_stack[curr_tree_depth];
  struct point_sphere_part_node * node = node_stack[curr_tree_depth];
  int skip_U = flags[curr_tree_depth] & U_FLAG;
  int skip_T = flags[curr_tree_depth] & T_FLAG;

  do {

    if (!skip_U) {

      flags[curr_tree_depth] |= U_FLAG;

      // angle + inc_angle >= M_PI_2
      if ((dot < angle.sin) | (angle.cos <= 0.0)) {

        if (node->flags & U_IS_LEAF) {

          angle = check_leaf_NNN(
            n, point_coordinates_xyz, (struct point_id_xyz *)(node->U),
            node->U_size, results, results_array_size, num_results, angle);

        } else {

          // traverse down one level
          ++curr_tree_depth;
          node = (struct point_sphere_part_node *)(node->U);
          dot = node->gc_norm_vector[0] * point_coordinates_xyz[0] +
                node->gc_norm_vector[1] * point_coordinates_xyz[1] +
                node->gc_norm_vector[2] * point_coordinates_xyz[2];
          dot_stack[curr_tree_depth] = dot;
          node_stack[curr_tree_depth] = node;
          flags[curr_tree_depth] = 0;
          // skip_U = 0;
          skip_T = 0;
          continue;
        }
      }
    }

    if (!skip_T) {

      flags[curr_tree_depth] = U_FLAG + T_FLAG;

      // angle - inc_angle < M_PI_2
      if ((dot > - angle.sin) || (angle.cos <= 0.0)) {

        if (node->flags & T_IS_LEAF) {

          angle = check_leaf_NNN(
            n, point_coordinates_xyz, (struct point_id_xyz *)(node->T),
            node->T_size, results, results_array_size, num_results, angle);

        } else {

          // traverse down one level
          ++curr_tree_depth;
          node = (struct point_sphere_part_node *)(node->T);
          dot = node->gc_norm_vector[0] * point_coordinates_xyz[0] +
                node->gc_norm_vector[1] * point_coordinates_xyz[1] +
                node->gc_norm_vector[2] * point_coordinates_xyz[2];
          dot_stack[curr_tree_depth] = dot;
          node_stack[curr_tree_depth] = node;
          flags[curr_tree_depth] = 0;
          skip_U = 0;
          // skip_T = 0;
          continue;
        }
      }
    }

    if (curr_tree_depth == 0) break;

    // go up one level in the tree

    curr_tree_depth--;
    dot = dot_stack[curr_tree_depth];
    node = node_stack[curr_tree_depth];
    skip_U = flags[curr_tree_depth] & U_FLAG;
    skip_T = flags[curr_tree_depth] & T_FLAG;

  } while (1);
}

void yac_point_sphere_part_search_NNN(struct point_sphere_part_search * search,
                                      size_t num_points,
                                      double (*coordinates_xyz)[3], size_t n,
                                      double ** cos_angles,
                                      size_t * cos_angles_array_size,
                                      double (**result_coordinates_xyz)[3],
                                      size_t * result_coordinates_xyz_array_size,
                                      size_t ** local_point_ids,
                                      size_t * local_point_ids_array_size,
                                      size_t * num_local_point_ids) {

  if (num_points == 0) return;

  if (cos_angles != NULL)
    ENSURE_ARRAY_SIZE(*cos_angles, *cos_angles_array_size, num_points * n);

  if (n == 1) {
    yac_point_sphere_part_search_NN(
      search, num_points, coordinates_xyz, (cos_angles!=NULL)?*cos_angles:NULL,
      result_coordinates_xyz, result_coordinates_xyz_array_size,
      local_point_ids, local_point_ids_array_size, num_local_point_ids);

    size_t total_num_local_points = 0;
    for (size_t i = 0; i < num_points; ++i)
      total_num_local_points += num_local_point_ids[i];

    if ((cos_angles != NULL) && (total_num_local_points > num_points)) {

      ENSURE_ARRAY_SIZE(*cos_angles, *cos_angles_array_size,
                        total_num_local_points);

      for (size_t i = num_points - 1, offset = total_num_local_points - 1;
           i != (size_t)-1; i--) {

        for (size_t j = 0; j < num_local_point_ids[i]; ++j, --offset)
          (*cos_angles)[offset] = (*cos_angles)[i];
      }
    }
    return;
  }


  if (search == NULL) {
    memset(num_local_point_ids, 0, num_points * sizeof(*num_local_point_ids));
    return;
  }

  struct point_sphere_part_node * base_node = &(search->base_node);

  size_t total_num_local_point_ids = 0;

  double * dot_stack = xmalloc(search->max_tree_depth * sizeof(*dot_stack));
  struct point_sphere_part_node ** node_stack =
    xmalloc(search->max_tree_depth * sizeof(*node_stack));
  int * flags = xmalloc(search->max_tree_depth * sizeof(*flags));

  struct point_id_xyz_angle * results = NULL;
  size_t results_array_size = 0;

  for (size_t i = 0; i < num_points; ++i) {

    struct point_sphere_part_node * curr_node = base_node;

    double * curr_coordinates_xyz = coordinates_xyz[i];

    size_t curr_tree_depth = 0;
    struct point_id_xyz * points = search->points;
    size_t curr_num_points = 0;

    // get the matching leaf for the current point
    do {

      double dot = curr_node->gc_norm_vector[0]*curr_coordinates_xyz[0] +
                   curr_node->gc_norm_vector[1]*curr_coordinates_xyz[1] +
                   curr_node->gc_norm_vector[2]*curr_coordinates_xyz[2];

      dot_stack[curr_tree_depth] = dot;
      node_stack[curr_tree_depth] = curr_node;
      flags[curr_tree_depth] = 0;

      // angle >= M_PI_2
      if (dot <= 0.0) {

        if (curr_node->U_size < n) {

          flags[curr_tree_depth] = U_FLAG + T_FLAG;
          curr_num_points = curr_node->U_size + curr_node->T_size;
          break;
        } else if (curr_node->flags & U_IS_LEAF) {

          flags[curr_tree_depth] = U_FLAG;
          curr_num_points = curr_node->U_size;
          break;
        } else {

          flags[curr_tree_depth] = U_FLAG;
          curr_node = curr_node->U;
        }

      } else {

        if (curr_node->T_size < n) {

          flags[curr_tree_depth] = U_FLAG + T_FLAG;
          curr_num_points = curr_node->U_size + curr_node->T_size;
          break;
        } else if (curr_node->flags & T_IS_LEAF) {

          points += curr_node->U_size;
          flags[curr_tree_depth] = T_FLAG;
          curr_num_points = curr_node->T_size;
          break;
        } else {

          points += curr_node->U_size;
          flags[curr_tree_depth] = T_FLAG;
          curr_node = curr_node->T;
        }
      }

      curr_tree_depth++;
    } while (1);

    assert(curr_num_points > 0);

    size_t num_results =
      initial_point_bnd_search_NNN(
        n, points, curr_num_points, curr_coordinates_xyz,
        &results, &results_array_size);

    // do a detailed search
    point_search_NNN(
      n, curr_coordinates_xyz, &results, &results_array_size, &num_results,
      dot_stack, node_stack, flags, curr_tree_depth);

    // extract the results
    ENSURE_ARRAY_SIZE(*local_point_ids, *local_point_ids_array_size,
                      total_num_local_point_ids + num_results);
    size_t * local_point_ids_ =
      (*local_point_ids) + total_num_local_point_ids;
    double * cos_angles_;
    if (cos_angles != NULL) {
      ENSURE_ARRAY_SIZE(*cos_angles, *cos_angles_array_size,
                        total_num_local_point_ids + num_results);
      cos_angles_ = (*cos_angles) + total_num_local_point_ids;
    } else {
      cos_angles_ = NULL;
    }
    double (*result_coordinates_xyz_)[3];
    if (result_coordinates_xyz != NULL) {
      ENSURE_ARRAY_SIZE(*result_coordinates_xyz,
                        *result_coordinates_xyz_array_size,
                        total_num_local_point_ids + num_results);
      result_coordinates_xyz_ =
        (*result_coordinates_xyz) + total_num_local_point_ids;
    } else {
      result_coordinates_xyz_ = NULL;
    }

    for (size_t j = 0; j < num_results; ++j) {

      local_point_ids_[j] = results[j].point.idx;
      if (cos_angles_ != NULL) cos_angles_[j] = results[j].cos_angle;
      if (result_coordinates_xyz_ != NULL) {
        result_coordinates_xyz_[j][0] = results[j].point.coordinates_xyz[0];
        result_coordinates_xyz_[j][1] = results[j].point.coordinates_xyz[1];
        result_coordinates_xyz_[j][2] = results[j].point.coordinates_xyz[2];
      }
    }

    num_local_point_ids[i] = num_results;
    total_num_local_point_ids += num_results;
  }

  free(results);
  free(flags);
  free(node_stack);
  free(dot_stack);
}

int yac_point_sphere_part_search_bnd_circle_contains_points(
  struct point_sphere_part_search * search, struct bounding_circle circle) {

  if (search == NULL) return 0;

  return point_check_bnd_circle(&(search->base_node), circle);
}

static void search_point(struct sphere_part_node * node,
                         double point[],
                         size_t ** overlap_cells,
                         size_t * overlap_cells_array_size,
                         size_t * num_overlap_cells,
                         struct overlaps * search_interval_tree_buffer,
                         double prev_gc_norm_vector[]) {

   double dot = point[0] * node->gc_norm_vector[0] +
                point[1] * node->gc_norm_vector[1] +
                point[2] * node->gc_norm_vector[2];

   // angle < M_PI_2
   if (dot > -yac_angle_tol) {

      if (node->flags & T_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->T_size);
         memcpy(*overlap_cells + *num_overlap_cells, node->T,
                node->T_size * sizeof(**overlap_cells));
         *num_overlap_cells += node->T_size;

      } else {
         search_point(node->T, point, overlap_cells,
                      overlap_cells_array_size, num_overlap_cells,
                      search_interval_tree_buffer, node->gc_norm_vector);
      }
   }

   // angle > M_PI_2
   if (dot < yac_angle_tol) {

      if (node->flags & U_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->U_size);
         memcpy(*overlap_cells + *num_overlap_cells, node->U,
                node->U_size * sizeof(**overlap_cells));
         *num_overlap_cells += node->U_size;

      } else {
         search_point(node->U, point, overlap_cells,
                      overlap_cells_array_size, num_overlap_cells,
                      search_interval_tree_buffer, node->gc_norm_vector);
      }
   }

   // fabs(angle - M_PI_2) <= (node->I_angle)
   // fabs(cos(angle)) <= sin(node->I_angle)
   if (fabs(dot) <= node->I_angle.sin) {

      if (node->flags & I_IS_INTERVAL_TREE) {

         // project the point onto the current great circle
         double GCp[3], bVp[3];
         crossproduct_kahan(node->gc_norm_vector, point, GCp);
         crossproduct_kahan(GCp, node->gc_norm_vector, bVp);

         struct interval search_interval;
         // if the projected point does not coincide with the norm vector of
         // the plane through the great circle
         if ((fabs(bVp[0]) > 1e-9) || (fabs(bVp[1]) > 1e-9) ||
             (fabs(bVp[2]) > 1e-9)) {
            normalise_vector(bVp);
            double base_angle = get_vector_angle(bVp, prev_gc_norm_vector);
            search_interval.left=base_angle;
            search_interval.right=base_angle;
         } else {
            search_interval.left = -M_PI;
            search_interval.right = M_PI;
         }

         search_interval_tree_buffer->num_overlaps = 0;

         yac_search_interval_tree(node->I.ivt.head_node, node->I.ivt.num_nodes,
                                  search_interval, search_interval_tree_buffer);

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells +
                           search_interval_tree_buffer->num_overlaps);

         for (size_t i = 0; i < search_interval_tree_buffer->num_overlaps;
              ++i) {
            (*overlap_cells)[(*num_overlap_cells)+i] =
               node->I.ivt.head_node[
                  search_interval_tree_buffer->overlap_iv[i]].value;
         }

         *num_overlap_cells += search_interval_tree_buffer->num_overlaps;

      } else {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->I_size);
         memcpy(*overlap_cells + *num_overlap_cells, node->I.list,
                node->I_size * sizeof(**overlap_cells));
         *num_overlap_cells += node->I_size;
      }
   }
}

void yac_bnd_sphere_part_search_do_point_search(
  struct bnd_sphere_part_search * search, coordinate_pointer coordinates_xyz,
  size_t count, size_t ** cells, size_t * num_cells_per_coordinate) {

  struct sphere_part_node * base_node = &(search->base_node);

  size_t * temp_search_results = NULL;
  size_t temp_search_results_array_size = 0;
  size_t num_temp_search_results = 0;

  struct overlaps search_interval_tree_buffer = {0, 0, NULL};

  memset(
    num_cells_per_coordinate, 0, count * sizeof(*num_cells_per_coordinate));
  size_t * cells_ = NULL;
  size_t cells_array_size = 0;
  size_t num_cells = 0;
  ENSURE_ARRAY_SIZE(cells_, cells_array_size, count);

  for (size_t i = 0; i < count; ++i) {

    double * curr_coordinates_xyz = &(coordinates_xyz[i][0]);

    num_temp_search_results = 0;

    double gc_norm_vector[3] = {0.0,0.0,1.0};

    search_point(base_node, curr_coordinates_xyz, &temp_search_results,
                 &temp_search_results_array_size, &num_temp_search_results,
                 &search_interval_tree_buffer, gc_norm_vector);

    ENSURE_ARRAY_SIZE(
      cells_, cells_array_size, num_cells + num_temp_search_results);

    memcpy(cells_ + num_cells, temp_search_results,
           num_temp_search_results * sizeof(*temp_search_results));
    num_cells_per_coordinate[i] = num_temp_search_results;
    num_cells += num_temp_search_results;
  }

  free(temp_search_results);
  free(search_interval_tree_buffer.overlap_iv);

   *cells = xrealloc(cells_, num_cells * sizeof(*cells_));
}

void yac_bnd_sphere_part_search_do_bnd_circle_search(
  struct bnd_sphere_part_search * search, struct bounding_circle * bnd_circles,
  size_t count, size_t ** cells, size_t * num_cells_per_bnd_circle) {

  struct sphere_part_node * base_node = &(search->base_node);

  size_t * temp_search_results = NULL;
  size_t temp_search_results_array_size = 0;
  size_t num_temp_search_results = 0;

  struct overlaps search_interval_tree_buffer = {0, 0, NULL};

  memset(
    num_cells_per_bnd_circle, 0, count * sizeof(*num_cells_per_bnd_circle));
  size_t * cells_ = NULL;
  size_t cells_array_size = 0;
  size_t num_cells = 0;
  ENSURE_ARRAY_SIZE(cells_, cells_array_size, count);

  for (size_t i = 0; i < count; ++i) {

    num_temp_search_results = 0;

    double gc_norm_vector[3] = {0.0,0.0,1.0};

    search_bnd_circle(base_node, bnd_circles[i], &temp_search_results,
                      &temp_search_results_array_size, &num_temp_search_results,
                      &search_interval_tree_buffer, gc_norm_vector);

    ENSURE_ARRAY_SIZE(
      cells_, cells_array_size, num_cells + num_temp_search_results);

    memcpy(cells_ + num_cells, temp_search_results,
           num_temp_search_results * sizeof(*temp_search_results));
    num_cells_per_bnd_circle[i] = num_temp_search_results;
    num_cells += num_temp_search_results;
  }

  free(temp_search_results);
  free(search_interval_tree_buffer.overlap_iv);

  *cells = xrealloc(cells_, num_cells * sizeof(*cells_));
}

static void free_sphere_part_tree (struct sphere_part_node tree) {

   // free I_list
   if (tree.flags & I_IS_INTERVAL_TREE)
      free(tree.I.ivt.head_node);

   if ((tree.flags & U_IS_LEAF) == 0) {
      free_sphere_part_tree(*(struct sphere_part_node*)(tree.U));
      free(tree.U);
   }

   if ((tree.flags & T_IS_LEAF) == 0) {
      free_sphere_part_tree(*(struct sphere_part_node*)(tree.T));
      free(tree.T);
   }
}

static void free_point_sphere_part_tree (struct point_sphere_part_node * tree) {

   if ((tree->flags & U_IS_LEAF) == 0) {
      free_point_sphere_part_tree(tree->U);
      free(tree->U);
   }

   if ((tree->flags & T_IS_LEAF) == 0) {
      free_point_sphere_part_tree(tree->T);
      free(tree->T);
   }
}

void yac_delete_point_sphere_part_search(
   struct point_sphere_part_search * search) {

   if (search == NULL) return;

   free_point_sphere_part_tree(&(search->base_node));
   free(search->points);
   free(search);
}

void yac_bnd_sphere_part_search_delete(struct bnd_sphere_part_search * search) {

  if (search == NULL) return;

  free_sphere_part_tree(search->base_node);
  free(search->ids);
  free(search);
}
