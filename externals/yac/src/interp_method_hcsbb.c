/**
 * @file interp_method_hcsbb.c
 * @brief implementation of the hybrid cubic spherical Bernstein-Bezier
 *        patch interpolation method
 * (see: Xiaoyu Liu, Larry L. Schumaker, Hybrid Bézier patches on sphere-like
 *      surfaces, Journal of Computational and Applied Mathematics, Volume 73,
 *      Issues 12, 1996, Pages 157-172, ISSN 0377-0427,
 *      https://doi.org/10.1016/0377-0427(96)00041-6)
 *
 * @copyright Copyright  (C)  2020 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#include "string.h"

#include "interp_method_hcsbb.h"
#include "yac_lapack_interface.h"
#include "geometry.h"
#include "ensure_array_size.h"

// the default values for HCSBB_GAUSS_NNN, HCSBB_D_NNN, and HCSBB_GAUSS_ORDER
// were determined empirically

#ifndef HCSBB_GAUSS_NNN
#define HCSBB_GAUSS_NNN (4)
#endif
#ifndef HCSBB_D_NNN
#define HCSBB_D_NNN (7)
#endif

#ifndef HCSBB_GAUSS_ORDER
#define HCSBB_GAUSS_ORDER (6)
#endif
// static size_t const num_gauss_points[] = {1,1,3,4,6,7,12,13,16,19,25,27,33};
// -> HCSBB_NUM_GAUSS_POINTS = num_gauss_points[HCSBB_GAUSS_ORDER]
#if HCSBB_GAUSS_ORDER == 0
#define HCSBB_NUM_GAUSS_POINTS (1)
#elif HCSBB_GAUSS_ORDER == 1
#define HCSBB_NUM_GAUSS_POINTS (1)
#elif HCSBB_GAUSS_ORDER == 2
#define HCSBB_NUM_GAUSS_POINTS (3)
#elif HCSBB_GAUSS_ORDER == 3
#define HCSBB_NUM_GAUSS_POINTS (4)
#elif HCSBB_GAUSS_ORDER == 4
#define HCSBB_NUM_GAUSS_POINTS (6)
#elif HCSBB_GAUSS_ORDER == 5
#define HCSBB_NUM_GAUSS_POINTS (7)
#elif HCSBB_GAUSS_ORDER == 6
#define HCSBB_NUM_GAUSS_POINTS (12)
#elif HCSBB_GAUSS_ORDER == 7
#define HCSBB_NUM_GAUSS_POINTS (13)
#elif HCSBB_GAUSS_ORDER == 8
#define HCSBB_NUM_GAUSS_POINTS (16)
#elif HCSBB_GAUSS_ORDER == 9
#define HCSBB_NUM_GAUSS_POINTS (19)
#else
#define HCSBB_NUM_GAUSS_POINTS (-1)
#error "interpolation_method_hcsbb: the specified gauss order is currently not supported."
#endif

// number spherical Bernstein-Bezier coefficients
// n_c = = ((degree + 2) * (degree + 1)) / 2
// sbb_degree = 3
#define HCSBB_NUM_SBB_COEFFS (10)

// This tolerence needs to be lower than the
// one being used to find matching cells.
#define SB_COORD_TOL (1.0e-6)

static size_t do_search_hcsbb(struct interp_method * method,
                              struct interp_grid * interp_grid,
                              size_t * tgt_points, size_t count,
                              struct interp_weights * weights);
static void delete_hcsbb(struct interp_method * method);

static struct interp_method_vtable
  interp_method_hcsbb_vtable = {
    .do_search = do_search_hcsbb,
    .delete = delete_hcsbb};

struct interp_method_hcsbb {

  struct interp_method_vtable * vtable;
  double value;
};

//! gauss integration point patch data (used to estimate derivatives)
//! (if num_src_points == 0, no valid gauss integration point patch could be
//!  generated)
struct gauss_nnn_patch {

  double bnd_triangle[3][3]; //!< bounding triangle

  size_t * src_points;
  size_t num_src_points;

  //! after init_gauss_nnn_patch:
  //!  w_g ([HCSBB_NUM_SBB_COEFFS][HCSBB_NUM_SBB_COEFFS])
  //! after compute_gauss_nnn_patch:
  //! w_g * w_nnn ([num_src_points][HCSBB_NUM_SBB_COEFFS])
  double * weights;
};

struct weight_vector_data {
  double weight;
  size_t point_id;
};

struct weight_vector_data_pos {
  struct weight_vector_data * weight;
  size_t pos;
};

struct weight_vector {

  struct weight_vector_data * data;
  size_t n;
};

enum interp_type_flag {
  ON_VERTEX = 0,   //!< there is a source point that exactly matches the target
                   //!< point
  ON_EDGE = 1,     //!< the target point is on an edge of the source point grid
  ON_TRIANGLE = 2, //!< the target point is within a triangle of the source
                   //!< point grid
};

struct vertex_interp_data {

  //! index of the associacted source point
  size_t src_point;
  //! coordinate of the source point
  double coordinate_xyz[3];
  //! gauss integration point patch used to estimate the derivative at the
  //! source point
  struct gauss_nnn_patch * gauss_nnn_patch;
};

struct edge_interp_data {

  //! c_300 = f(v1)
  //! c_030 = f(v2)
  //! c_210 = f(v1) + D_12 p(v1)/3
  //! c_120 = f(v2) + D_21 p(v2)/3
  //! we do not explicitly store c_300 and c_030
  //! the remaining coefficients are zero on the edge
  struct weight_vector c[2]; // c_210, c_120;

  //! derivative data of both vertices
  struct vertex_interp_data * vertex_d_data[2];

  //! derivative data of edge middle point
  struct {
    //! point at which a derivative is to be estimated
    double coordinate_xyz[3];
    //! gauss integration point patch used to estimate the derivative
    struct gauss_nnn_patch * gauss_nnn_patch;
  } middle_d_data;

  //! spherical barycentric coordinate of the middle point
  //! (b0 and b1 are identical)
  double sb_coord_middle_point;

  //! vector that is a tangent of the sphere in the middle points of the edge
  //! and is perpendicular to the edge
  double g[3];
};

struct triangle_interp_data {

  //! edges of the triangle (edge[0] is the on opposing corner src_points[0])
  struct edge_interp_data * edges[3];

  //! weights for alpha values used to compute c_111
  struct weight_vector c_111_a[3];
};

struct tgt_point_search_data {

  size_t tgt_point;
  double coordinate_xyz[3];

  int is_masked;
  size_t src_points[3]; // sorted by global ids
  double sb_coords[3];

  //! depending on type flag interp_data contains a point to the required data
  enum interp_type_flag type_flag;
  union {
    struct edge_interp_data * edge;
    struct triangle_interp_data * triangle;
    size_t idx;
  } interp_data;
};

struct bnd_triangle_reorder {
  double coords[3][3]; // this need to be the first element in this struct
  size_t reorder_idx;
};

// compute the spherical barycentric coordinates of the given vertices with
// respect to the given triangle
static void compute_sb_coords(
  double * sb_coords, size_t num_vertices, double triangle[3][3]) {

  double A[3][3];
  lapack_int n = 3, nrhs = (lapack_int) num_vertices, lda = n, ldx = n, ipiv[3];
  memcpy(A, triangle, sizeof(A));

  // for a vertex v the spherical barycentric coordinates b are defined as
  // follows
  // A * b = v
  // where: A is the matrix consisting of the vertex coordinates of the three
  // corners of the triangle
  // we compute b by solving this linear system using LAPACK
  YAC_ASSERT(
    !LAPACKE_dgesv(
      LAPACK_COL_MAJOR, n, nrhs, &A[0][0], lda, ipiv, sb_coords, ldx),
    "ERROR: internal error (could not solve linear 3x3 system)")
}

static inline int compare_size_t(const void * a, const void * b) {

  size_t const * a_ = a, * b_ = b;

  return (*a_ > *b_) - (*b_ > *a_);
}

static inline int compare_2_size_t(const void * a, const void * b) {

  size_t const * a_ = a, * b_ = b;
  int ret;
  if ((ret = (a_[0] > b_[0]) - (b_[0] > a_[0]))) return ret;
  else return (a_[1] > b_[1]) - (b_[1] > a_[1]);
}

static inline int compare_3_size_t(const void * a, const void * b) {

  size_t const * a_ = a, * b_ = b;
  int ret;
  if ((ret = (a_[0] > b_[0]) - (b_[0] > a_[0]))) return ret;
  if ((ret = (a_[1] > b_[1]) - (b_[1] > a_[1]))) return ret;
  else return (a_[2] > b_[2]) - (b_[2] > a_[2]);
}

static inline int compare_HCSBB_D_NNN_size_t(const void * a, const void * b) {

  size_t const * a_ = a, * b_ = b;
  int ret;
  for (size_t i = 0; i < HCSBB_D_NNN; ++i)
    if ((ret = (a_[i] > b_[i]) - (b_[i] > a_[i])))
      return ret;
  return 0;
}

// Go through the source field triangles for all target points and check for
// duplicated vertices, edges, and triangles. Some computation for each
// require source field vertex, edge, and triangle can be done independent of
// the actual target point, and hence needs to be done only once.
static void get_unique_interp_data(
  struct tgt_point_search_data * tgt_point_data, size_t num_edges,
  size_t num_triangles,
  size_t **unique_vertices, size_t * num_unique_vertices,
  size_t (**unique_edge_to_unique_vertex)[2], size_t * num_unique_edges,
  size_t (**unique_triangle_to_unique_edge)[3], size_t * num_unique_triangles) {

  size_t (*edge_vertices)[2] =
    xmalloc((num_edges + 3 * num_triangles) * sizeof(*edge_vertices));
  size_t (*triangle_vertices)[3] =
    xmalloc(num_triangles * sizeof(*triangle_vertices));

  // get all edges and triangles
  for (size_t i = 0; i < num_edges; ++i, ++tgt_point_data)
    memcpy(edge_vertices[i], tgt_point_data->src_points, 2 * sizeof(size_t));
  for (size_t i = 0; i < num_triangles; ++i, ++tgt_point_data) {
    edge_vertices[3*i+num_edges+0][0] = tgt_point_data->src_points[0];
    edge_vertices[3*i+num_edges+0][1] = tgt_point_data->src_points[1];
    edge_vertices[3*i+num_edges+1][0] = tgt_point_data->src_points[0];
    edge_vertices[3*i+num_edges+1][1] = tgt_point_data->src_points[2];
    edge_vertices[3*i+num_edges+2][0] = tgt_point_data->src_points[1];
    edge_vertices[3*i+num_edges+2][1] = tgt_point_data->src_points[2];
    memcpy(
      triangle_vertices[i], tgt_point_data->src_points, 3 * sizeof(size_t));
  }

  // sort edge end triangle data
  qsort(edge_vertices, num_edges + 3 * num_triangles,
        sizeof(*edge_vertices), compare_2_size_t);
  qsort(triangle_vertices, num_triangles,
        sizeof(*triangle_vertices), compare_3_size_t);

  *num_unique_edges = num_edges + 3 * num_triangles;
  *num_unique_triangles = num_triangles;

  // remove duplicated edges and vertices
  yac_remove_duplicates_size_t_2(edge_vertices, num_unique_edges);
  yac_remove_duplicates_size_t_3(triangle_vertices, num_unique_triangles);

  // get all unique vertices
  *num_unique_vertices = 2 * (*num_unique_edges);
  size_t * vertices = xmalloc((*num_unique_vertices) * sizeof(*vertices));
  memcpy(
    vertices, &(edge_vertices[0][0]), *num_unique_vertices * sizeof(*vertices));
  qsort(vertices, *num_unique_vertices, sizeof(*vertices), compare_size_t);
  yac_remove_duplicates_size_t(vertices, num_unique_vertices);
  *unique_vertices =
    xrealloc(vertices, *num_unique_vertices * sizeof(*vertices));

  size_t (*triangle_edges)[3] =
    xmalloc(3 * (*num_unique_triangles) * sizeof(*triangle_edges));
  for (size_t i = 0; i < *num_unique_triangles; ++i) {
    // edge across first vertex
    triangle_edges[3*i+0][0] = triangle_vertices[i][1];
    triangle_edges[3*i+0][1] = triangle_vertices[i][2];
    triangle_edges[3*i+0][2] = 3*i+0;
    // edge across second vertex
    triangle_edges[3*i+1][0] = triangle_vertices[i][0];
    triangle_edges[3*i+1][1] = triangle_vertices[i][2];
    triangle_edges[3*i+1][2] = 3*i+1;
    // edge across third vertex
    triangle_edges[3*i+2][0] = triangle_vertices[i][0];
    triangle_edges[3*i+2][1] = triangle_vertices[i][1];
    triangle_edges[3*i+2][2] = 3*i+2;
  }

  // sort triangle edges
  qsort(triangle_edges, 3 * (*num_unique_triangles),
        sizeof(*triangle_edges), compare_2_size_t);

  // look up triangle edges in the unique edges
  *unique_triangle_to_unique_edge =
    xmalloc(*num_unique_triangles * sizeof(**unique_triangle_to_unique_edge));
  for (size_t i = 0, j = 0; i < 3 * (*num_unique_triangles); ++i) {
    size_t * curr_edge = triangle_edges[i];
    while (compare_2_size_t(curr_edge, edge_vertices[j])) ++j;
    (&(*unique_triangle_to_unique_edge)[0][0])[curr_edge[2]] = j;
  }
  free(triangle_edges);

  tgt_point_data -= num_edges + num_triangles;
  for (size_t i = 0, j = 0; i < num_edges; ++i, ++tgt_point_data) {
    while (compare_2_size_t(edge_vertices[j], tgt_point_data->src_points)) ++j;
    YAC_ASSERT(
      j < *num_unique_edges, "ERROR(get_unique_interp_data): edge ids missmatch")
    tgt_point_data->interp_data.idx = j;
  }
  for (size_t i = 0, j = 0; i < num_triangles; ++i, ++tgt_point_data) {
    while (compare_3_size_t(triangle_vertices[j], tgt_point_data->src_points))
      ++j;
    YAC_ASSERT(
      j < *num_unique_triangles,
      "ERROR(get_unique_interp_data): triangle ids missmatch")
    tgt_point_data->interp_data.idx = j;
  }
  free(triangle_vertices);

  // generate unique_edge_to_unique_vertex
  *unique_edge_to_unique_vertex =
    xmalloc(*num_unique_edges * sizeof(**unique_edge_to_unique_vertex));
  size_t * reorder_idx =
    xmalloc(2 * (*num_unique_edges) * sizeof(*reorder_idx));
  for (size_t i = 0; i < 2 * (*num_unique_edges); ++i) reorder_idx[i] = i;
  yac_quicksort_index_size_t_size_t(
    &(edge_vertices[0][0]), 2 * (*num_unique_edges), reorder_idx);
  for (size_t i = 0, j = 0; i < 2 * (*num_unique_edges); ++i) {
    size_t curr_vertex = (&(edge_vertices[0][0]))[i];
    while (curr_vertex != (*unique_vertices)[j]) ++j;
    (&(*unique_edge_to_unique_vertex[0][0]))[reorder_idx[i]] = j;
  }
  free(edge_vertices);
  free(reorder_idx);
}

static void compute_gauss_nnn_patch(
  struct gauss_nnn_patch * gauss_nnn_patch,
  double (*gauss_points)[3], size_t * nnn_search_results,
  double w_nnn[][HCSBB_NUM_GAUSS_POINTS], size_t * src_point_buffer,
  yac_int * global_id_buffer, const_coordinate_pointer src_point_coords,
  const_yac_int_pointer src_point_global_ids) {

  // get all unique nnn_search result points
  memcpy(src_point_buffer, nnn_search_results,
         HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS * sizeof(*src_point_buffer));
  qsort(src_point_buffer, HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS,
        sizeof(*src_point_buffer), compare_size_t);
  size_t num_unique_result_points = HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS;
  yac_remove_duplicates_size_t(src_point_buffer, &num_unique_result_points);

  // sort result points by their global ids
  for (size_t i = 0; i < num_unique_result_points; ++i)
    global_id_buffer[i] = src_point_global_ids[src_point_buffer[i]];
  yac_quicksort_index_yac_int_size_t(
    global_id_buffer, num_unique_result_points, src_point_buffer);

  gauss_nnn_patch->src_points =
    xmalloc(num_unique_result_points * sizeof(*src_point_buffer));
  memcpy(gauss_nnn_patch->src_points, src_point_buffer,
         num_unique_result_points * sizeof(*src_point_buffer));
  gauss_nnn_patch->num_src_points = num_unique_result_points;

  // compute w_nnn
  {
    for (size_t i = 0; i < num_unique_result_points; ++i)
      for (size_t j = 0; j < HCSBB_NUM_GAUSS_POINTS; ++j)
        w_nnn[i][j] = 0.0;

    double inv_distances[HCSBB_GAUSS_NNN];
    size_t reorder_idx[HCSBB_GAUSS_NNN];
    for (size_t i = 0; i < HCSBB_NUM_GAUSS_POINTS; ++i) {

      size_t * curr_result_points = nnn_search_results + i * HCSBB_GAUSS_NNN;
      double * curr_gauss_points = gauss_points[i];
      double distance_sum = 0.0;
      size_t j;

      // for all results of the current request point
      for (j = 0; j < HCSBB_GAUSS_NNN; ++j) {

        double distance =
          get_vector_angle(
            curr_gauss_points,
            (double*)(src_point_coords[curr_result_points[j]]));

        // if src and tgt point are at the same location
        if (distance < yac_angle_tol) break;

        double inv_distance = 1.0 / (distance * distance);
        inv_distances[j] = inv_distance;
        distance_sum += inv_distance;
      }

      // if there is no source point that exactly matches the target point
      if (j == HCSBB_GAUSS_NNN) {

        yac_int curr_global_ids[HCSBB_GAUSS_NNN];
        for (size_t k = 0; k < HCSBB_GAUSS_NNN; ++k)
          curr_global_ids[k] = src_point_global_ids[curr_result_points[k]];
        for (size_t l = 0; l < HCSBB_GAUSS_NNN; ++l) reorder_idx[l] = l;

        // sort current results by their global ids
        yac_quicksort_index_yac_int_size_t(
          curr_global_ids, HCSBB_GAUSS_NNN, reorder_idx);

        distance_sum = 1.0 / distance_sum;
        for (size_t k = 0, l = 0; k < HCSBB_GAUSS_NNN; ++k) {
          yac_int curr_global_id = curr_global_ids[k];
          while (curr_global_id != global_id_buffer[l]) ++l;
          w_nnn[l][i] = inv_distances[reorder_idx[k]] * distance_sum;
        }

      } else {

        yac_int curr_global_id = src_point_global_ids[curr_result_points[j]];
        size_t l = 0;
        while (curr_global_id != global_id_buffer[l]) ++l;
        w_nnn[l][i] = 1.0;
      }
    }
  }

  double (*w_g)[HCSBB_NUM_GAUSS_POINTS] =
    (double (*)[HCSBB_NUM_GAUSS_POINTS])(gauss_nnn_patch->weights);

  // compute and store w_g * w_nnn
  double (*weights)[HCSBB_NUM_SBB_COEFFS] =
    xmalloc(num_unique_result_points * sizeof(*weights));
  gauss_nnn_patch->weights = (double*)weights;
  for (size_t i = 0; i < num_unique_result_points; ++i) {
    for (size_t j = 0; j < HCSBB_NUM_SBB_COEFFS; ++j) {
      double accu = 0.0;
      for (size_t l = 0; l < HCSBB_NUM_GAUSS_POINTS; ++l)
        accu += w_g[j][l] * w_nnn[i][l];
      weights[i][j] = accu;
    }
  }

  free(w_g);
}

// computes the directional derivatives of the spherical Bernstein basis
// polynomials for degree = 3 in the given direction
// the polynomials B_ijk are ordered in the following order:
// B_003, B_012, B_021, B_030, B_102, B_111, B_120, B_201, B_210, B_300
static void compute_d_sbb_polynomials_3d(
  double bnd_triangle[3][3], double direction[3], double coordinate_xyz[3],
  double d_sbb_polynomials[HCSBB_NUM_SBB_COEFFS]) {

  double sb_coord_p[3] = {
      coordinate_xyz[0], coordinate_xyz[1], coordinate_xyz[2]};
  double sb_coord_d[3] = {direction[0], direction[1], direction[2]};

  // compute the spherical barycentric coordinates of the point and the
  // direction
  compute_sb_coords(&(sb_coord_p[0]), 1, bnd_triangle);
  compute_sb_coords(&(sb_coord_d[0]), 1, bnd_triangle);

#define POW_0(x) (1.0)
#define POW_1(x) ((x))
#define POW_2(x) ((x)*(x))
#define POW_3(x) ((x)*(x)*(x))
#define FAC_0 (1.0)
#define FAC_1 (1.0)
#define FAC_2 (2.0)
#define FAC_3 (6.0)

  d_sbb_polynomials[0] =
    // sb_coord_d[0] * (0.0 * FAC_3 / (FAC_0 * FAC_0 * FAC_3)) *
    // (POW_0(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_3(sb_coord_p[2])) +
    // sb_coord_d[1] * (0.0 * FAC_3 / (FAC_0 * FAC_0 * FAC_3)) *
    // (POW_0(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_3(sb_coord_p[2])) +
    sb_coord_d[2] * (3.0 * FAC_3 / (FAC_0 * FAC_0 * FAC_3)) *
    (POW_0(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_2(sb_coord_p[2]));
  d_sbb_polynomials[1] =
    // sb_coord_d[0] * (0.0 * FAC_3 / (FAC_0 * FAC_1 * FAC_2)) *
    // (POW_0(sb_coord_p[0]) * POW_1(sb_coord_p[1]) * POW_2(sb_coord_p[2])) +
    sb_coord_d[1] * (1.0 * FAC_3 / (FAC_0 * FAC_1 * FAC_2)) *
    (POW_0(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_2(sb_coord_p[2])) +
    sb_coord_d[2] * (2.0 * FAC_3 / (FAC_0 * FAC_1 * FAC_2)) *
    (POW_0(sb_coord_p[0]) * POW_1(sb_coord_p[1]) * POW_1(sb_coord_p[2]));
  d_sbb_polynomials[2] =
    // sb_coord_d[0] * (0.0 * FAC_3 / (FAC_0 * FAC_2 * FAC_1)) *
    // (POW_0(sb_coord_p[0]) * POW_2(sb_coord_p[1]) * POW_1(sb_coord_p[2])) +
    sb_coord_d[1] * (2.0 * FAC_3 / (FAC_0 * FAC_2 * FAC_1)) *
    (POW_0(sb_coord_p[0]) * POW_1(sb_coord_p[1]) * POW_1(sb_coord_p[2])) +
    sb_coord_d[2] * (1.0 * FAC_3 / (FAC_0 * FAC_2 * FAC_1)) *
    (POW_0(sb_coord_p[0]) * POW_2(sb_coord_p[1]) * POW_0(sb_coord_p[2]));
  d_sbb_polynomials[3] =
    // sb_coord_d[0] * (0.0 * FAC_3 / (FAC_0 * FAC_3 * FAC_0)) *
    // (POW_0(sb_coord_p[0]) * POW_3(sb_coord_p[1]) * POW_0(sb_coord_p[2])) +
    sb_coord_d[1] * (3.0 * FAC_3 / (FAC_0 * FAC_3 * FAC_0)) *
    (POW_0(sb_coord_p[0]) * POW_2(sb_coord_p[1]) * POW_0(sb_coord_p[2]));
    // sb_coord_d[2] * (0.0 * FAC_3 / (FAC_0 * FAC_3 * FAC_0)) *
    // (POW_0(sb_coord_p[0]) * POW_3(sb_coord_p[1]) * POW_0(sb_coord_p[2]));
  d_sbb_polynomials[4] =
    sb_coord_d[0] * (1.0 * FAC_3 / (FAC_1 * FAC_0 * FAC_2)) *
    (POW_0(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_2(sb_coord_p[2])) +
    // sb_coord_d[1] * (0.0 * FAC_3 / (FAC_1 * FAC_0 * FAC_2)) *
    // (POW_1(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_2(sb_coord_p[2])) +
    sb_coord_d[2] * (2.0 * FAC_3 / (FAC_1 * FAC_0 * FAC_2)) *
    (POW_1(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_1(sb_coord_p[2]));
  d_sbb_polynomials[5] =
    sb_coord_d[0] * (1.0 * FAC_3 / (FAC_1 * FAC_1 * FAC_1)) *
    (POW_0(sb_coord_p[0]) * POW_1(sb_coord_p[1]) * POW_1(sb_coord_p[2])) +
    sb_coord_d[1] * (1.0 * FAC_3 / (FAC_1 * FAC_1 * FAC_1)) *
    (POW_1(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_1(sb_coord_p[2])) +
    sb_coord_d[2] * (1.0 * FAC_3 / (FAC_1 * FAC_1 * FAC_1)) *
    (POW_1(sb_coord_p[0]) * POW_1(sb_coord_p[1]) * POW_0(sb_coord_p[2]));
  d_sbb_polynomials[6] =
    sb_coord_d[0] * (1.0 * FAC_3 / (FAC_1 * FAC_2 * FAC_0)) *
    (POW_0(sb_coord_p[0]) * POW_2(sb_coord_p[1]) * POW_0(sb_coord_p[2])) +
    sb_coord_d[1] * (2.0 * FAC_3 / (FAC_1 * FAC_2 * FAC_0)) *
    (POW_1(sb_coord_p[0]) * POW_1(sb_coord_p[1]) * POW_0(sb_coord_p[2]));
    // sb_coord_d[2] * (0.0 * FAC_3 / (FAC_1 * FAC_2 * FAC_0)) *
    // (POW_1(sb_coord_p[0]) * POW_2(sb_coord_p[1]) * POW_0(sb_coord_p[2]));
  d_sbb_polynomials[7] =
    sb_coord_d[0] * (2.0 * FAC_3 / (FAC_2 * FAC_0 * FAC_1)) *
    (POW_1(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_1(sb_coord_p[2])) +
    // sb_coord_d[1] * (0.0 * FAC_3 / (FAC_2 * FAC_0 * FAC_1)) *
    // (POW_2(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_1(sb_coord_p[2])) +
    sb_coord_d[2] * (1.0 * FAC_3 / (FAC_2 * FAC_0 * FAC_1)) *
    (POW_2(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_0(sb_coord_p[2]));
  d_sbb_polynomials[8] =
    sb_coord_d[0] * (2.0 * FAC_3 / (FAC_2 * FAC_1 * FAC_0)) *
    (POW_1(sb_coord_p[0]) * POW_1(sb_coord_p[1]) * POW_0(sb_coord_p[2])) +
    sb_coord_d[1] * (1.0 * FAC_3 / (FAC_2 * FAC_1 * FAC_0)) *
    (POW_2(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_0(sb_coord_p[2]));
    // sb_coord_d[2] * (0.0 * FAC_3 / (FAC_2 * FAC_1 * FAC_0)) *
    // (POW_2(sb_coord_p[0]) * POW_1(sb_coord_p[1]) * POW_0(sb_coord_p[2]));
  d_sbb_polynomials[9] =
    sb_coord_d[0] * (3.0 * FAC_3 / (FAC_3 * FAC_0 * FAC_0)) *
    (POW_2(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_0(sb_coord_p[2]));
    // sb_coord_d[1] * (0.0 * FAC_3 / (FAC_3 * FAC_0 * FAC_0)) *
    // (POW_3(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_0(sb_coord_p[2])) +
    // sb_coord_d[2] * (0.0 * FAC_3 / (FAC_3 * FAC_0 * FAC_0)) *
    // (POW_3(sb_coord_p[0]) * POW_0(sb_coord_p[1]) * POW_0(sb_coord_p[2]));

#undef POW_0
#undef POW_1
#undef POW_2
#undef POW_3
#undef FAC_0
#undef FAC_1
#undef FAC_2
#undef FAC_3
}

static void get_derivative_weights(
  struct gauss_nnn_patch * gauss_nnn_patch, double coordinate_xyz[3],
  double direction[3], struct weight_vector_data * weights,
  double mult_factor) {

  // compute w_p
  double d_sbb_vertex_polynomials[HCSBB_NUM_SBB_COEFFS];
  compute_d_sbb_polynomials_3d(
    gauss_nnn_patch->bnd_triangle, direction, coordinate_xyz,
    d_sbb_vertex_polynomials);

  // compute w_p * (w_g * w_nnn)
  size_t * src_points = gauss_nnn_patch->src_points;
  size_t num_src_points = gauss_nnn_patch->num_src_points;
  double * d_data_weights = gauss_nnn_patch->weights;
  for (size_t i = 0, k = 0; i < num_src_points; ++i) {

    double weight = 0.0;
    for (size_t j = 0; j < HCSBB_NUM_SBB_COEFFS; ++j, ++k)
      weight += d_sbb_vertex_polynomials[j] * d_data_weights[k];

    weights[i].weight = weight * mult_factor;
    weights[i].point_id = src_points[i];
  }
}

static void compute_edge_coefficients(
  struct edge_interp_data * edge_data, size_t num_edges) {

  size_t edge_weight_vector_data_buffer_size = 0;
  // count total number of weight_vector_data elements required
  for (size_t i = 0; i < num_edges; ++i)
    edge_weight_vector_data_buffer_size +=
      2 + edge_data[i].vertex_d_data[0]->gauss_nnn_patch->num_src_points +
      edge_data[i].vertex_d_data[1]->gauss_nnn_patch->num_src_points;

  // 2 -> two corner weights
  // 2*(HCSBB_NUM_GAUSS_POINTS*HCSBB_GAUSS_NNN+1)
  //   -> maximum number of weights per edge coefficient
  struct weight_vector_data * curr_edge_weight_vector_data =
    (edge_weight_vector_data_buffer_size > 0)?
      xmalloc(edge_weight_vector_data_buffer_size *
              sizeof(*curr_edge_weight_vector_data)):NULL;

  for (size_t i = 0; i < num_edges; ++i) {

    struct edge_interp_data * curr_edge = edge_data + i;
    struct vertex_interp_data * curr_vertex_d_data[2] = {
      curr_edge->vertex_d_data[0], curr_edge->vertex_d_data[1]};

    double direction[3] = {
      curr_vertex_d_data[1]->coordinate_xyz[0] -
      curr_vertex_d_data[0]->coordinate_xyz[0],
      curr_vertex_d_data[1]->coordinate_xyz[1] -
      curr_vertex_d_data[0]->coordinate_xyz[1],
      curr_vertex_d_data[1]->coordinate_xyz[2] -
      curr_vertex_d_data[0]->coordinate_xyz[2]};

    // c_210 and c_120
    for (size_t j = 0; j < 2; ++j) {

      size_t n = curr_vertex_d_data[j]->gauss_nnn_patch->num_src_points + 1;

      get_derivative_weights(
        curr_vertex_d_data[j]->gauss_nnn_patch,
        curr_vertex_d_data[j]->coordinate_xyz, direction,
        curr_edge_weight_vector_data, 1.0/3.0);

      curr_edge_weight_vector_data[n-1].weight = 1.0;
      curr_edge_weight_vector_data[n-1].point_id =
        curr_edge->vertex_d_data[j]->src_point;
      curr_edge->c[j].data = curr_edge_weight_vector_data;
      curr_edge->c[j].n = n;

      direction[0] *= -1, direction[1] *= -1, direction[2] *= -1;

      curr_edge_weight_vector_data += n;
    }
  }
}

static int compare_weight_vector_data_pos_weight(
  void const * a, void const * b) {

  struct weight_vector_data const * weight_a =
    ((struct weight_vector_data_pos const *)a)->weight;
  struct weight_vector_data const * weight_b =
    ((struct weight_vector_data_pos const *)b)->weight;

  int ret = (weight_a->point_id > weight_b->point_id) -
            (weight_a->point_id < weight_b->point_id);
  if (ret) return ret;
  double abs_weight_a = fabs(weight_a->weight);
  double abs_weight_b = fabs(weight_b->weight);
  ret = (abs_weight_a > abs_weight_b) - (abs_weight_a < abs_weight_b);
  if (ret) return ret;
  return (weight_a->weight > weight_b->weight) -
         (weight_a->weight < weight_b->weight);
}

static int compare_weight_vector_data_point(
  void  const * a, void const * b) {

  struct weight_vector_data const * weight_a =
    (struct weight_vector_data const *)a;
  struct weight_vector_data const * weight_b =
    (struct weight_vector_data const *)b;

  return (weight_a->point_id > weight_b->point_id) -
         (weight_a->point_id < weight_b->point_id);
}

static int compare_weight_vector_data_pos_pos(
  void const * a, void const * b) {

  return ((struct weight_vector_data_pos const *)a)->pos -
         ((struct weight_vector_data_pos const *)b)->pos;
}

static void compact_weight_vector_data(
  struct weight_vector_data * weights, size_t * n,
  struct weight_vector_data_pos * buffer) {

  size_t n_ = *n;

  if (n_ <= 1) return;

  for (size_t i = 0; i < n_; ++i) {
    buffer[i].weight = weights + i;
    buffer[i].pos = i;
  }

  qsort(buffer, n_, sizeof(*buffer), compare_weight_vector_data_pos_weight);

  size_t new_n = 1;
  struct weight_vector_data_pos * prev_weight_data_pos = buffer;
  struct weight_vector_data * prev_weight_data = prev_weight_data_pos->weight;
  struct weight_vector_data_pos * curr_weight_data_pos = buffer + 1;
  for (size_t i = 1; i < n_; ++i, ++curr_weight_data_pos) {

    struct weight_vector_data * curr_weight_data = curr_weight_data_pos->weight;
    if (!compare_weight_vector_data_point(prev_weight_data, curr_weight_data)) {
      prev_weight_data->weight += curr_weight_data->weight;
      if (prev_weight_data_pos->pos > curr_weight_data_pos->pos)
        prev_weight_data_pos->pos = curr_weight_data_pos->pos;
    } else {
      ++new_n;
      ++prev_weight_data_pos;
      *prev_weight_data_pos = *curr_weight_data_pos;
      prev_weight_data = prev_weight_data_pos->weight;
    }
  }

  if (n_ == new_n) return;

  qsort(buffer, new_n, sizeof(*buffer), compare_weight_vector_data_pos_pos);

  for (size_t i = 0; i < new_n; ++i) weights[i] = *(buffer[i].weight);

  *n = new_n;
}

static void compute_triangle_coefficients(
  struct triangle_interp_data * triangle_data, size_t num_triangles) {

  // maximum number of coefficients per c_111 alpha value
  // 6 edge coefficients (HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS + 1)
  // 2 corner coefficients (1)
  // 1 middle point coefficient (HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS)
  struct weight_vector_data * weight_vector_data_buffer =
    xmalloc((6 * (HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS + 1) + 2 * 1 +
             HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS) *
             sizeof(*weight_vector_data_buffer));
  struct weight_vector_data_pos * weight_vector_data_pos_buffer =
    xmalloc((6 * (HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS + 1) + 2 * 1 +
             HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS) *
             sizeof(*weight_vector_data_pos_buffer));

  struct triangle_interp_data * curr_triangle = triangle_data;
  for (size_t i = 0; i < num_triangles; ++i, ++curr_triangle) {

    // z_w0/3 = B_020(w0) * (b_0(g0)*c_120 + b_1(g0)*c_030 + b_2(g0)*c_021) +
    //          B_011(w0) * (b_0(g0)*a_1     b_1(g0)*c_021 + b_2(g0)*c_012) +
    //          B_002(w0) * (b_0(g0)*c_102 + b_1(g0)*c_012 + b_2(g0)*c_003)
    // w0 = ||v_1+v_2||  -> middle point of edge 0
    // g0 = v_2 x v_1    -> vector, which is a tangent of the sphere in w0 and
    //                      perpendicular to w0
    // z_w0 = D_g0 f(w0) -> directional derivative of f in w0 in the direction
    //                      of g0
    // =>
    // a_0 = (D_g f(w0)/3 -
    //        B_020(w0) * (b_0(g0)*c_120 + b_1(g0)*c_030 + b_2(g0)*c_021) -
    //        B_011(w0) * (                b_1(g0)*c_021 + b_2(g0)*c_012) -
    //        B_002(w0) * (b_0(g0)*c_102 + b_1(g0)*c_012 + b_2(g0)*c_003)) /
    //       (B_011(w0) * b_0(g0))
    // a_1 = (D_g f(w1)/3 -
    //        B_200(w1) * (b_0(g1)*c_300 + b_1(g1)*c_210 + b_2(g1)*c_201) -
    //        B_101(w1) * (b_0(g1)*c_201 +                 b_2(g1)*c_102) -
    //        B_002(w1) * (b_0(g1)*c_102 + b_1(g1)*c_012 + b_2(g1)*c_003)) /
    //       (B_101(w1) * b_1(g1))
    // a_2 = (D_g f(w2)/3 -
    //        B_200(w2) * (b_0(g2)*c_300 + b_1(g2)*c_210 + b_2(g2)*c_201) -
    //        B_110(w2) * (b_0(g2)*c_210 + b_1(g2)*c_120                ) -
    //        B_020(w2) * (b_0(g2)*c_120 + b_1(g2)*c_030 + b_2(g2)*c_021)) /
    //       (B_110(w2) * b_2(g2))

    struct edge_interp_data ** curr_edges = curr_triangle->edges;
    struct vertex_interp_data * vertex_d_data[3] = {
      curr_edges[2]->vertex_d_data[0],
      curr_edges[2]->vertex_d_data[1],
      curr_edges[1]->vertex_d_data[1]};
    size_t src_points[3] = {
      vertex_d_data[0]->src_point,
      vertex_d_data[1]->src_point,
      vertex_d_data[2]->src_point};

    // remark the corners of the triangle are sorted by global ids, the same is
    // true for the edge data
    struct weight_vector_data c_3_data[3] = {
      {.weight = 1.0, .point_id = src_points[0]},
      {.weight = 1.0, .point_id = src_points[1]},
      {.weight = 1.0, .point_id = src_points[2]}};
    struct weight_vector c_3[3] = {{.data = &(c_3_data[0]), .n = 1},
                                   {.data = &(c_3_data[1]), .n = 1},
                                   {.data = &(c_3_data[2]), .n = 1}};
    struct weight_vector * c_300 = &(c_3[0]);
    struct weight_vector * c_030 = &(c_3[1]);
    struct weight_vector * c_003 = &(c_3[2]);
    struct weight_vector * c_021 = &(curr_edges[0]->c[0]);
    struct weight_vector * c_012 = &(curr_edges[0]->c[1]);
    struct weight_vector * c_201 = &(curr_edges[1]->c[0]);
    struct weight_vector * c_102 = &(curr_edges[1]->c[1]);
    struct weight_vector * c_210 = &(curr_edges[2]->c[0]);
    struct weight_vector * c_120 = &(curr_edges[2]->c[1]);

    double corner_coordinate_xyz[3][3] = {
      {vertex_d_data[0]->coordinate_xyz[0],
       vertex_d_data[0]->coordinate_xyz[1],
       vertex_d_data[0]->coordinate_xyz[2]},
      {vertex_d_data[1]->coordinate_xyz[0],
       vertex_d_data[1]->coordinate_xyz[1],
       vertex_d_data[1]->coordinate_xyz[2]},
      {vertex_d_data[2]->coordinate_xyz[0],
       vertex_d_data[2]->coordinate_xyz[1],
       vertex_d_data[2]->coordinate_xyz[2]}};

    // vector which are tangents of the sphere in the middle points of the edges
    // and are perpendicular to the edges
    double * g[3] = {
      &(curr_edges[0]->g[0]), &(curr_edges[1]->g[0]), &(curr_edges[2]->g[0])};

    // spherical barycentric coordinates of g
    double b_g[3][3] = {{g[0][0], g[0][1], g[0][2]},
                        {g[1][0], g[1][1], g[1][2]},
                        {g[2][0], g[2][1], g[2][2]}};
    compute_sb_coords(&(b_g[0][0]), 3, corner_coordinate_xyz);

    // computing quadratic spherical Bernstein-Bezier polynomials for the edge
    // middle points (e.g. B_020(w0)) using w0 as an example
    //   I. w0 is in the middle of the edge opposing corner 0
    //  II. for the barycentric coordinates if w0 we can say the following
    //    a) b_0(w0) = 0 -> because it is on the edge opposing v0
    //    b) b_1(w0) = b_2(w0) -> because it is in the middle between v1 and v2
    // III. B_ijk(p) = d!/(i!*j!*k!) * b_0(p)^i * b_1(p)^j * b_2(p)^k
    //      -> for all (i + j + k) = d
    //      -> p = w0 and d = 2 in our case
    //
    // from II.a) it follows: B_200 = B_110 = B_101 = 0
    // from II.b) it follows: B_020 = B_002 = 2 * B_011 = b_1(w0) = b_2(w0) = bw

    double bw[3] = {
      curr_edges[0]->sb_coord_middle_point *
      curr_edges[0]->sb_coord_middle_point,
      curr_edges[1]->sb_coord_middle_point *
      curr_edges[1]->sb_coord_middle_point,
      curr_edges[2]->sb_coord_middle_point *
      curr_edges[2]->sb_coord_middle_point};

#define COPY_D_DATA(edge_idx) \
  { \
    get_derivative_weights( \
      curr_edges[edge_idx]->middle_d_data.gauss_nnn_patch, \
      curr_edges[edge_idx]->middle_d_data.coordinate_xyz, g[edge_idx], \
      weight_vector_data_buffer + n, 1.0 / 3.0); \
    n += curr_edges[edge_idx]->middle_d_data.gauss_nnn_patch->num_src_points; \
  }
#define COPY_COEFF(c_ijk, edge_idx, b_idx, bw_factor) \
  { \
    struct weight_vector_data * curr_buffer = weight_vector_data_buffer + n; \
    size_t curr_n = c_ijk->n; \
    memcpy(curr_buffer, c_ijk->data, curr_n * sizeof(*curr_buffer)); \
    n += curr_n; \
    double factor = - bw[edge_idx] * b_g[edge_idx][b_idx] * bw_factor; \
    for (size_t i_ = 0; i_ < curr_n; ++i_) curr_buffer[i_].weight *= factor; \
  }
#define MULT_COEFF(edge_idx) \
  { \
    double factor = 1.0 / (2.0 * bw[edge_idx] * b_g[edge_idx][edge_idx]); \
    for (size_t i_ = 0; i_ < n; ++i_) \
      weight_vector_data_buffer[i_].weight *= factor; \
  }
#define COMPACT_WEIGHTS(edge_idx) \
  { \
    compact_weight_vector_data( \
      weight_vector_data_buffer, &n, weight_vector_data_pos_buffer); \
    curr_triangle->c_111_a[edge_idx].data = \
      xmalloc(n * sizeof(*(curr_triangle->c_111_a[edge_idx].data))); \
    memcpy(curr_triangle->c_111_a[edge_idx].data, weight_vector_data_buffer, \
           n * sizeof(*weight_vector_data_buffer)); \
    curr_triangle->c_111_a[edge_idx].n = n; \
  }

    // a_0 = (D_g f(w0)/3 -
    //        B_020(w0) * (b_0(g0)*c_120 + b_1(g0)*c_030 + b_2(g0)*c_021) -
    //        B_011(w0) * (                b_1(g0)*c_021 + b_2(g0)*c_012) -
    //        B_002(w0) * (b_0(g0)*c_102 + b_1(g0)*c_012 + b_2(g0)*c_003)) /
    //       (B_011(w0) * b_0(g0))
    size_t n = 0;
    COPY_D_DATA(0); // D_g f(w0)/3
    COPY_COEFF(c_120, 0, 0, 1.0); // - B_020(w0) * b_0(g0) * c_120
    COPY_COEFF(c_030, 0, 1, 1.0); // - B_020(w0) * b_1(g0) * c_030
    COPY_COEFF(c_021, 0, 2, 1.0); // - B_020(w0) * b_2(g0) * c_021
    COPY_COEFF(c_021, 0, 1, 2.0); // - B_011(w0) * b_1(g0) * c_021
    COPY_COEFF(c_012, 0, 2, 2.0); // - B_011(w0) * b_2(g0) * c_012
    COPY_COEFF(c_102, 0, 0, 1.0); // - B_002(w0) * b_0(g0) * c_102
    COPY_COEFF(c_012, 0, 1, 1.0); // - B_002(w0) * b_1(g0) * c_012
    COPY_COEFF(c_003, 0, 2, 1.0); // - B_002(w0) * b_2(g0) * c_003
    MULT_COEFF(0); // (...) / (B_011(w0) * b_0(g0))
    COMPACT_WEIGHTS(0); // a_0 = ...

    // a_1 = (D_g f(w1)/3 -
    //        B_200(w1) * (b_0(g1)*c_300 + b_1(g1)*c_210 + b_2(g1)*c_201) -
    //        B_101(w1) * (b_0(g1)*c_201 +                 b_2(g1)*c_102) -
    //        B_002(w1) * (b_0(g1)*c_102 + b_1(g1)*c_012 + b_2(g1)*c_003)) /
    //       (B_101(w1) * b_1(g1))
    n = 0;
    COPY_D_DATA(1); // D_g f(w1)/3
    COPY_COEFF(c_300, 1, 0, 1.0); // - B_200(w1) * b_0(g1) * c_300
    COPY_COEFF(c_210, 1, 1, 1.0); // - B_200(w1) * b_1(g1) * c_210
    COPY_COEFF(c_201, 1, 2, 1.0); // - B_200(w1) * b_2(g1) * c_201
    COPY_COEFF(c_201, 1, 0, 2.0); // - B_101(w1) * b_0(g1) * c_201
    COPY_COEFF(c_102, 1, 2, 2.0); // - B_101(w1) * b_2(g1) * c_102
    COPY_COEFF(c_102, 1, 0, 1.0); // - B_002(w1) * b_0(g1) * c_102
    COPY_COEFF(c_012, 1, 1, 1.0); // - B_002(w1) * b_1(g1) * c_012
    COPY_COEFF(c_003, 1, 2, 1.0); // - B_002(w1) * b_2(g1) * c_003
    MULT_COEFF(1); // (...) / (B_101(w1) * b_1(g1))
    COMPACT_WEIGHTS(1); // a_1 = ...

    // a_2 = (D_g f(w2)/3 -
    //        B_200(w2) * (b_0(g2)*c_300 + b_1(g2)*c_210 + b_2(g2)*c_201) -
    //        B_110(w2) * (b_0(g2)*c_210 + b_1(g2)*c_120                ) -
    //        B_020(w2) * (b_0(g2)*c_120 + b_1(g2)*c_030 + b_2(g2)*c_021)) /
    //       (B_110(w2) * b_2(g2))
    n = 0;
    COPY_D_DATA(2); // D_g f(w2)/3
    COPY_COEFF(c_300, 2, 0, 1.0); // - B_200(w2) * b_0(g2) * c_300
    COPY_COEFF(c_210, 2, 1, 1.0); // - B_200(w2) * b_1(g2) * c_210
    COPY_COEFF(c_201, 2, 2, 1.0); // - B_200(w2) * b_2(g2) * c_201
    COPY_COEFF(c_210, 2, 0, 2.0); // - B_110(w2) * b_0(g2) * c_210
    COPY_COEFF(c_120, 2, 1, 2.0); // - B_110(w2) * b_1(g2) * c_120
    COPY_COEFF(c_120, 2, 0, 1.0); // - B_020(w2) * b_0(g2) * c_120
    COPY_COEFF(c_030, 2, 1, 1.0); // - B_020(w2) * b_1(g2) * c_030
    COPY_COEFF(c_021, 2, 2, 1.0); // - B_020(w2) * b_2(g2) * c_021
    MULT_COEFF(2); // (...) / (B_110(w2) * b_2(g2))
    COMPACT_WEIGHTS(2); // a_2 = ...

#undef COMPACT_WEIGHTS
#undef MULT_COEFF
#undef COPY_COEFF
#undef COPY_D_DATA
  }
  free(weight_vector_data_buffer);
  free(weight_vector_data_pos_buffer);
}

static size_t get_max_num_weights(
  struct tgt_point_search_data * tgt_point_data) {

  YAC_ASSERT(
    (tgt_point_data->type_flag == ON_VERTEX) ||
    (tgt_point_data->type_flag == ON_EDGE) ||
    (tgt_point_data->type_flag == ON_TRIANGLE),
    "ERROR(get_max_num_weights): invalid type_flag")
  switch(tgt_point_data->type_flag) {

    case (ON_VERTEX):
      return 1;
    case (ON_EDGE):
      return 2 + tgt_point_data->interp_data.edge->c[0].n +
             tgt_point_data->interp_data.edge->c[1].n;
    default:
    case (ON_TRIANGLE): {
      struct triangle_interp_data * triangle_data =
        tgt_point_data->interp_data.triangle;
      return 3 +
             triangle_data->edges[0]->c[0].n + triangle_data->edges[0]->c[1].n +
             triangle_data->edges[1]->c[0].n + triangle_data->edges[1]->c[1].n +
             triangle_data->edges[2]->c[0].n + triangle_data->edges[2]->c[1].n +
             triangle_data->c_111_a[0].n +
             triangle_data->c_111_a[1].n +
             triangle_data->c_111_a[2].n;
    }
  };
}

static void compute_hcsbb_weights_vertex(
  struct tgt_point_search_data * tgt_point_data,
  struct weight_vector_data * weights, size_t * n) {

  weights[0].weight = 1.0;
  weights[0].point_id = tgt_point_data->src_points[0];
  *n = 1;
}

static void compute_hcsbb_weights_edge(
  struct tgt_point_search_data * tgt_point_data,
  struct weight_vector_data * weights, size_t * n) {

  // spherical barycentric coordinates of the target point
  double * sb_coords = tgt_point_data->sb_coords;

  // compute the spherical Bernstein base polynomials for the given vertices
  // (since the target point is on the edge, we only need a subset of the
  //  polynomials, because the rest are zero anyway)
  // B_ijk(p) = d!/(i!*j!*k!) * b_0(p)^i * b_1(p)^j * b_2(p)^k
  double B_300 = sb_coords[0] * sb_coords[0] * sb_coords[0];
  double B_030 = sb_coords[1] * sb_coords[1] * sb_coords[1];
  double B_210 = 3.0 * sb_coords[0] * sb_coords[0] * sb_coords[1];
  double B_120 = 3.0 * sb_coords[0] * sb_coords[1] * sb_coords[1];

  // c_300 = f(v1)
  // c_030 = f(v2)
  // c_210 and c_120
  struct weight_vector * c = &(tgt_point_data->interp_data.edge->c[0]);

  // f(p) = SUM(c_ijk * B_ijk) // for 300, 030, 210, and 120

  size_t * src_points = tgt_point_data->src_points;

  // c_300 * B_300
  // weights[0].weight = B_300 + B_210;
  weights[0].weight = B_300;
  weights[0].point_id = src_points[0];
  // c_030 * B_030
  // weights[1].weight = B_030 + B_120;
  weights[1].weight = B_030;
  weights[1].point_id = src_points[1];
  // c_210 * B_210
  struct weight_vector_data * curr_weights = weights + 2;
  memcpy(curr_weights, c[0].data, c[0].n * sizeof(*curr_weights));
  for (size_t i = 0; i < c[0].n; ++i) curr_weights[i].weight *= B_210;
  // c_120 * B_120
  curr_weights += c[0].n;
  memcpy(curr_weights, c[1].data, c[1].n * sizeof(*curr_weights));
  for (size_t i = 0; i < c[1].n; ++i) curr_weights[i].weight *= B_120;

  *n = 2 + c[0].n + c[1].n;
}

static void evaluate_blending_function(
  double * sb_coords, double * A) {

  double b[3] = {sb_coords[1]*sb_coords[1] * sb_coords[2]*sb_coords[2],
                 sb_coords[0]*sb_coords[0] * sb_coords[2]*sb_coords[2],
                 sb_coords[0]*sb_coords[0] * sb_coords[1]*sb_coords[1]};
  double temp = 1.0 / (b[0] + b[1] + b[2]);

  for (size_t i = 0; i < 3; ++i) {

    if (fabs(b[i]) < 1e-9) A[i] = 0.0;
    else                   A[i] = b[i] * temp;
  }
}

static void compute_hcsbb_weights_triangle(
  struct tgt_point_search_data * tgt_point_data,
  struct weight_vector_data * weights, size_t * n) {

  // spherical barycentric coordinates of the target point
  double * sb_coords = tgt_point_data->sb_coords;

  struct edge_interp_data * curr_edges[3] = {
    tgt_point_data->interp_data.triangle->edges[0],
    tgt_point_data->interp_data.triangle->edges[1],
    tgt_point_data->interp_data.triangle->edges[2]};
    struct vertex_interp_data * vertex_d_data[3] = {
      curr_edges[2]->vertex_d_data[0],
      curr_edges[2]->vertex_d_data[1],
      curr_edges[1]->vertex_d_data[1]};

  // compute the spherical Bernstein base polynomials for the given vertices
  // B_ijk(p) = d!/(i!*j!*k!) * b_0(p)^i * b_1(p)^j * b_2(p)^k
  double B_vertex[3] =
    {sb_coords[0] * sb_coords[0] * sb_coords[0],  // B_300
     sb_coords[1] * sb_coords[1] * sb_coords[1],  // B_030
     sb_coords[2] * sb_coords[2] * sb_coords[2]}; // B_003
  double B_edge[3][2] =
    {{3.0 * sb_coords[1] * sb_coords[1] * sb_coords[2],   // B_021
      3.0 * sb_coords[1] * sb_coords[2] * sb_coords[2]},  // B_012
     {3.0 * sb_coords[0] * sb_coords[0] * sb_coords[2],   // B_201
      3.0 * sb_coords[0] * sb_coords[2] * sb_coords[2]},  // B_102
     {3.0 * sb_coords[0] * sb_coords[0] * sb_coords[1],   // B_210
      3.0 * sb_coords[0] * sb_coords[1] * sb_coords[1]}}; // B_120
  double B_111 = 6.0 * sb_coords[0] * sb_coords[1] * sb_coords[2];

  // f(p) = SUM(c_ijk * B_ijk) // for all (i+j+k) = 3

  // the three vertices coefficients (300, 030, and 003)
  for (size_t i = 0; i < 3; ++i) {
    weights[i].weight = B_vertex[i];
    weights[i].point_id = vertex_d_data[i]->src_point;
  }

  size_t n_ = 3;

  // the six edges coefficients
  struct weight_vector * c_ijk[3][2] = {
    {&(curr_edges[0]->c[0]),  // c_021
     &(curr_edges[0]->c[1])}, // c_012
    {&(curr_edges[1]->c[0]),  // c_201
     &(curr_edges[1]->c[1])}, // c_102
    {&(curr_edges[2]->c[0]),  // c_210
     &(curr_edges[2]->c[1])}  // c_120
  };
  for (size_t edge_idx = 0; edge_idx < 3; ++edge_idx) {
    for (size_t i = 0; i < 2; ++i) {
      struct weight_vector_data * curr_weights = weights + n_;
      size_t curr_n = c_ijk[edge_idx][i]->n;
      n_ += curr_n;
      memcpy(curr_weights, c_ijk[edge_idx][i]->data,
             c_ijk[edge_idx][i]->n * sizeof(*curr_weights));
      double B = B_edge[edge_idx][i];
      for (size_t j = 0; j < curr_n; ++j)
        curr_weights[j].weight *= B;
    }
  }

  // the middle vertex
  // c_111 = SUM(c_111_a[i] * A[i]) // for i = 0..2
  // where a is the blending function
  double A[3];
  evaluate_blending_function(sb_coords, &(A[0]));

  struct weight_vector * c_111_a =
    tgt_point_data->interp_data.triangle->c_111_a;
  for (size_t i = 0; i < 3; ++i) {
    struct weight_vector_data * curr_weights = weights + n_;
    size_t curr_n = c_111_a[i].n;
    n_ += curr_n;
    memcpy(curr_weights, c_111_a[i].data, curr_n * sizeof(*curr_weights));
    double B_111_x_A_i = B_111 * A[i];
    for (size_t j = 0; j < curr_n; ++j)
      curr_weights[j].weight *= B_111_x_A_i;
  }

  *n = n_;
}

static void compute_hcsbb_weights(
  struct tgt_point_search_data * tgt_point_data,
  struct weight_vector_data * weights, size_t * n,
  struct weight_vector_data_pos * compact_buffer) {

  YAC_ASSERT(
    (tgt_point_data->type_flag == ON_VERTEX) ||
    (tgt_point_data->type_flag == ON_EDGE) ||
    (tgt_point_data->type_flag == ON_TRIANGLE),
    "ERROR(compute_hcsbb_weights): invalid type_flag")

  switch(tgt_point_data->type_flag) {

    case (ON_VERTEX):
      compute_hcsbb_weights_vertex(tgt_point_data, weights, n);
      break;
    case (ON_EDGE):
      compute_hcsbb_weights_edge(tgt_point_data, weights, n);
      break;
    default:
    case (ON_TRIANGLE):
      compute_hcsbb_weights_triangle(tgt_point_data, weights, n);
      break;
  }

  compact_weight_vector_data(weights, n, compact_buffer);
}

static void compute_weights(
  struct tgt_point_search_data * tgt_point_data, size_t num_tgt_points,
  struct edge_interp_data * edge_data, size_t num_edges,
  struct triangle_interp_data * triangle_data, size_t num_triangles,
  struct weight_vector_data ** weights, size_t ** num_weights_per_tgt,
  size_t * total_num_weights) {

  //--------------------------------------------------------------------------
  // prepare patch data for the computation of the weights for each individual
  // target point
  //--------------------------------------------------------------------------

  // compute edge coefficients
  compute_edge_coefficients(edge_data, num_edges);

  // compute triangle coefficients
  compute_triangle_coefficients(triangle_data, num_triangles);

  struct weight_vector_data * weight_vector_data = NULL;
  size_t weight_vector_data_array_size = 0;
  size_t total_num_weights_ = 0;
  size_t * num_weights_per_tgt_ =
    xmalloc(num_tgt_points * sizeof(*num_weights_per_tgt_));
  struct weight_vector_data_pos * compact_buffer =
    xmalloc((3 * 1 + // c_300, c_030, c_003
             3 * (1 + HCSBB_NUM_GAUSS_POINTS * HCSBB_GAUSS_NNN) + // c_012, c_021, c_102, ...
             3 * (HCSBB_NUM_GAUSS_POINTS * HCSBB_GAUSS_NNN + 2 +
                  6 * (1 + HCSBB_NUM_GAUSS_POINTS * HCSBB_GAUSS_NNN))) * // c_111_a[3]
             sizeof(*compact_buffer));

  //--------------------------------------
  // compute weights for the target points
  //--------------------------------------

  for (size_t i = 0; i < (num_tgt_points + 127) / 128; ++i) {

    size_t curr_num_tgt_points = MIN(128, num_tgt_points - i * 128);
    struct tgt_point_search_data * curr_tgt_point_data =
      tgt_point_data + i * 128;
    size_t * curr_num_weights_per_tgt = num_weights_per_tgt_ + i * 128;

    size_t curr_max_num_weights = 0;
    for (size_t j = 0; j < curr_num_tgt_points; ++j)
      curr_max_num_weights += get_max_num_weights(curr_tgt_point_data + j);

    ENSURE_ARRAY_SIZE(weight_vector_data, weight_vector_data_array_size,
                      total_num_weights_ + curr_max_num_weights);

    for (size_t j = 0; j < curr_num_tgt_points; ++j) {
      compute_hcsbb_weights(
        curr_tgt_point_data + j, weight_vector_data + total_num_weights_,
        curr_num_weights_per_tgt + j, compact_buffer);
      total_num_weights_ += curr_num_weights_per_tgt[j];
    }
  }

  free(compact_buffer);

  *weights =
    xrealloc(
      weight_vector_data, total_num_weights_ * sizeof(*weight_vector_data));
  *num_weights_per_tgt = num_weights_per_tgt_;
  *total_num_weights = total_num_weights_;
}

static inline int compare_tgt_point_search_data(
  const void * a, const void * b) {

  struct tgt_point_search_data const * a_ = a;
  struct tgt_point_search_data const * b_ = b;

  int ret;
  if ((ret = a_->is_masked - b_->is_masked)) return ret;
  if (a_->is_masked) return 0;
  if ((ret = (int)(a_->type_flag) - (int)(b_->type_flag))) return ret;
  for (int i = 0; i <= (int)a_->type_flag; ++i)
    if ((ret = (a_->src_points[i] > b_->src_points[i]) -
               (a_->src_points[i] < b_->src_points[i]))) return ret;
  return (a_->tgt_point > b_->tgt_point) - (a_->tgt_point < b_->tgt_point);
}

static void generate_bnd_triangle(
  size_t * points, const_coordinate_pointer coords,
  double bnd_triangle[][3]) {

  double point_coords[HCSBB_D_NNN][3];

  for (size_t i = 0; i < HCSBB_D_NNN; ++i)
    memcpy(point_coords[i], coords[points[i]], 3 * sizeof(double));

  yac_compute_bnd_triangle(
    &(point_coords[0][0]), HCSBB_D_NNN, bnd_triangle, 16);
}

static int compare_bnd_triangles(void const * a, void const * b) {

  double const * tri_a = a,  * tri_b = b;
  int ret = 0;

  for (size_t i = 0; (!ret) && (i < 9); ++i) {
    double a_ = tri_a[i], b_ = tri_b[i];
    if (fabs(a_ - b_) > 1e-9) ret = (a_ > b_) - (a_ < b_);
  }
  return ret;
}

static void generate_gauss_legendre_points(
  double gauss_vertices[HCSBB_NUM_GAUSS_POINTS][3],
  double gauss_sb_coords[HCSBB_NUM_GAUSS_POINTS][3],
  double bnd_triangle[3][3]) {

#if HCSBB_GAUSS_ORDER == 3
  static double const abscissa[HCSBB_NUM_GAUSS_POINTS][3] = {
    {1.0-0.333333333333333-0.333333333333333,0.333333333333333,0.333333333333333},
    {1.0-0.200000000000000-0.200000000000000,0.200000000000000,0.200000000000000},
    {1.0-0.600000000000000-0.200000000000000,0.600000000000000,0.200000000000000},
    {1.0-0.200000000000000-0.600000000000000,0.200000000000000,0.600000000000000}
  };
#elif HCSBB_GAUSS_ORDER == 4
  static double const abscissa[HCSBB_NUM_GAUSS_POINTS][3] = {
    {1.0-0.091576213509771-0.091576213509771,0.091576213509771,0.091576213509771},
    {1.0-0.816847572980459-0.091576213509771,0.816847572980459,0.091576213509771},
    {1.0-0.091576213509771-0.816847572980459,0.091576213509771,0.816847572980459},
    {1.0-0.445948490915965-0.445948490915965,0.445948490915965,0.445948490915965},
    {1.0-0.108103018168070-0.445948490915965,0.108103018168070,0.445948490915965},
    {1.0-0.445948490915965-0.108103018168070,0.445948490915965,0.108103018168070}
  };
#elif HCSBB_GAUSS_ORDER == 5
  static double const abscissa[HCSBB_NUM_GAUSS_POINTS][3] = {
    {1.0-0.333333333333333-0.333333333333333,0.333333333333333,0.333333333333333},
    {1.0-0.101286507323456-0.101286507323456,0.101286507323456,0.101286507323456},
    {1.0-0.797426985353087-0.101286507323456,0.797426985353087,0.101286507323456},
    {1.0-0.101286507323456-0.797426985353087,0.101286507323456,0.797426985353087},
    {1.0-0.470142064105115-0.470142064105115,0.470142064105115,0.470142064105115},
    {1.0-0.059715871789770-0.470142064105115,0.059715871789770,0.470142064105115},
    {1.0-0.470142064105115-0.059715871789770,0.470142064105115,0.059715871789770}
  };
#elif HCSBB_GAUSS_ORDER == 6
  static double const abscissa[HCSBB_NUM_GAUSS_POINTS][3] = {
    {1.0-0.063089014491502-0.063089014491502,0.063089014491502,0.063089014491502},
    {1.0-0.873821971016996-0.063089014491502,0.873821971016996,0.063089014491502},
    {1.0-0.063089014491502-0.873821971016996,0.063089014491502,0.873821971016996},
    {1.0-0.249286745170910-0.249286745170910,0.249286745170910,0.249286745170910},
    {1.0-0.501426509658179-0.249286745170910,0.501426509658179,0.249286745170910},
    {1.0-0.249286745170910-0.501426509658179,0.249286745170910,0.501426509658179},
    {1.0-0.310352451033785-0.053145049844816,0.310352451033785,0.053145049844816},
    {1.0-0.053145049844816-0.310352451033785,0.053145049844816,0.310352451033785},
    {1.0-0.636502499121399-0.053145049844816,0.636502499121399,0.053145049844816},
    {1.0-0.053145049844816-0.636502499121399,0.053145049844816,0.636502499121399},
    {1.0-0.636502499121399-0.310352451033785,0.636502499121399,0.310352451033785},
    {1.0-0.310352451033785-0.636502499121399,0.310352451033785,0.636502499121399}
  };
#elif HCSBB_GAUSS_ORDER == 7
  static double const abscissa[HCSBB_NUM_GAUSS_POINTS][3] = {
    {1.0-0.333333333333333-0.333333333333333,0.333333333333333,0.333333333333333},
    {1.0-0.260345966079038-0.260345966079038,0.260345966079038,0.260345966079038},
    {1.0-0.479308067841923-0.260345966079038,0.479308067841923,0.260345966079038},
    {1.0-0.260345966079038-0.479308067841923,0.260345966079038,0.479308067841923},
    {1.0-0.065130102902216-0.065130102902216,0.065130102902216,0.065130102902216},
    {1.0-0.869739794195568-0.065130102902216,0.869739794195568,0.065130102902216},
    {1.0-0.065130102902216-0.869739794195568,0.065130102902216,0.869739794195568},
    {1.0-0.312865496004874-0.048690315425316,0.312865496004874,0.048690315425316},
    {1.0-0.048690315425316-0.312865496004874,0.048690315425316,0.312865496004874},
    {1.0-0.638444188569809-0.048690315425316,0.638444188569809,0.048690315425316},
    {1.0-0.048690315425316-0.638444188569809,0.048690315425316,0.638444188569809},
    {1.0-0.638444188569809-0.312865496004874,0.638444188569809,0.312865496004874},
    {1.0-0.312865496004874-0.638444188569809,0.312865496004874,0.638444188569809}
  };
#elif HCSBB_GAUSS_ORDER == 8
  static double const abscissa[HCSBB_NUM_GAUSS_POINTS][3] = {
    {1.0-0.333333333333333-0.333333333333333,0.333333333333333,0.333333333333333},
    {1.0-0.081414823414554-0.459292588292723,0.081414823414554,0.459292588292723},
    {1.0-0.459292588292723-0.081414823414554,0.459292588292723,0.081414823414554},
    {1.0-0.459292588292723-0.459292588292723,0.459292588292723,0.459292588292723},
    {1.0-0.658861384496480-0.170569307751760,0.658861384496480,0.170569307751760},
    {1.0-0.170569307751760-0.658861384496480,0.170569307751760,0.658861384496480},
    {1.0-0.170569307751760-0.170569307751760,0.170569307751760,0.170569307751760},
    {1.0-0.898905543365938-0.050547228317031,0.898905543365938,0.050547228317031},
    {1.0-0.050547228317031-0.898905543365938,0.050547228317031,0.898905543365938},
    {1.0-0.050547228317031-0.050547228317031,0.050547228317031,0.050547228317031},
    {1.0-0.008394777409958-0.263112829634638,0.008394777409958,0.263112829634638},
    {1.0-0.263112829634638-0.008394777409958,0.263112829634638,0.008394777409958},
    {1.0-0.008394777409958-0.728492392955404,0.008394777409958,0.728492392955404},
    {1.0-0.728492392955404-0.008394777409958,0.728492392955404,0.008394777409958},
    {1.0-0.263112829634638-0.728492392955404,0.263112829634638,0.728492392955404},
    {1.0-0.728492392955404-0.263112829634638,0.728492392955404,0.263112829634638}
  };
#else
  static double const abscissa[1][3];
#error "interpolation_method_hcsbb: the specified gauss order is currently not supported."
#endif

  for (size_t i = 0; i < HCSBB_NUM_GAUSS_POINTS; ++i) {

    gauss_vertices[i][0] = bnd_triangle[0][0] * abscissa[i][0] +
                           bnd_triangle[1][0] * abscissa[i][1] +
                           bnd_triangle[2][0] * abscissa[i][2];
    gauss_vertices[i][1] = bnd_triangle[0][1] * abscissa[i][0] +
                           bnd_triangle[1][1] * abscissa[i][1] +
                           bnd_triangle[2][1] * abscissa[i][2];
    gauss_vertices[i][2] = bnd_triangle[0][2] * abscissa[i][0] +
                           bnd_triangle[1][2] * abscissa[i][1] +
                           bnd_triangle[2][2] * abscissa[i][2];
    double scale = 1.0 / sqrt(gauss_vertices[i][0] * gauss_vertices[i][0] +
                              gauss_vertices[i][1] * gauss_vertices[i][1] +
                              gauss_vertices[i][2] * gauss_vertices[i][2]);
    gauss_vertices[i][0] *= scale;
    gauss_vertices[i][1] *= scale;
    gauss_vertices[i][2] *= scale;
    gauss_sb_coords[i][0] = abscissa[i][0] * scale;
    gauss_sb_coords[i][1] = abscissa[i][1] * scale;
    gauss_sb_coords[i][2] = abscissa[i][2] * scale;
  }
}

// computes the spherical Bernstein basis polynomials
// for degree = 3 the polynomials B_ijk are ordered in the following order:
// B_003, B_012, B_021, B_030, B_102, B_111, B_120, B_201, B_210, B_300
static void compute_sbb_polynomials_gauss_3d(
  double sb_coords[HCSBB_NUM_GAUSS_POINTS][3],
  double sbb_polynomials[HCSBB_NUM_SBB_COEFFS][HCSBB_NUM_GAUSS_POINTS]) {

#define POW_0(x) (1.0)
#define POW_1(x) ((x))
#define POW_2(x) ((x)*(x))
#define POW_3(x) ((x)*(x)*(x))
#define FAC_0 (1.0)
#define FAC_1 (1.0)
#define FAC_2 (2.0)
#define FAC_3 (6.0)

  for (size_t vertex_idx = 0; vertex_idx < HCSBB_NUM_GAUSS_POINTS; ++vertex_idx)
    sbb_polynomials[0][vertex_idx] =
      (FAC_3 / (FAC_0 * FAC_0 * FAC_3)) * (POW_0(sb_coords[vertex_idx][0]) *
                                           POW_0(sb_coords[vertex_idx][1]) *
                                           POW_3(sb_coords[vertex_idx][2]));
  for (size_t vertex_idx = 0; vertex_idx < HCSBB_NUM_GAUSS_POINTS; ++vertex_idx)
    sbb_polynomials[1][vertex_idx] =
      (FAC_3 / (FAC_0 * FAC_1 * FAC_2)) * (POW_0(sb_coords[vertex_idx][0]) *
                                           POW_1(sb_coords[vertex_idx][1]) *
                                           POW_2(sb_coords[vertex_idx][2]));
  for (size_t vertex_idx = 0; vertex_idx < HCSBB_NUM_GAUSS_POINTS; ++vertex_idx)
    sbb_polynomials[2][vertex_idx] =
      (FAC_3 / (FAC_0 * FAC_2 * FAC_1)) * (POW_0(sb_coords[vertex_idx][0]) *
                                           POW_2(sb_coords[vertex_idx][1]) *
                                           POW_1(sb_coords[vertex_idx][2]));
  for (size_t vertex_idx = 0; vertex_idx < HCSBB_NUM_GAUSS_POINTS; ++vertex_idx)
    sbb_polynomials[3][vertex_idx] =
      (FAC_3 / (FAC_0 * FAC_3 * FAC_0)) * (POW_0(sb_coords[vertex_idx][0]) *
                                           POW_3(sb_coords[vertex_idx][1]) *
                                           POW_0(sb_coords[vertex_idx][2]));
  for (size_t vertex_idx = 0; vertex_idx < HCSBB_NUM_GAUSS_POINTS; ++vertex_idx)
    sbb_polynomials[4][vertex_idx] =
      (FAC_3 / (FAC_1 * FAC_0 * FAC_2)) * (POW_1(sb_coords[vertex_idx][0]) *
                                           POW_0(sb_coords[vertex_idx][1]) *
                                           POW_2(sb_coords[vertex_idx][2]));
  for (size_t vertex_idx = 0; vertex_idx < HCSBB_NUM_GAUSS_POINTS; ++vertex_idx)
    sbb_polynomials[5][vertex_idx] =
      (FAC_3 / (FAC_1 * FAC_1 * FAC_1)) * (POW_1(sb_coords[vertex_idx][0]) *
                                           POW_1(sb_coords[vertex_idx][1]) *
                                           POW_1(sb_coords[vertex_idx][2]));
  for (size_t vertex_idx = 0; vertex_idx < HCSBB_NUM_GAUSS_POINTS; ++vertex_idx)
    sbb_polynomials[6][vertex_idx] =
      (FAC_3 / (FAC_1 * FAC_2 * FAC_0)) * (POW_1(sb_coords[vertex_idx][0]) *
                                           POW_2(sb_coords[vertex_idx][1]) *
                                           POW_0(sb_coords[vertex_idx][2]));
  for (size_t vertex_idx = 0; vertex_idx < HCSBB_NUM_GAUSS_POINTS; ++vertex_idx)
    sbb_polynomials[7][vertex_idx] =
      (FAC_3 / (FAC_2 * FAC_0 * FAC_1)) * (POW_2(sb_coords[vertex_idx][0]) *
                                           POW_0(sb_coords[vertex_idx][1]) *
                                           POW_1(sb_coords[vertex_idx][2]));
  for (size_t vertex_idx = 0; vertex_idx < HCSBB_NUM_GAUSS_POINTS; ++vertex_idx)
    sbb_polynomials[8][vertex_idx] =
      (FAC_3 / (FAC_2 * FAC_1 * FAC_0)) * (POW_2(sb_coords[vertex_idx][0]) *
                                           POW_1(sb_coords[vertex_idx][1]) *
                                           POW_0(sb_coords[vertex_idx][2]));
  for (size_t vertex_idx = 0; vertex_idx < HCSBB_NUM_GAUSS_POINTS; ++vertex_idx)
    sbb_polynomials[9][vertex_idx] =
      (FAC_3 / (FAC_3 * FAC_0 * FAC_0)) * (POW_3(sb_coords[vertex_idx][0]) *
                                           POW_0(sb_coords[vertex_idx][1]) *
                                           POW_0(sb_coords[vertex_idx][2]));

#undef POW_0
#undef POW_1
#undef POW_2
#undef POW_3
#undef FAC_0
#undef FAC_1
#undef FAC_2
#undef FAC_3
}

static void inverse(double A[HCSBB_NUM_SBB_COEFFS][HCSBB_NUM_SBB_COEFFS]) {

// LAPACKE_dsytrf_work and LAPACKE_dsytri might not be available (see yac_lapack_interface.h).
#ifdef YAC_LAPACK_NO_DSYTR
  lapack_int ipiv[HCSBB_NUM_SBB_COEFFS + 1] = {0};
  double work[HCSBB_NUM_SBB_COEFFS*HCSBB_NUM_SBB_COEFFS] = {0};

  YAC_ASSERT(
    !LAPACKE_dgetrf(
      LAPACK_COL_MAJOR, (lapack_int) HCSBB_NUM_SBB_COEFFS,
      (lapack_int) HCSBB_NUM_SBB_COEFFS, &A[0][0],
      (lapack_int) HCSBB_NUM_SBB_COEFFS, ipiv), "internal ERROR: dgetrf")

  YAC_ASSERT(
    !LAPACKE_dgetri_work(
      LAPACK_COL_MAJOR, (lapack_int) HCSBB_NUM_SBB_COEFFS, &A[0][0],
      (lapack_int) HCSBB_NUM_SBB_COEFFS, ipiv, work,
      (lapack_int) (HCSBB_NUM_SBB_COEFFS * HCSBB_NUM_SBB_COEFFS)),
      "internal ERROR: dgetri")
#else
  lapack_int ipiv[HCSBB_NUM_SBB_COEFFS] = {0};
  double work[HCSBB_NUM_SBB_COEFFS] = {0};

  YAC_ASSERT(
    !LAPACKE_dsytrf_work(
      LAPACK_COL_MAJOR, 'L', (lapack_int) HCSBB_NUM_SBB_COEFFS , &A[0][0],
      (lapack_int) HCSBB_NUM_SBB_COEFFS, ipiv, work,
      (lapack_int) HCSBB_NUM_SBB_COEFFS), "internal ERROR: dsytrf")

  YAC_ASSERT(
    !LAPACKE_dsytri_work(
      LAPACK_COL_MAJOR, 'L', (lapack_int) HCSBB_NUM_SBB_COEFFS , &A[0][0],
      (lapack_int) HCSBB_NUM_SBB_COEFFS, ipiv, work),
      "internal ERROR: dsytri")

  for (size_t i = 0; i < HCSBB_NUM_SBB_COEFFS; ++i)
    for (size_t j = i + 1; j < HCSBB_NUM_SBB_COEFFS; ++j)
      A[j][i] = A[i][j];
#endif
}

static void init_gauss_nnn_patch(
  struct gauss_nnn_patch * gauss_nnn_patch, double (*bnd_triangle)[3],
  double (*gauss_vertices)[3],
  double gauss_sbb_polynomials[HCSBB_NUM_SBB_COEFFS][HCSBB_NUM_GAUSS_POINTS],
  double C[HCSBB_NUM_SBB_COEFFS][HCSBB_NUM_SBB_COEFFS]) {

  // compute bounding triangle
  memcpy(&(gauss_nnn_patch->bnd_triangle[0][0]), bnd_triangle,
         3 * sizeof(*bnd_triangle));

  double gauss_sb_coords[HCSBB_NUM_GAUSS_POINTS][3];

  // compute gauss points
  generate_gauss_legendre_points(gauss_vertices, gauss_sb_coords, bnd_triangle);

  double (*w_g)[HCSBB_NUM_GAUSS_POINTS] =
    xmalloc(HCSBB_NUM_SBB_COEFFS * sizeof(w_g[0]));
  gauss_nnn_patch->weights = (double*)w_g;

  // compute w_g
  {
    // compute the spherical Bernstein basis polynomials for the
    // gauss points (P_g)
    compute_sbb_polynomials_gauss_3d(gauss_sb_coords, gauss_sbb_polynomials);

    // P_g^T * P_g
    // (actually this is a triangular matrix -> computing halve would be enough)
    for (size_t i = 0; i < HCSBB_NUM_SBB_COEFFS; ++i) {
      for (size_t j = 0; j < HCSBB_NUM_SBB_COEFFS; ++j) {
        double accu = 0;
        for (size_t k = 0; k < HCSBB_NUM_GAUSS_POINTS; ++k) {
          accu += gauss_sbb_polynomials[i][k] * gauss_sbb_polynomials[j][k];
        }
        C[i][j] = accu;
      }
    }

    // (P_p^T * P_p)^-1
    inverse(C);

    // (P_p^T * P_p)^-1 * P_p^T
    for (size_t i = 0; i < HCSBB_NUM_SBB_COEFFS; ++i) {
      for (size_t j = 0; j < HCSBB_NUM_GAUSS_POINTS; ++j) {
        double accu = 0;
        for (size_t k = 0; k < HCSBB_NUM_SBB_COEFFS; ++k) {
          accu += C[i][k] * gauss_sbb_polynomials[k][j];
        }
        w_g[i][j] = accu;
      }
    }
  }
}

static void free_gauss_nnn_patches(
  struct gauss_nnn_patch * gauss_nnn_patches, size_t num_gauss_nnn_patches) {

  for (size_t i = 0; i < num_gauss_nnn_patches; ++i) {

    free(gauss_nnn_patches[i].src_points);
    free(gauss_nnn_patches[i].weights);
  }

  free(gauss_nnn_patches);
}

// This routine generates a triangle around each coordinate. The generated
// triangle contains at around HCSBB_D_NNN source points, this ensures that the
// size of the triangle scales with the resolution of the source grid in the
// area around the coordinate.
static void generate_d_data_triangle(
  size_t count, coordinate_pointer coords, struct interp_grid * interp_grid,
  size_t * num_unique_triangles_, double (**unique_bnd_triangles_)[3][3],
  size_t ** coord_to_triangle) {

  // for all points at which we need to estimate the derivative, do a nearest
  // neighbour search
  size_t * nnn_search_results =
    xmalloc((HCSBB_D_NNN + 1) * count * sizeof(*nnn_search_results));
  yac_interp_grid_do_nnn_search_src(
    interp_grid, coords, count, HCSBB_D_NNN, nnn_search_results);

  // sort the results for each search coordinate by the global ids of the
  // result ids (we do not actually care about the distance to the search
  // point (we also add a reorder index)
  const_yac_int_pointer src_global_ids =
    yac_interp_grid_get_src_field_global_ids(interp_grid, 0);
  for (size_t i = 0, j = count - 1; i < count; ++i, --j) {
    // get global_ids for the current search results
    size_t * curr_nnn_search_results =
      nnn_search_results + j * HCSBB_D_NNN;
    yac_int nnn_search_result_global_ids[HCSBB_D_NNN];
    for (int k = 0; k < HCSBB_D_NNN; ++k)
      nnn_search_result_global_ids[k] =
        src_global_ids[curr_nnn_search_results[k]];
    // sort results of current point
    yac_quicksort_index_yac_int_size_t(
      nnn_search_result_global_ids, HCSBB_D_NNN, curr_nnn_search_results);
    // make space for reorder index
    memmove(nnn_search_results + j * (HCSBB_D_NNN + 1),
            curr_nnn_search_results,
            HCSBB_D_NNN * sizeof(*nnn_search_results));
    // set reorder index
    nnn_search_results[j * (HCSBB_D_NNN + 1) + HCSBB_D_NNN] = j;
  }

  // sort all results in order to be able to finde duplicated results
  qsort(nnn_search_results, count,
        (HCSBB_D_NNN + 1) * sizeof(*nnn_search_results),
        compare_HCSBB_D_NNN_size_t);

  size_t * coord_to_search_results =
    xmalloc(count * sizeof(*coord_to_search_results));
  size_t num_unique_search_results = 0;
  size_t dummy_search_results[HCSBB_D_NNN + 1];
  for (size_t i = 0; i < HCSBB_D_NNN + 1; ++i)
    dummy_search_results[i] = SIZE_MAX;
  for (size_t i = 0, * prev_search_result = dummy_search_results;
       i < count; ++i) {

    size_t * curr_search_result = nnn_search_results + i * (HCSBB_D_NNN + 1);
    size_t coord_idx = curr_search_result[HCSBB_D_NNN];
    if (compare_HCSBB_D_NNN_size_t(prev_search_result, curr_search_result)) {
      prev_search_result = curr_search_result;
      ++num_unique_search_results;
    }
    coord_to_search_results[coord_idx] = num_unique_search_results - 1;
  }

  // Remark: we have to call yac_interp_grid_get_src_field_coords after the call
  // to yac_interp_grid_do_nnn_search_src, because the source grid might have been
  // changed by the search
  const_coordinate_pointer src_point_coords =
    yac_interp_grid_get_src_field_coords(interp_grid, 0);

  // generate bnd_triangles for all unique search results
  struct bnd_triangle_reorder * bnd_triangles_reorder =
    xmalloc(num_unique_search_results * sizeof(*bnd_triangles_reorder));
  for (size_t i = 0, j = 0, * prev_search_result = dummy_search_results;
       i < count; ++i) {

    size_t * curr_search_result = nnn_search_results + i * (HCSBB_D_NNN + 1);
    if (compare_HCSBB_D_NNN_size_t(prev_search_result, curr_search_result)) {
      generate_bnd_triangle(
        curr_search_result, src_point_coords, bnd_triangles_reorder[j].coords);
      bnd_triangles_reorder[j].reorder_idx = j;
      prev_search_result = curr_search_result;
      ++j;
    }
  }
  free(nnn_search_results);

  // sort generated bounding triangles
  qsort(bnd_triangles_reorder, num_unique_search_results,
        sizeof(*bnd_triangles_reorder), compare_bnd_triangles);

  // determine number of unique bnd_triangles
  size_t * search_result_to_bnd_triangle =
    xmalloc(num_unique_search_results * sizeof(*search_result_to_bnd_triangle));
  double (*unique_bnd_triangles)[3][3] = (double(*)[3][3])bnd_triangles_reorder;
  size_t num_unique_triangles = 0;
  double dummy_bnd_triangle[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};
  double (*prev_bnd_triangle)[3] = dummy_bnd_triangle;
  for (size_t i = 0; i < num_unique_search_results; ++i) {
    double (*curr_bnd_triangle)[3] = bnd_triangles_reorder[i].coords;
    size_t reorder_idx = bnd_triangles_reorder[i].reorder_idx;
    if (compare_bnd_triangles(prev_bnd_triangle, curr_bnd_triangle)) {
      memmove(unique_bnd_triangles[num_unique_triangles], curr_bnd_triangle,
              sizeof(*unique_bnd_triangles));
      ++num_unique_triangles;
      prev_bnd_triangle = curr_bnd_triangle;
    }
    search_result_to_bnd_triangle[reorder_idx] = num_unique_triangles - 1;
  }

  // convert coord_to_search_results to coord_to_triangle
  for (size_t i = 0; i < count; ++i)
    coord_to_search_results[i] =
      search_result_to_bnd_triangle[coord_to_search_results[i]];
  free(search_result_to_bnd_triangle);

  *num_unique_triangles_ = num_unique_triangles;
  *unique_bnd_triangles_ =
    xrealloc(unique_bnd_triangles,
             num_unique_triangles * sizeof(*unique_bnd_triangles));
  *coord_to_triangle =
    xrealloc(coord_to_search_results, count * sizeof(*coord_to_search_results));
}

// To estimate the derivate at a given point we first build triangle around the
// point. The size of this triangle scales with the source grid resolution.
// For each triangle we generate a number of supporting points, which are
// interpolated using a simple nearest neighbour approach. Using these
// supporting points we can interpolate the source field based on
// Bernstein basis polynomials. Based on these polynomials, we can generate the
// weights to compute derivatives of points.
static void generate_derivative_data(
  size_t num_vertices, struct vertex_interp_data * vertex_data,
  size_t num_edges, struct edge_interp_data * edge_data,
  struct interp_grid * interp_grid, size_t * num_gauss_nnn_patches_,
  struct gauss_nnn_patch ** gauss_nnn_patches_) {

  // coordinates at which we will need to estimate a derivative
  coordinate_pointer coords =
    xmalloc((num_vertices + num_edges) * sizeof(*coords));

  // get coordinates of all vertices and all edge middle points
  for (size_t i = 0; i < num_vertices; ++i)
    memcpy(coords[i], vertex_data[i].coordinate_xyz, sizeof(*coords));
  for (size_t i = 0; i < num_edges; ++i)
    memcpy(
      coords[i + num_vertices], edge_data[i].middle_d_data.coordinate_xyz,
      sizeof(*coords));

  size_t num_gauss_nnn_patches;
  double (*bnd_triangles)[3][3];
  size_t * coord_to_gauss_nnn_patch;
  generate_d_data_triangle(
    num_vertices + num_edges, coords, interp_grid,
    &num_gauss_nnn_patches, &bnd_triangles, &coord_to_gauss_nnn_patch);
  free(coords);

  struct gauss_nnn_patch * gauss_nnn_patches =
    xmalloc(num_gauss_nnn_patches * sizeof(*gauss_nnn_patches));
  double (*gauss_points)[HCSBB_NUM_GAUSS_POINTS][3] =
    xmalloc(num_gauss_nnn_patches * sizeof(*gauss_points));
  double (*gauss_sbb_polynomials_buffer)[HCSBB_NUM_GAUSS_POINTS] =
    xmalloc(HCSBB_NUM_SBB_COEFFS * sizeof(*gauss_sbb_polynomials_buffer));
  double (*work_buffer)[HCSBB_NUM_SBB_COEFFS] =
    xmalloc(HCSBB_NUM_SBB_COEFFS * sizeof(*work_buffer));

  *num_gauss_nnn_patches_ = num_gauss_nnn_patches;
  *gauss_nnn_patches_ = gauss_nnn_patches;

  // initialise all gauss patches
  for (size_t i = 0; i < num_gauss_nnn_patches; ++i)
    init_gauss_nnn_patch(
      gauss_nnn_patches + i, bnd_triangles[i],
      gauss_points[i], gauss_sbb_polynomials_buffer, work_buffer);
  free(work_buffer);
  free(gauss_sbb_polynomials_buffer);
  free(bnd_triangles);

  for (size_t i = 0; i < num_vertices; ++i)
    vertex_data[i].gauss_nnn_patch =
      gauss_nnn_patches + coord_to_gauss_nnn_patch[i];
  for (size_t i = 0; i < num_edges; ++i)
    edge_data[i].middle_d_data.gauss_nnn_patch =
      gauss_nnn_patches + coord_to_gauss_nnn_patch[i + num_vertices];
  free(coord_to_gauss_nnn_patch);

  // do an NNN search for all gauss integration point
  size_t * nnn_search_results =
    xmalloc(HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS * num_gauss_nnn_patches *
            sizeof(*nnn_search_results));
  yac_interp_grid_do_nnn_search_src(
    interp_grid, (coordinate_pointer)&(gauss_points[0][0][0]),
    HCSBB_NUM_GAUSS_POINTS * num_gauss_nnn_patches,
    HCSBB_GAUSS_NNN, nnn_search_results);

  // generate the nnn weights for the interpolation of the gauss integration
  // points which are required for the derivative estimates
  double (*w_nnn_buffer)[HCSBB_NUM_GAUSS_POINTS] =
    xmalloc(HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS * sizeof(*w_nnn_buffer));
  size_t * src_point_buffer =
    xmalloc(
      HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS * sizeof(*src_point_buffer));
  yac_int * global_id_buffer =
    xmalloc(
      HCSBB_GAUSS_NNN * HCSBB_NUM_GAUSS_POINTS * sizeof(*global_id_buffer));
  {
    const_coordinate_pointer src_point_coords =
      yac_interp_grid_get_src_field_coords(interp_grid, 0);
    const_yac_int_pointer src_point_global_ids =
      yac_interp_grid_get_src_field_global_ids(interp_grid, 0);
    for (size_t i = 0; i < num_gauss_nnn_patches; ++i)
      compute_gauss_nnn_patch(
        gauss_nnn_patches + i, gauss_points[i],
        nnn_search_results + i * HCSBB_NUM_GAUSS_POINTS * HCSBB_GAUSS_NNN,
        w_nnn_buffer, src_point_buffer, global_id_buffer, src_point_coords,
        src_point_global_ids);
  }
  free(nnn_search_results);
  free(global_id_buffer);
  free(src_point_buffer);
  free(w_nnn_buffer);
  free(gauss_points);
}

static int check_polygon(
  double * check_coord, int num_vertices, size_t const * vertices,
  size_t vertex_start_idx, const_yac_int_pointer global_ids,
  const_coordinate_pointer point_coords,
  size_t * triangle, double * sb_coords) {

  // number of triangle the current source cell can be split into
  size_t num_triangles = (size_t)(num_vertices - 2);
  size_t triangle_indices[num_triangles][3];
  yac_triangulate_cell_indices(
    vertices, num_vertices, vertex_start_idx, &(triangle_indices[0]));

  int ret = 0;
  for (size_t i = 0; (!ret) && (i < num_triangles); ++i) {

    triangle[0] = triangle_indices[i][0];
    triangle[1] = triangle_indices[i][1];
    triangle[2] = triangle_indices[i][2];
    yac_int triangle_ids[3] =
      {global_ids[triangle[0]],
       global_ids[triangle[1]],
       global_ids[triangle[2]]};

    // this is done to ensure bit-reproducibility
    yac_quicksort_index_yac_int_size_t(triangle_ids, 3, triangle);

    double triangle_coords[3][3];
    for (int j = 0; j < 3; ++j)
      memcpy(
        triangle_coords[j], point_coords[triangle[j]], 3 * sizeof(double));

    sb_coords[0] = check_coord[0];
    sb_coords[1] = check_coord[1];
    sb_coords[2] = check_coord[2];
    compute_sb_coords(sb_coords, 1, triangle_coords);

    // if the point is in the current triangle
    ret = (sb_coords[0] > -SB_COORD_TOL) &&
          (sb_coords[1] > -SB_COORD_TOL) &&
          (sb_coords[2] > -SB_COORD_TOL);
  }

  return ret;
}

static int get_matching_vertex_triangle(
  double * tgt_coord, size_t src_cell,
  struct const_basic_grid_data * src_grid_data,
  const_coordinate_pointer src_point_coords,
  size_t * src_triangle_vertices, double * sb_coords) {

  size_t const * src_vertices =
    src_grid_data->cell_to_vertex +
    src_grid_data->cell_to_vertex_offsets[src_cell];
  int num_src_vertices = src_grid_data->num_vertices_per_cell[src_cell];
  yac_int const * vertex_ids = src_grid_data->ids[CORNER];

  // Find start point for triangulation of current source cell. This ensures
  // that the result is decomposition independent.
  size_t lowest_vertex_id_idx = 0;
  {
    yac_int lowest_vertex_id = vertex_ids[src_vertices[0]];
    for (int i = 1; i < num_src_vertices; ++i) {
      if (vertex_ids[src_vertices[i]] < lowest_vertex_id) {
        lowest_vertex_id = vertex_ids[src_vertices[i]];
        lowest_vertex_id_idx = i;
      }
    }
  }

  return
    check_polygon(
      tgt_coord, num_src_vertices, src_vertices, lowest_vertex_id_idx,
      vertex_ids, src_point_coords, src_triangle_vertices, sb_coords);
}

static int get_matching_aux_cell_triangle(
  double * tgt_coord, size_t src_cell,
  struct const_basic_grid_data * src_grid_data,
  const_coordinate_pointer src_point_coords,
  size_t * vertex_to_cell, size_t * vertex_to_cell_offsets,
  int * num_cells_per_vertex, size_t * src_triangle, double * sb_coords) {

  size_t num_vertices = src_grid_data->num_vertices_per_cell[src_cell];
  size_t const * curr_vertices =
    src_grid_data->cell_to_vertex +
    src_grid_data->cell_to_vertex_offsets[src_cell];
  yac_int const * cell_ids = src_grid_data->ids[CELL];

  int ret = 0;

  // for all auxiliary grid cells
  for (size_t i = 0; (!ret) && (i < num_vertices); ++i) {

    size_t curr_vertex = curr_vertices[i];

    int curr_num_src_cells = num_cells_per_vertex[curr_vertex];
    if (curr_num_src_cells == 0) continue;

    size_t * curr_src_cells =
      vertex_to_cell + vertex_to_cell_offsets[curr_vertex];

    // the auxiliary grid cells are already constructed in a way that ensures
    // that they are decomposition independent => start index == 0
    ret =
      check_polygon(
        tgt_coord, curr_num_src_cells, curr_src_cells, 0, cell_ids,
        src_point_coords, src_triangle, sb_coords);
  }

  return ret;
}

static struct tgt_point_search_data * init_tgt_point_data (
  size_t num_tgt_points, size_t * selected_tgt, size_t * tgt_points,
  coordinate_pointer tgt_coords, size_t * src_cells,
  struct interp_grid * interp_grid) {

  enum yac_location src_location =
    yac_interp_grid_get_src_field_location(interp_grid, 0);

  struct tgt_point_search_data * tgt_point_data =
    xmalloc(num_tgt_points * sizeof(*tgt_point_data));

  // if the source data is located at cells, we need to have information on the
  // auxiliary grid (cell points are the vertices of this grid)
  size_t * vertex_to_cell = NULL;
  size_t * vertex_to_cell_offsets = NULL;
  int * num_cells_per_vertex = NULL;
  if (src_location == CELL)
    yac_interp_grid_get_aux_grid_src(
      interp_grid, src_cells, num_tgt_points,
      &vertex_to_cell, &vertex_to_cell_offsets, &num_cells_per_vertex);

  // Remark: you have to get the pointers after the call to
  // yac_dist_grid_get_aux_grid_cells, because the routine might change
  // the underlying memory
  struct const_basic_grid_data * src_grid_data =
    yac_interp_grid_get_basic_grid_data_src(interp_grid);
  const_coordinate_pointer src_point_coords =
    yac_interp_grid_get_src_field_coords(interp_grid, 0);
  const_int_pointer src_mask =
    yac_interp_grid_get_src_field_mask(interp_grid, 0);

  for (size_t i = 0; i < num_tgt_points; ++i, ++tgt_point_data) {

    size_t tgt_idx = selected_tgt[i];

    double * tgt_coord = tgt_coords[tgt_idx];
    tgt_point_data->tgt_point = tgt_points[tgt_idx];
    memcpy(
      tgt_point_data->coordinate_xyz, tgt_coord, 3 * sizeof(*tgt_coord));

  // For the target point we need to find a triangle of source points connected
  // by great circles. The triangulation of the source cells needs to be
  // identical for all source cells on all processes
    size_t src_triangle_vertices[3];
    double sb_coords[3];
    int successfull;
    YAC_ASSERT(
      (src_location == CORNER) || (src_location == CELL),
      "ERROR(init_tgt_point_search_data): invalid source location")
    if (src_location == CORNER) {
      successfull =
        get_matching_vertex_triangle(
          tgt_coord, src_cells[i], src_grid_data, src_point_coords,
          src_triangle_vertices, sb_coords);
    } else {
      successfull =
        get_matching_aux_cell_triangle(
          tgt_coord, src_cells[i], src_grid_data, src_point_coords,
          vertex_to_cell, vertex_to_cell_offsets, num_cells_per_vertex,
          src_triangle_vertices, sb_coords);
    }

    if (successfull) {

      int tmp_type_flag = -1;
      int is_masked = 0;
      for (size_t k = 0; k < 3; ++k) {
        if (fabs(sb_coords[k]) > SB_COORD_TOL) {
          tmp_type_flag++;
          tgt_point_data->src_points[tmp_type_flag] = src_triangle_vertices[k];
          tgt_point_data->sb_coords[tmp_type_flag] = sb_coords[k];
          // check the mask for the source vertices
          if (src_mask != NULL) is_masked |= !src_mask[src_triangle_vertices[k]];
        }
      }
      tgt_point_data->type_flag = (enum interp_type_flag)tmp_type_flag;
      tgt_point_data->is_masked = is_masked;
    } else {
      tgt_point_data->is_masked = 1;
    }
  }
  free(num_cells_per_vertex);
  free(vertex_to_cell_offsets);
  free(vertex_to_cell);

  return tgt_point_data - num_tgt_points;
}

static struct vertex_interp_data * init_vertex_interp_data(
  size_t num_vertices, size_t * vertices,
  struct interp_grid * interp_grid) {

  const_coordinate_pointer src_point_coords =
    yac_interp_grid_get_src_field_coords(interp_grid, 0);

  struct vertex_interp_data * vertex_data =
    xmalloc(num_vertices * sizeof(*vertex_data));

  for (size_t i = 0; i < num_vertices; ++i) {
    size_t vertex_idx = vertices[i];
    memcpy(vertex_data[i].coordinate_xyz,
           src_point_coords[vertex_idx], sizeof(*src_point_coords));
    vertex_data[i].src_point = vertex_idx;
  }

  return vertex_data;
}

static struct edge_interp_data * init_edge_interp_data(
  size_t num_edges, size_t (*edge_to_vertex_d_data)[2],
  struct vertex_interp_data * vertex_data) {

  struct edge_interp_data * edge_data =
    xmalloc(num_edges * sizeof(*edge_data));

  for (size_t i = 0; i < num_edges; ++i) {

    // link edge to derivative data of its two vertices
    edge_data[i].vertex_d_data[0] = vertex_data + edge_to_vertex_d_data[i][0];
    edge_data[i].vertex_d_data[1] = vertex_data + edge_to_vertex_d_data[i][1];

    double * vertex_coordinates_xyz[2] = {
      edge_data[i].vertex_d_data[0]->coordinate_xyz,
      edge_data[i].vertex_d_data[1]->coordinate_xyz};

    // compute middle point of edge
    double * middle_point = edge_data[i].middle_d_data.coordinate_xyz;
    middle_point[0] = vertex_coordinates_xyz[0][0] +
                      vertex_coordinates_xyz[1][0];
    middle_point[1] = vertex_coordinates_xyz[0][1] +
                      vertex_coordinates_xyz[1][1];
    middle_point[2] = vertex_coordinates_xyz[0][2] +
                      vertex_coordinates_xyz[1][2];
    double sb_coord =
      1.0 / sqrt(middle_point[0] * middle_point[0] +
                 middle_point[1] * middle_point[1] +
                 middle_point[2] * middle_point[2]);
    middle_point[0] *= sb_coord;
    middle_point[1] *= sb_coord;
    middle_point[2] *= sb_coord;
    edge_data[i].sb_coord_middle_point = sb_coord;

    // compute vector that is orthogonal to the edge and a tangent on the
    // sphere in the middle point
    crossproduct_d(vertex_coordinates_xyz[0], vertex_coordinates_xyz[1],
                   &(edge_data[i].g[0]));
  }

  return edge_data;
}

static struct triangle_interp_data * init_triangle_interp_data(
  size_t num_triangles, size_t (*triangle_to_edge_data)[3],
  struct edge_interp_data * edge_data) {

  struct triangle_interp_data * triangle_data =
    xmalloc(num_triangles * sizeof(*triangle_data));

  for (size_t i = 0; i < num_triangles; ++i) {

    // link triangle to the data of its edges
    triangle_data[i].edges[0] = edge_data + triangle_to_edge_data[i][0];
    triangle_data[i].edges[1] = edge_data + triangle_to_edge_data[i][1];
    triangle_data[i].edges[2] = edge_data + triangle_to_edge_data[i][2];
  }

  return triangle_data;
}

static size_t do_search_hcsbb (struct interp_method * method,
                               struct interp_grid * interp_grid,
                               size_t * tgt_points, size_t count,
                               struct interp_weights * weights) {

  YAC_ASSERT(
    yac_interp_grid_get_num_src_fields(interp_grid) == 1,
    "ERROR(do_search_hcsbb): invalid number of source fields")

  YAC_ASSERT(
    (yac_interp_grid_get_src_field_location(interp_grid, 0) == CORNER) ||
    (yac_interp_grid_get_src_field_location(interp_grid, 0) == CELL),
    "ERROR(do_search_hcsbb): unsupported source field location type")

  // get coordinates of target points
  coordinate_pointer tgt_coords = xmalloc(count * sizeof(*tgt_coords));
  yac_interp_grid_get_tgt_coordinates(
    interp_grid, tgt_points, count, tgt_coords);

  size_t * size_t_buffer = xmalloc(2 * count * sizeof(*size_t_buffer));
  size_t * src_cells = size_t_buffer;
  size_t * reorder_idx = size_t_buffer + count;

  // get matching source cells for all target points
  // (this search assumes that all edges of the grid are great circles)
  yac_interp_grid_do_points_search_gc(
    interp_grid, tgt_coords, count, src_cells);

  // sort target points, for which we found a source cell, to the beginning of
  // the array
  for (size_t i = 0; i < count; ++i) reorder_idx[i] = i;
  yac_quicksort_index_size_t_size_t(src_cells, count, reorder_idx);

  // count the number of target for which a valid source cell was found
  size_t temp_result_count = 0;
  for (temp_result_count = 0; temp_result_count < count; ++temp_result_count)
    if (src_cells[temp_result_count] == SIZE_MAX) break;

  // initialise search data for all target point, for which a matching source
  // cell was found
  struct tgt_point_search_data * tgt_point_data =
    init_tgt_point_data(
      temp_result_count, reorder_idx, tgt_points, tgt_coords, src_cells,
      interp_grid);
  free(tgt_coords);

  // sort target point data (first by whether they are marked or not, second
  // by interp_type, and third by vertex indices)
  qsort(tgt_point_data, temp_result_count, sizeof(*tgt_point_data),
        compare_tgt_point_search_data);

  // move target points which cannot be interpolated to the end of the list
  for (size_t i = temp_result_count; i < count; ++i)
    reorder_idx[i] = tgt_points[reorder_idx[i]];
  memcpy(tgt_points + temp_result_count, reorder_idx + temp_result_count,
         (count - temp_result_count) * sizeof(*tgt_points));
  // put remaining target point to the front of the list (can also contain
  // failed target points)
  for (size_t i = 0; i < temp_result_count; ++i)
    tgt_points[i] = tgt_point_data[i].tgt_point;
  free(size_t_buffer);

  // count the number of target points that can be interpolated
  size_t result_count = 0;
  size_t result_counts[3] = {0,0,0};
  for (result_count = 0; result_count < temp_result_count; ++result_count) {
    if (tgt_point_data[result_count].is_masked) break;
    else result_counts[(int)(tgt_point_data[result_count].type_flag)]++;
  }

  // get all edges and triangles required for the interpolation
  // (each triangle requires three edges, usually each triangle edges is shared
  //  between two triangles)
  size_t * unique_vertices;
  size_t (*unique_edge_to_unique_vertex)[2];
  size_t (*unique_triangle_to_unique_edge)[3];
  size_t num_unique_vertices, num_unique_edges, num_unique_triangles;
  get_unique_interp_data(
    tgt_point_data + result_counts[0], result_counts[1], result_counts[2],
    &unique_vertices, &num_unique_vertices,
    &unique_edge_to_unique_vertex, &num_unique_edges,
    &unique_triangle_to_unique_edge, &num_unique_triangles);

  // initialise data for all unique vertices, edges and triangles
  struct vertex_interp_data * vertex_data =
    init_vertex_interp_data(
      num_unique_vertices, unique_vertices, interp_grid);
  struct edge_interp_data * edge_data =
    init_edge_interp_data(
      num_unique_edges, unique_edge_to_unique_vertex, vertex_data);
  struct triangle_interp_data * triangle_data =
    init_triangle_interp_data(
      num_unique_triangles, unique_triangle_to_unique_edge, edge_data);
  free(unique_vertices);
  free(unique_edge_to_unique_vertex);
  free(unique_triangle_to_unique_edge);

  // link tgt_point_data to the edge_data and triangle_data
  for (size_t i = result_counts[0];
       i < result_counts[0] + result_counts[1]; ++i)
    tgt_point_data[i].interp_data.edge =
      edge_data + tgt_point_data[i].interp_data.idx;
  for (size_t i = result_counts[0] + result_counts[1];
       i < result_counts[0] + result_counts[1] + result_counts[2]; ++i)
    tgt_point_data[i].interp_data.triangle =
      triangle_data + tgt_point_data[i].interp_data.idx;

  // for the interpolation we will need to estimate derivatives at all required
  // vertices and the middle points of all edges
  size_t num_gauss_nnn_patches;
  struct gauss_nnn_patch * gauss_nnn_patches; // have to free this later
  generate_derivative_data(
    num_unique_vertices, vertex_data, num_unique_edges, edge_data, interp_grid,
    &num_gauss_nnn_patches, &gauss_nnn_patches);

  // compute the actual weights
  struct weight_vector_data * weight_vector;
  size_t * num_weights_per_tgt, total_num_weights;
  compute_weights(tgt_point_data, result_count,
                  edge_data, num_unique_edges,
                  triangle_data, num_unique_triangles,
                  &weight_vector, &num_weights_per_tgt, &total_num_weights);
  free(vertex_data);
  if (num_unique_edges > 0) free(edge_data->c[0].data);
  free(edge_data);
  for (size_t i = 0; i < num_unique_triangles; ++i)
    for (size_t j = 0; j < 3; ++j) free(triangle_data[i].c_111_a[j].data);
  free(triangle_data);
  free(tgt_point_data);

  // free memory associated with the gauss_nnn_patches
  free_gauss_nnn_patches(gauss_nnn_patches, num_gauss_nnn_patches);

  size_t * src_points = xmalloc(total_num_weights * sizeof(*src_points));
  double * weight_data = xmalloc(total_num_weights * sizeof(*weight_data));
  for (size_t i = 0; i < total_num_weights; ++i) {
    weight_data[i] = weight_vector[i].weight;
    src_points[i] = weight_vector[i].point_id;
  }
  free(weight_vector);

  struct remote_point * srcs =
    yac_interp_grid_get_src_remote_points(
      interp_grid, 0, src_points, total_num_weights);

  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(
        interp_grid, tgt_points, result_count),
    .count = result_count};

  // store weights
  yac_interp_weights_add_wsum(
    weights, &tgts, num_weights_per_tgt, srcs, weight_data);

  free(tgts.data);
  free(srcs);
  free(weight_data);
  free(src_points);
  free(num_weights_per_tgt);

  return result_count;
}

struct interp_method * yac_interp_method_hcsbb_new() {

  struct interp_method_hcsbb * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_hcsbb_vtable;

  return (struct interp_method*)method;
}

static void delete_hcsbb(struct interp_method * method) {
  free(method);
}
