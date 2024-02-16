/**
 * @file test_common.h
 *
 * @copyright Copyright  (C)  2018 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
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
#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include "grid_cell.h"
#include "utils.h"

/** \example test_common.c
 * A test for test-internal functions - not relevant for the external user.
 */

#define TO_POINTER(a) (to_pointer(&a, sizeof(a)))

struct grid_cell generate_cell_deg(
  double * lon, double * lat, enum yac_edge_type * edge_type,
  size_t num_corners);
struct grid_cell generate_cell_rad(
  double * lon, double * lat, enum yac_edge_type * edge_type,
  size_t num_corners);
struct grid_cell generate_cell_3d(
  coordinate_pointer coords, enum yac_edge_type * edge_type,
  size_t num_corners);
struct grid_cell generate_cell(
  double * lon, double * lat, enum yac_edge_type * edge_type,
  size_t num_corners);
int intersect(enum yac_edge_type edge_type_a,
              double lon_a, double lat_a, double lon_b, double lat_b,
              enum yac_edge_type edge_type_b,
              double lon_c, double lat_c, double lon_d, double lat_d,
              double * intersection);
void * to_pointer(void * data, size_t size_data);
int double_are_equal(double a, double b);
int double_are_unequal(double a, double b);
void set_even_io_rank_list(MPI_Comm comm);
void clear_yac_io_env();

#endif // TEST_COMMON_H

