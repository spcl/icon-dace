/**
 * @file test_interp_method_parallel.c
 *
 * @copyright Copyright  (C)  2019 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *
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
#include <stdlib.h>

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "read_icon_grid.h"
#include "interp_method.h"
#include "interp_method_fixed.h"
#include "yac_mpi.h"

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

int main(void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  set_even_io_rank_list(MPI_COMM_WORLD);

  char * filenames[] ={"icon_grid_R02B02.nc", "icon_grid_R02B03_G.nc"};
  struct basic_grid_data grid_data[2];

  for (int i = 0; i < 2; ++i)
    grid_data[i] =
      read_icon_grid_information_parallel2(
        filenames[i], MPI_COMM_WORLD);

  struct yac_basic_grid * grids[2] =
    {yac_basic_grid_new(filenames[0], grid_data[0]),
     yac_basic_grid_new(filenames[1], grid_data[1])};

  struct dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

  struct interp_field src_fields[] =
    {{.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct interp_field tgt_field =
    {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, filenames[0], filenames[1],
                        num_src_fields, src_fields, tgt_field);

  struct interp_method * method_stack[] =
    {yac_interp_method_fixed_new(-1.0), NULL};

  struct interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  struct interpolation * interpolation_src =
    yac_interp_weights_get_interpolation(
      weights, MAPPING_ON_SRC, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
  struct interpolation * interpolation_tgt =
    yac_interp_weights_get_interpolation(
      weights, MAPPING_ON_TGT, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

  yac_interpolation_delete(interpolation_src);
  yac_interpolation_delete(interpolation_tgt);

  yac_interp_weights_delete(weights);

  yac_interp_method_delete(method_stack);

  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  for (int i = 0; i < 2; ++i) yac_basic_grid_delete(grids[i]);

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
