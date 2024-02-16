/**
 * @file interp_stack_config.h
 * @brief Structs and interfaces to defined interpolation stack configurations
 *
 * @copyright Copyright  (C)  2022 Moritz Hanke <hanke@dkrz.de>
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

#ifndef INTERP_STACK_CONFIG_H
#define INTERP_STACK_CONFIG_H

#include "interp_method_avg.h"
#include "interp_method_check.h"
#include "interp_method_conserv.h"
#include "interp_method_creep.h"
#include "interp_method_nnn.h"
#include "interp_method_fixed.h"
#include "interp_method_file.h"
#include "interp_method_hcsbb.h"
#include "interp_method_spmap.h"
#include "interp_method_callback.h"

#define MAX_ROUTINE_NAME_LENGTH (256)
#define MAX_FILE_NAME_LENGTH (512)

enum yac_interpolation_list {
   UNDEFINED             = 0,
   AVERAGE               = 1,
   N_NEAREST_NEIGHBOR    = 2,
   CONSERVATIVE          = 3,
   SOURCE_TO_TARGET_MAP  = 4,
   FIXED_VALUE           = 5,
   USER_FILE             = 6,
   CHECK                 = 7,
   BERNSTEIN_BEZIER      = 8,
   RADIAL_BASIS_FUNCTION = 9,
   CREEP                 = 10,
   USER_CALLBACK         = 11,
};

/** \example test_interp_stack_config.c
 * Tests for interpolation stack interface routines.
 */

struct yac_interp_stack_config;
union yac_interp_stack_config_entry;

struct yac_interp_stack_config * yac_interp_stack_config_new();
void yac_interp_stack_config_delete(
  struct yac_interp_stack_config * interp_stack_config);
struct yac_interp_stack_config * yac_interp_stack_config_copy(
  struct yac_interp_stack_config * interp_stack);

void yac_interp_stack_config_add_average(
  struct yac_interp_stack_config * interp_stack_config,
  enum yac_interp_avg_weight_type reduction_type, int partial_coverage);
void yac_interp_stack_config_add_nnn(
  struct yac_interp_stack_config * interp_stack_config,
  enum yac_interp_nnn_weight_type type, size_t n, double scale);
void yac_interp_stack_config_add_conservative(
  struct yac_interp_stack_config * interp_stack_config,
  int order, int enforced_conserv, int partial_coverage,
  enum yac_interp_method_conserv_normalisation normalisation);
void yac_interp_stack_config_add_spmap(
  struct yac_interp_stack_config * interp_stack_config,
  double spread_distance, double max_search_distance,
  enum yac_interp_spmap_weight_type weight_type);
void yac_interp_stack_config_add_hcsbb(
  struct yac_interp_stack_config * interp_stack_config);
void yac_interp_stack_config_add_user_file(
  struct yac_interp_stack_config * interp_stack_config,
  char const * filename, char const * src_grid_name,
  char const * tgt_grid_name);
void yac_interp_stack_config_add_fixed(
  struct yac_interp_stack_config * interp_stack_config, double value);
void yac_interp_stack_config_add_check(
  struct yac_interp_stack_config * interp_stack_config,
  char const * constructor_key, char const * do_search_key);
void yac_interp_stack_config_add_creep(
  struct yac_interp_stack_config * interp_stack_config, int creep_distance);
void yac_interp_stack_config_add_user_callback(
  struct yac_interp_stack_config * interp_stack_config,
  char const * func_compute_weights_key);

int yac_interp_stack_config_compare(void const * a, void const * b);
struct interp_method ** yac_interp_stack_config_generate(
  struct yac_interp_stack_config * interp_stack);

size_t yac_interp_stack_config_get_pack_size(
  struct yac_interp_stack_config * interp_stack, MPI_Comm comm);
void yac_interp_stack_config_pack(
  struct yac_interp_stack_config * interp_stack,
  void * buffer, int buffer_size, int * position, MPI_Comm comm);
struct yac_interp_stack_config * yac_interp_stack_config_unpack(
  void * buffer, int buffer_size, int * position, MPI_Comm comm);

size_t yac_interp_stack_config_get_size(
  struct yac_interp_stack_config * interp_stack);
union yac_interp_stack_config_entry const *
  yac_interp_stack_config_get_entry(
    struct yac_interp_stack_config * interp_stack,
    size_t interp_stack_idx);

enum yac_interpolation_list yac_interp_stack_config_entry_get_type(
  union yac_interp_stack_config_entry const * interp_stack_entry);

void yac_interp_stack_config_entry_get_average(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  enum yac_interp_avg_weight_type * reduction_type,
  int * partial_coverage);
void yac_interp_stack_config_entry_get_nnn(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  enum yac_interp_nnn_weight_type * type, size_t * n, double * scale);
void yac_interp_stack_config_entry_get_conservative(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  int * order, int * enforced_conserv, int * partial_coverage,
  enum yac_interp_method_conserv_normalisation * normalisation);
void yac_interp_stack_config_entry_get_spmap(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  double * spread_distance, double * max_search_distance,
  enum yac_interp_spmap_weight_type * weight_type);
void yac_interp_stack_config_entry_get_user_file(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** filename, char const ** src_grid_name,
  char const ** tgt_grid_name);
void yac_interp_stack_config_entry_get_fixed(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  double * value);
void yac_interp_stack_config_entry_get_check(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** constructor_key, char const ** do_search_key);
void yac_interp_stack_config_entry_get_creep(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  int * creep_distance);
void yac_interp_stack_config_entry_get_user_callback(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** func_compute_weights_key);

#endif // INTERP_STACK_CONFIG_H
