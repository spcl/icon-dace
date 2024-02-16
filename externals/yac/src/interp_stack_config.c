/**
 * @file interp_stack_config.c
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>

#include "yac_mpi.h"
#include "interp_stack_config.h"

union yac_interp_stack_config_entry {
  // - every interpolation type needs to have the entry type as first
  //   entry in order to be able to compare them using memcmp
  // - no pointers are allowed here in order to be comparable
  struct {
    enum yac_interpolation_list type;
  } general;
  struct {
    enum yac_interpolation_list type;
    enum yac_interp_avg_weight_type reduction_type;
    int partial_coverage;
  } average;
  struct {
    enum yac_interpolation_list type;
    struct yac_nnn_config config;
  } n_nearest_neighbor,
    radial_basis_function;
  struct {
    enum yac_interpolation_list type;
    int order;
    int enforced_conserv;
    int partial_coverage;
    enum yac_interp_method_conserv_normalisation normalisation;
  } conservative;
  struct {
    enum yac_interpolation_list type;
    double spread_distance;
    double max_search_distance;
    enum yac_interp_spmap_weight_type weight_type;
  } spmap;
  struct {
    enum yac_interpolation_list type;
  } hcsbb;
  struct {
    enum yac_interpolation_list type;
    char filename[MAX_FILE_NAME_LENGTH];
    char src_grid_name[MAX_FILE_NAME_LENGTH];
    char tgt_grid_name[MAX_FILE_NAME_LENGTH];
  } user_file;
  struct {
    enum yac_interpolation_list type;
    double value;
  } fixed;
  struct {
    enum yac_interpolation_list type;
    char constructor_key[MAX_ROUTINE_NAME_LENGTH];
    char do_search_key[MAX_ROUTINE_NAME_LENGTH];
  } check;
  struct {
    enum yac_interpolation_list type;
    int creep_distance;
  } creep;
  struct {
    enum yac_interpolation_list type;
    char func_compute_weights_key[MAX_ROUTINE_NAME_LENGTH];
  } user_callback;
};

struct yac_interp_stack_config {
  union yac_interp_stack_config_entry * config;
  size_t size;
};

struct yac_interp_stack_config * yac_interp_stack_config_copy(
  struct yac_interp_stack_config * interp_stack) {

  struct yac_interp_stack_config * interp_stack_copy =
    xmalloc(1 * sizeof(*interp_stack_copy));
  interp_stack_copy->size = interp_stack->size;
  interp_stack_copy->config =
    xmalloc(interp_stack_copy->size * sizeof(*(interp_stack_copy->config)));
  memcpy(
    interp_stack_copy->config, interp_stack->config,
    interp_stack_copy->size * sizeof(*(interp_stack_copy->config)));
  return interp_stack_copy;
}

static void inline check_interpolation_type(
  enum yac_interpolation_list type, char const * routine) {

  YAC_ASSERT_F(
    (type == AVERAGE) ||
    (type == RADIAL_BASIS_FUNCTION) ||
    (type == N_NEAREST_NEIGHBOR) ||
    (type == CONSERVATIVE) ||
    (type == SOURCE_TO_TARGET_MAP) ||
    (type == FIXED_VALUE) ||
    (type == BERNSTEIN_BEZIER) ||
    (type == USER_FILE) ||
    (type == CHECK) ||
    (type == CREEP) ||
    (type == USER_CALLBACK) ||
    (type == UNDEFINED),
    "ERROR(%s): invalid interpolation type", routine)
}

static int yac_interp_stack_config_entry_compare(
  void const * a_, void const * b_) {

  union yac_interp_stack_config_entry const * a =
   (union yac_interp_stack_config_entry const *)a_;
  union yac_interp_stack_config_entry const * b =
   (union yac_interp_stack_config_entry const *)b_;

  check_interpolation_type(
    a->general.type, "yac_interp_stack_config_entry_compare");
  check_interpolation_type(
    b->general.type, "yac_interp_stack_config_entry_compare");

  if (a->general.type != b->general.type)
    return (a->general.type > b->general.type) -
           (a->general.type < b->general.type);

  YAC_ASSERT(
    (a->general.type == AVERAGE) ||
    (a->general.type == N_NEAREST_NEIGHBOR) ||
    (a->general.type == CONSERVATIVE) ||
    (a->general.type == SOURCE_TO_TARGET_MAP) ||
    (a->general.type == FIXED_VALUE) ||
    (a->general.type == USER_FILE) ||
    (a->general.type == CHECK) ||
    (a->general.type == BERNSTEIN_BEZIER) ||
    (a->general.type == RADIAL_BASIS_FUNCTION) ||
    (a->general.type == CREEP) ||
    (a->general.type == USER_CALLBACK),
    "ERROR(yac_interp_stack_config_entry_compare): invalid interpolation type")

  switch(a->general.type) {
    default:
    case (BERNSTEIN_BEZIER):
      return 0;
    case (AVERAGE): {
      if (a->average.reduction_type != b->average.reduction_type)
        return (a->average.reduction_type > b->average.reduction_type) -
               (a->average.reduction_type < b->average.reduction_type);
      return (a->average.partial_coverage > b->average.partial_coverage) -
             (a->average.partial_coverage < b->average.partial_coverage);
    }
    case (RADIAL_BASIS_FUNCTION):
    case (N_NEAREST_NEIGHBOR): {
      if (a->n_nearest_neighbor.config.n !=
          b->n_nearest_neighbor.config.n)
        return (a->n_nearest_neighbor.config.n >
                b->n_nearest_neighbor.config.n) -
               (a->n_nearest_neighbor.config.n <
                b->n_nearest_neighbor.config.n);
      if (a->n_nearest_neighbor.config.type !=
          b->n_nearest_neighbor.config.type)
        return (a->n_nearest_neighbor.config.type >
                b->n_nearest_neighbor.config.type) -
               (a->n_nearest_neighbor.config.type <
                b->n_nearest_neighbor.config.type);
      switch (a->n_nearest_neighbor.config.type) {
        case (NNN_GAUSS):
          return (a->n_nearest_neighbor.config.data.gauss_scale >
                  b->n_nearest_neighbor.config.data.gauss_scale) -
                 (a->n_nearest_neighbor.config.data.gauss_scale <
                  b->n_nearest_neighbor.config.data.gauss_scale);
        case (NNN_RBF):
          return (a->n_nearest_neighbor.config.data.rbf_scale >
                  b->n_nearest_neighbor.config.data.rbf_scale) -
                 (a->n_nearest_neighbor.config.data.rbf_scale <
                  b->n_nearest_neighbor.config.data.rbf_scale);
        default:
          return 0;
      };
    }
    case (CONSERVATIVE): {
      if (a->conservative.order !=
          b->conservative.order)
        return (a->conservative.order >
                b->conservative.order) -
               (a->conservative.order <
                b->conservative.order);
      if (a->conservative.enforced_conserv !=
          b->conservative.enforced_conserv)
        return (a->conservative.enforced_conserv >
                b->conservative.enforced_conserv) -
               (a->conservative.enforced_conserv <
                b->conservative.enforced_conserv);
      if (a->conservative.partial_coverage !=
          b->conservative.partial_coverage)
        return (a->conservative.partial_coverage >
                b->conservative.partial_coverage) -
               (a->conservative.partial_coverage <
                b->conservative.partial_coverage);
      return (a->conservative.normalisation >
              b->conservative.normalisation) -
              (a->conservative.normalisation <
              b->conservative.normalisation);
    }
    case (SOURCE_TO_TARGET_MAP): {
      if (fabs(a->spmap.spread_distance -
               b->spmap.spread_distance) > yac_angle_tol)
        return (a->spmap.spread_distance >
                b->spmap.spread_distance) -
               (a->spmap.spread_distance <
                b->spmap.spread_distance);
      if (fabs(a->spmap.max_search_distance -
               b->spmap.max_search_distance) > yac_angle_tol)
        return (a->spmap.max_search_distance >
                b->spmap.max_search_distance) -
               (a->spmap.max_search_distance <
                b->spmap.max_search_distance);
      return (a->spmap.weight_type >
              b->spmap.weight_type) -
             (a->spmap.weight_type <
              b->spmap.weight_type);
    }
    case (FIXED_VALUE): {
      return (a->fixed.value >
              b->fixed.value) -
             (a->fixed.value <
              b->fixed.value);
    }
    case (USER_FILE): {
      int ret;
      if ((ret = strcmp(a->user_file.filename, b->user_file.filename)))
        return ret;
      if ((ret = strcmp(a->user_file.src_grid_name, b->user_file.src_grid_name)))
        return ret;
      return strcmp(a->user_file.tgt_grid_name, b->user_file.tgt_grid_name);
    }
    case (CHECK): {
      int ret;
      if ((ret = strncmp(a->check.constructor_key, b->check.constructor_key,
                         sizeof(a->check.constructor_key)))) return ret;
      return
        strncmp(
          a->check.do_search_key, b->check.do_search_key,
          sizeof(a->check.do_search_key));
    };
    case (CREEP): {
      return
        (a->creep.creep_distance > b->creep.creep_distance) -
        (a->creep.creep_distance < b->creep.creep_distance);
    }
    case (USER_CALLBACK): {
      return strcmp(a->user_callback.func_compute_weights_key,
                    b->user_callback.func_compute_weights_key);
    }
  };
}

int yac_interp_stack_config_compare(void const * a_, void const * b_) {

  struct yac_interp_stack_config const * a =
   (struct yac_interp_stack_config const *)a_;
  struct yac_interp_stack_config const * b =
   (struct yac_interp_stack_config const *)b_;

  int ret;
  if ((ret = (a->size > b->size) - (a->size < b->size))) return ret;

  size_t stack_size = a->size;
  for (size_t method_idx = 0; method_idx < stack_size; ++method_idx)
    if ((ret =
           yac_interp_stack_config_entry_compare(
             &(a->config[method_idx]), &(b->config[method_idx]))))
      return ret;
  return 0;
}

struct interp_method ** yac_interp_stack_config_generate(
  struct yac_interp_stack_config * interp_stack) {

  size_t interp_stack_size = interp_stack->size;
  struct interp_method ** method_stack =
    xmalloc((interp_stack_size + 1) * sizeof(*method_stack));
  method_stack[interp_stack_size] = NULL;

  for (size_t i = 0; i < interp_stack_size; ++i) {

    check_interpolation_type(
      interp_stack->config[i].general.type, "yac_interp_stack_config_generate");
    YAC_ASSERT(
      interp_stack->config[i].general.type != UNDEFINED,
      "ERROR(yac_interp_stack_config_generate): "
      "unsupported interpolation method")
    switch((int)(interp_stack->config[i].general.type)) {
      default:
      case(AVERAGE): {
        enum yac_interp_avg_weight_type weight_type =
          interp_stack->config[i].average.reduction_type;
        int partial_coverage =
          (int)interp_stack->config[i].average.partial_coverage;
        method_stack[i] =
          yac_interp_method_avg_new(weight_type, partial_coverage);
        break;
      }
      case(CONSERVATIVE): {
        int order =
          interp_stack->config[i].conservative.order;
        int enforced_conserv =
          interp_stack->config[i].conservative.enforced_conserv;
        int partial_coverage =
          interp_stack->config[i].conservative.partial_coverage;
        enum yac_interp_method_conserv_normalisation normalisation =
          interp_stack->config[i].conservative.normalisation;

        method_stack[i] =
          yac_interp_method_conserv_new(
            order, enforced_conserv, partial_coverage, normalisation);
        break;
      }
      case(FIXED_VALUE): {
        double fixed_value =
          interp_stack->config[i].fixed.value;
        method_stack[i] = yac_interp_method_fixed_new(fixed_value);
        break;
      }
      case(USER_FILE): {
        char const * weight_file_name =
          interp_stack->config[i].user_file.filename;
        char const * src_grid_name =
          interp_stack->config[i].user_file.src_grid_name;
        char const * tgt_grid_name =
          interp_stack->config[i].user_file.tgt_grid_name;
        method_stack[i] =
          yac_interp_method_file_new(
            weight_file_name, src_grid_name, tgt_grid_name);
        break;
      }
      case(CHECK): {
        func_constructor constructor_callback;
        void * constructor_user_data;
        func_do_search do_search_callback;
        void * do_search_user_data;

        yac_interp_method_check_get_constructor_callback(
          interp_stack->config[i].check.constructor_key,
          &constructor_callback, &constructor_user_data);
        yac_interp_method_check_get_do_search_callback(
          interp_stack->config[i].check.do_search_key,
          &do_search_callback, &do_search_user_data);

        method_stack[i] =
          yac_interp_method_check_new(
            constructor_callback, constructor_user_data,
            do_search_callback, do_search_user_data);
        break;
      }
      case (N_NEAREST_NEIGHBOR): {
        method_stack[i] =
          yac_interp_method_nnn_new(
            interp_stack->config[i].
              n_nearest_neighbor.config);
        break;
      }
      case (BERNSTEIN_BEZIER): {
        method_stack[i] =
          yac_interp_method_hcsbb_new();
        break;
      }
      case (RADIAL_BASIS_FUNCTION): {
        method_stack[i] =
          yac_interp_method_nnn_new(
            interp_stack->config[i].
              radial_basis_function.config);
        break;
      }
      case (SOURCE_TO_TARGET_MAP): {
        method_stack[i] =
          yac_interp_method_spmap_new(
            interp_stack->config[i].spmap.spread_distance,
            interp_stack->config[i].spmap.max_search_distance,
            interp_stack->config[i].spmap.weight_type);
        break;
      }
      case (CREEP): {
        method_stack[i] =
          yac_interp_method_creep_new(
            interp_stack->config[i].creep.creep_distance);
        break;
      }
      case(USER_CALLBACK): {
        yac_func_compute_weights compute_weights_callback;
        void * user_data;
        yac_interp_method_callback_get_compute_weights_callback(
          interp_stack->config[i].user_callback.func_compute_weights_key,
          &compute_weights_callback, &user_data);

        method_stack[i] =
          yac_interp_method_callback_new(compute_weights_callback, user_data);
        break;
      }
    };
  }

  return method_stack;
}

static size_t yac_interp_stack_config_get_string_pack_size(
  char const * string, MPI_Comm comm) {

  int strlen_pack_size, string_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &strlen_pack_size), comm);

  YAC_ASSERT(
    string != NULL, "ERROR(yac_interp_stack_config_get_string_pack_size): "
    "string is NULL");

  yac_mpi_call(
    MPI_Pack_size(
      (int)(strlen(string)), MPI_CHAR, comm, &string_pack_size), comm);

  return (size_t)strlen_pack_size + (size_t)string_pack_size;
}

static size_t yac_interp_stack_config_get_entry_pack_size(
  union yac_interp_stack_config_entry * entry, MPI_Comm comm) {

  int int_pack_size, dbl_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &int_pack_size), comm);
  yac_mpi_call(MPI_Pack_size(1, MPI_DOUBLE, comm, &dbl_pack_size), comm);

  check_interpolation_type(
    entry->general.type,
    "yac_interp_stack_config_get_entry_pack_size");
  YAC_ASSERT(
    entry->general.type != UNDEFINED,
    "ERROR(yac_interp_stack_config_get_entry_pack_size): "
    "invalid interpolation type")
  switch (entry->general.type) {
    default:
    case (AVERAGE):
      return (size_t)int_pack_size + // type
             (size_t)int_pack_size + // reduction_type
             (size_t)int_pack_size;  // partial_coverage
    case (RADIAL_BASIS_FUNCTION):
    case (N_NEAREST_NEIGHBOR):
      return (size_t)int_pack_size + // type
             (size_t)int_pack_size + // weight_type
             (size_t)int_pack_size + // n
             (size_t)dbl_pack_size;  // scale
    case (CONSERVATIVE):
      return (size_t)int_pack_size + // type
             (size_t)int_pack_size + // order
             (size_t)int_pack_size + // enforced_conserv
             (size_t)int_pack_size + // partial_coverage
             (size_t)int_pack_size;  // normalisation
    case (SOURCE_TO_TARGET_MAP):
      return (size_t)int_pack_size + // type
             (size_t)dbl_pack_size + // spread_distance
             (size_t)dbl_pack_size + // max_search_distance
             (size_t)int_pack_size;  // weight_type
    case (FIXED_VALUE):
      return (size_t)int_pack_size + // type
             (size_t)dbl_pack_size;  // value
    case (BERNSTEIN_BEZIER):
      return (size_t)int_pack_size;  // type
    case (USER_FILE):
      return (size_t)int_pack_size + // type
             yac_interp_stack_config_get_string_pack_size(
               entry->user_file.filename, comm) +
             yac_interp_stack_config_get_string_pack_size(
               entry->user_file.src_grid_name, comm) +
             yac_interp_stack_config_get_string_pack_size(
               entry->user_file.tgt_grid_name, comm);
    case (CHECK):
      return (size_t)int_pack_size + // type
             yac_interp_stack_config_get_string_pack_size(
               entry->check.constructor_key, comm) + // constructor_key
             yac_interp_stack_config_get_string_pack_size(
               entry->check.do_search_key, comm);    // do_search_key
    case (CREEP):
      return (size_t)int_pack_size + // type
             (size_t)int_pack_size;  // creep_distance
    case (USER_CALLBACK):
      return (size_t)int_pack_size + // type
             yac_interp_stack_config_get_string_pack_size(
               entry->user_callback.func_compute_weights_key, comm);
               // func_compute_weights_key
  }
}

size_t yac_interp_stack_config_get_pack_size(
  struct yac_interp_stack_config * interp_stack, MPI_Comm comm) {

  int size_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &size_pack_size), comm);

  size_t config_pack_size = 0;

  for (size_t i = 0; i < interp_stack->size; ++i)
    config_pack_size +=
      yac_interp_stack_config_get_entry_pack_size(
        interp_stack->config + i, comm);

  return (size_t)size_pack_size + config_pack_size;
}

static void yac_interp_stack_config_pack_string(
  char const * string, void * buffer, int buffer_size, int * position,
  MPI_Comm comm) {

  size_t len = (string == NULL)?0:strlen(string);

  YAC_ASSERT(
    len <= INT_MAX, "ERROR(yac_interp_stack_config_pack_string): string too long")

  int len_int = (int)len;

  yac_mpi_call(
    MPI_Pack(
      &len_int, 1, MPI_INT, buffer, buffer_size, position, comm), comm);

  if (len > 0)
    yac_mpi_call(
      MPI_Pack(
        string, len_int, MPI_CHAR, buffer, buffer_size, position, comm),
      comm);
}

static void yac_interp_stack_config_pack_entry(
  union yac_interp_stack_config_entry * entry,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  int type = (int)(entry->general.type);
  yac_mpi_call(
    MPI_Pack(&type, 1, MPI_INT, buffer, buffer_size, position, comm), comm);

  check_interpolation_type(
    entry->general.type,
    "yac_interp_stack_config_pack_entry");
  YAC_ASSERT(
    entry->general.type != UNDEFINED,
    "ERROR(yac_interp_stack_config_pack_entry): "
    "invalid interpolation type")
  switch (entry->general.type) {
    default:
    case (AVERAGE): {
      int reduction_type = (int)(entry->average.reduction_type);
      yac_mpi_call(
        MPI_Pack(
          &reduction_type, 1, MPI_INT, buffer, buffer_size, position, comm),
        comm);
      yac_mpi_call(
        MPI_Pack(
          &(entry->average.partial_coverage), 1, MPI_INT,
          buffer, buffer_size, position, comm), comm);
      break;
    }
    case (RADIAL_BASIS_FUNCTION):
    case (N_NEAREST_NEIGHBOR): {
      YAC_ASSERT(
        entry->n_nearest_neighbor.config.n <= INT_MAX,
        "ERROR(yac_interp_stack_config_pack_entry): "
        "n_nearest_neighbor.config.n bigger than INT_MAX")
      int type = (int)(entry->n_nearest_neighbor.config.type);
      yac_mpi_call(
        MPI_Pack(
          &type, 1, MPI_INT, buffer, buffer_size, position, comm), comm);
      int n = (int)(entry->n_nearest_neighbor.config.n);
      yac_mpi_call(
        MPI_Pack(
          &n, 1, MPI_INT, buffer, buffer_size, position, comm), comm);
      yac_mpi_call(
        MPI_Pack(
          &(entry->n_nearest_neighbor.config.data.rbf_scale), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
      break;
    }
    case (CONSERVATIVE): {
      yac_mpi_call(
        MPI_Pack(
          &(entry->conservative.order), 1, MPI_INT,
          buffer, buffer_size, position, comm), comm);
      yac_mpi_call(
        MPI_Pack(
          &(entry->conservative.enforced_conserv), 1, MPI_INT,
          buffer, buffer_size, position, comm), comm);
      yac_mpi_call(
        MPI_Pack(
          &(entry->conservative.partial_coverage), 1, MPI_INT,
          buffer, buffer_size, position, comm), comm);
      int normalisation = (int)(entry->conservative.normalisation);
      yac_mpi_call(
        MPI_Pack(
          &normalisation, 1, MPI_INT, buffer, buffer_size, position, comm),
        comm);
      break;
    }
    case (SOURCE_TO_TARGET_MAP): {
      yac_mpi_call(
        MPI_Pack(
          &(entry->spmap.spread_distance), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
      yac_mpi_call(
        MPI_Pack(
          &(entry->spmap.max_search_distance), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
      int weight_type = (int)(entry->spmap.weight_type);
      yac_mpi_call(
        MPI_Pack(
          &weight_type, 1, MPI_INT, buffer, buffer_size, position, comm),
        comm);
      break;
    }
    case (FIXED_VALUE): {
      yac_mpi_call(
        MPI_Pack(
          &(entry->fixed.value), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
      break;
    }
    case (BERNSTEIN_BEZIER):
      break;
    case (USER_FILE): {
      yac_interp_stack_config_pack_string(
        entry->user_file.filename, buffer, buffer_size, position, comm);
      yac_interp_stack_config_pack_string(
        entry->user_file.src_grid_name, buffer, buffer_size, position, comm);
      yac_interp_stack_config_pack_string(
        entry->user_file.tgt_grid_name, buffer, buffer_size, position, comm);
      break;
    }
    case (CHECK): {
      yac_interp_stack_config_pack_string(
        entry->check.constructor_key, buffer, buffer_size, position, comm);
      yac_interp_stack_config_pack_string(
        entry->check.do_search_key, buffer, buffer_size, position, comm);
      break;
    }
    case (CREEP): {
      yac_mpi_call(
        MPI_Pack(
          &(entry->creep.creep_distance), 1, MPI_INT,
          buffer, buffer_size, position, comm), comm);
      break;
    }
    case (USER_CALLBACK): {
      yac_interp_stack_config_pack_string(
        entry->user_callback.func_compute_weights_key,
        buffer, buffer_size, position, comm);
      break;
    }
  }
}

void yac_interp_stack_config_pack(
  struct yac_interp_stack_config * interp_stack,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  int stack_size = (int)(interp_stack->size);
  yac_mpi_call(
    MPI_Pack(
      &stack_size, 1, MPI_INT,
      buffer, buffer_size, position, comm), comm);

  for (size_t i = 0; i < interp_stack->size; ++i)
    yac_interp_stack_config_pack_entry(
      interp_stack->config + i, buffer, buffer_size, position, comm);
}

static void yac_interp_stack_config_unpack_n_string(
  void * buffer, int buffer_size, int * position,
  char * string, int max_string_len, MPI_Comm comm) {

  int string_len;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &string_len, 1, MPI_INT, comm), comm);

  YAC_ASSERT(
    string_len >= 0,
    "ERROR(yac_interp_stack_config_unpack_n_string): invalid string length")

  YAC_ASSERT(
    string_len < max_string_len,
    "ERROR(yac_interp_stack_config_unpack_n_string): string length to long")

  if (string_len > 0)
    yac_mpi_call(
      MPI_Unpack(
        buffer, buffer_size, position, string, string_len, MPI_CHAR, comm),
      comm);
  string[string_len] = '\0';
}

static void yac_interp_stack_config_unpack_entry(
  void * buffer, int buffer_size, int * position,
  union yac_interp_stack_config_entry * entry, MPI_Comm comm) {

  int type;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &type, 1, MPI_INT, comm), comm);

  entry->general.type = (enum yac_interpolation_list)type;

  check_interpolation_type(
    entry->general.type,
    "yac_interp_stack_config_unpack_entry");
  YAC_ASSERT(
    entry->general.type != UNDEFINED,
    "ERROR(yac_interp_stack_config_unpack_entry): "
    "invalid interpolation type")
  switch (type) {
    default:
    case (AVERAGE): {
      int reduction_type;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &reduction_type, 1, MPI_INT, comm),
        comm);
      entry->average.reduction_type =
        (enum yac_interp_avg_weight_type)reduction_type;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &(entry->average.partial_coverage),
          1, MPI_INT, comm), comm);
      break;
    }
    case (RADIAL_BASIS_FUNCTION):
    case (N_NEAREST_NEIGHBOR): {
      int type;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &type, 1, MPI_INT, comm), comm);
      entry->n_nearest_neighbor.config.type =
        (enum yac_interp_nnn_weight_type)type;
      int n;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &n, 1, MPI_INT, comm), comm);
      YAC_ASSERT(
        n >= 0,
        "ERROR(yac_interp_stack_config_unpack_entry): "
        "invalid n_nearest_neighbor.config.n")
      entry->n_nearest_neighbor.config.n = (size_t)n;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->n_nearest_neighbor.config.data.rbf_scale),
          1, MPI_DOUBLE, comm), comm);
      break;
    }
    case (CONSERVATIVE): {
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->conservative.order), 1, MPI_INT, comm), comm);
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->conservative.enforced_conserv), 1, MPI_INT, comm), comm);
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->conservative.partial_coverage), 1, MPI_INT, comm), comm);
      int normalisation;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &normalisation, 1, MPI_INT, comm),
        comm);
      entry->conservative.normalisation =
        (enum yac_interp_method_conserv_normalisation)normalisation;
      break;
    }
    case (SOURCE_TO_TARGET_MAP): {
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->spmap.spread_distance), 1, MPI_DOUBLE, comm), comm);
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->spmap.max_search_distance), 1, MPI_DOUBLE, comm), comm);
      int weight_type;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &weight_type, 1, MPI_INT, comm),
        comm);
      entry->spmap.weight_type =
        (enum yac_interp_spmap_weight_type)weight_type;
      break;
    }
    case (FIXED_VALUE): {
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->fixed.value), 1, MPI_DOUBLE, comm), comm);
      break;
    }
    case (BERNSTEIN_BEZIER):
      break;
    case (USER_FILE): {
      yac_interp_stack_config_unpack_n_string(
        buffer, buffer_size, position,
        entry->user_file.filename, MAX_FILE_NAME_LENGTH, comm);
      yac_interp_stack_config_unpack_n_string(
        buffer, buffer_size, position,
        entry->user_file.src_grid_name, MAX_FILE_NAME_LENGTH, comm);
      yac_interp_stack_config_unpack_n_string(
        buffer, buffer_size, position,
        entry->user_file.tgt_grid_name, MAX_FILE_NAME_LENGTH, comm);
      break;
    }
    case (CHECK): {
      yac_interp_stack_config_unpack_n_string(
        buffer, buffer_size, position,
        entry->check.constructor_key, MAX_ROUTINE_NAME_LENGTH, comm);
      yac_interp_stack_config_unpack_n_string(
        buffer, buffer_size, position,
        entry->check.do_search_key, MAX_ROUTINE_NAME_LENGTH, comm);
      break;
    }
    case (CREEP): {
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->creep.creep_distance), 1, MPI_INT, comm), comm);
      break;
    }
    case (USER_CALLBACK): {
      yac_interp_stack_config_unpack_n_string(
        buffer, buffer_size, position,
        entry->user_callback.func_compute_weights_key,
        MAX_ROUTINE_NAME_LENGTH, comm);
      break;
    }
  }
}

struct yac_interp_stack_config * yac_interp_stack_config_unpack(
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  int stack_size;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &stack_size, 1, MPI_INT, comm), comm);

  YAC_ASSERT(
    stack_size >= 0,
    "ERROR(yac_interp_stack_config_unpack_interp_stack): invalid stack size")

  struct yac_interp_stack_config * interp_stack =
    yac_interp_stack_config_new();

  interp_stack->size = (size_t)stack_size;
  interp_stack->config =
    xmalloc((size_t)stack_size * sizeof(*interp_stack->config));

  for (int i = 0; i < stack_size; ++i)
    yac_interp_stack_config_unpack_entry(
      buffer, buffer_size, position, interp_stack->config + (size_t)i, comm);

  return interp_stack;
}

struct yac_interp_stack_config * yac_interp_stack_config_new() {

  struct yac_interp_stack_config * interp_stack_config =
    xmalloc(1 * sizeof(*interp_stack_config));
  interp_stack_config->config = NULL;
  interp_stack_config->size = 0;

  return interp_stack_config;
}
void yac_interp_stack_config_delete(
  struct yac_interp_stack_config * interp_stack_config) {
  free(interp_stack_config->config);
  free(interp_stack_config);
}

static union yac_interp_stack_config_entry *
  yac_interp_stack_config_add_entry(
    struct yac_interp_stack_config * interp_stack_config) {

  interp_stack_config->size++;
  interp_stack_config->config =
    xrealloc(
      interp_stack_config->config,
      interp_stack_config->size * sizeof(*(interp_stack_config->config)));

  return interp_stack_config->config + (interp_stack_config->size - 1);
}

void yac_interp_stack_config_add_average(
  struct yac_interp_stack_config * interp_stack_config,
  enum yac_interp_avg_weight_type reduction_type, int partial_coverage) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->average.type = AVERAGE;
  entry->average.reduction_type = reduction_type;
  entry->average.partial_coverage = partial_coverage;
}

void yac_interp_stack_config_add_nnn(
  struct yac_interp_stack_config * interp_stack_config,
  enum yac_interp_nnn_weight_type type, size_t n, double scale) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  if (type == NNN_RBF) {
    entry->radial_basis_function.type = RADIAL_BASIS_FUNCTION;
    entry->radial_basis_function.config =
      (struct yac_nnn_config){
        .type = NNN_RBF, .n = n, .data.rbf_scale = scale};
  } else {
    entry->n_nearest_neighbor.type = N_NEAREST_NEIGHBOR;
    entry->n_nearest_neighbor.config =
      (struct yac_nnn_config){
        .type = type, .n = n, .data.rbf_scale = scale};
  }
}

void yac_interp_stack_config_add_conservative(
  struct yac_interp_stack_config * interp_stack_config,
  int order, int enforced_conserv, int partial_coverage,
  enum yac_interp_method_conserv_normalisation normalisation) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->conservative.type = CONSERVATIVE;
  entry->conservative.order = order;
  entry->conservative.enforced_conserv = enforced_conserv;
  entry->conservative.partial_coverage = partial_coverage;
  entry->conservative.normalisation = normalisation;
}

void yac_interp_stack_config_add_spmap(
  struct yac_interp_stack_config * interp_stack_config,
  double spread_distance, double max_search_distance,
  enum yac_interp_spmap_weight_type weight_type) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->spmap.type = SOURCE_TO_TARGET_MAP;
  entry->spmap.spread_distance = spread_distance;
  entry->spmap.max_search_distance = max_search_distance;
  entry->spmap.weight_type = weight_type;
}

void yac_interp_stack_config_add_hcsbb(
  struct yac_interp_stack_config * interp_stack_config) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->hcsbb.type = BERNSTEIN_BEZIER;
}

static void check_string(
  char const * string, char const * file, int line, char const * routine,
  char const * variable) {

  YAC_ASSERT_F(
    string != NULL, "ERROR(%s:%d:%s): %s is NULL",
    file, line, routine, variable)
  YAC_ASSERT_F(
    strlen(string) < MAX_FILE_NAME_LENGTH, "ERROR(%s:%d:%s): %s is too long",
    file, line, routine, variable)
}

void yac_interp_stack_config_add_user_file(
  struct yac_interp_stack_config * interp_stack_config,
  char const * filename, char const * src_grid_name,
  char const * tgt_grid_name) {

  check_string(
    filename, __FILE__, __LINE__,
    "yac_interp_stack_config_add_user_file", "filename");
  check_string(
    src_grid_name, __FILE__, __LINE__,
    "yac_interp_stack_config_add_user_file", "src_grid_name");
  check_string(
    tgt_grid_name, __FILE__, __LINE__,
    "yac_interp_stack_config_add_user_file", "tgt_grid_name");

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->user_file.type = USER_FILE;
  strcpy(entry->user_file.filename, filename);
  strcpy(entry->user_file.src_grid_name, src_grid_name);
  strcpy(entry->user_file.tgt_grid_name, tgt_grid_name);
}

void yac_interp_stack_config_add_fixed(
  struct yac_interp_stack_config * interp_stack_config, double value) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->fixed.type = FIXED_VALUE;
  entry->fixed.value = value;
}

void yac_interp_stack_config_add_check(
  struct yac_interp_stack_config * interp_stack_config,
  char const * constructor_key, char const * do_search_key) {

  YAC_ASSERT_F(
    !constructor_key || strlen(constructor_key) < MAX_ROUTINE_NAME_LENGTH,
    "ERROR(yac_interp_stack_config_add_check): "
    "constructor_key name \"%s\" is too long "
    "(has to be smaller than %d)", constructor_key, MAX_ROUTINE_NAME_LENGTH);
  YAC_ASSERT_F(
    !do_search_key || strlen(do_search_key) < MAX_ROUTINE_NAME_LENGTH,
    "ERROR(yac_interp_stack_config_add_check): "
    "do_search_key name \"%s\" is too long "
    "(has to be smaller than %d)", do_search_key, MAX_ROUTINE_NAME_LENGTH);

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->check.type = CHECK;
  if (constructor_key) {
    strcpy(entry->check.constructor_key, constructor_key);
  } else {
    memset(entry->check.constructor_key, '\0',
           sizeof(entry->check.constructor_key));
  }
  if (do_search_key) {
    strcpy(entry->check.do_search_key, do_search_key);
  } else {
    memset(entry->check.do_search_key, '\0',
           sizeof(entry->check.do_search_key));
  }
}

void yac_interp_stack_config_add_creep(
  struct yac_interp_stack_config * interp_stack_config, int creep_distance) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->creep.type = CREEP;
  entry->creep.creep_distance = creep_distance;
}

void yac_interp_stack_config_add_user_callback(
  struct yac_interp_stack_config * interp_stack_config,
  char const * func_compute_weights_key) {

  check_string(
    func_compute_weights_key, __FILE__, __LINE__,
    "yac_interp_stack_config_add_user_callback",
    "func_compute_weights_key");

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->user_callback.type = USER_CALLBACK;
  strcpy(
    entry->user_callback.func_compute_weights_key, func_compute_weights_key);
}

size_t yac_interp_stack_config_get_size(
  struct yac_interp_stack_config * interp_stack) {

  return interp_stack->size;
}

union yac_interp_stack_config_entry const *
  yac_interp_stack_config_get_entry(
    struct yac_interp_stack_config * interp_stack,
    size_t interp_stack_idx) {

  YAC_ASSERT(
    interp_stack_idx < interp_stack->size,
    "ERROR(yac_interp_stack_config_get_entry): "
    "invalid interpolation stack index");

  return interp_stack->config + interp_stack_idx;
}


enum yac_interpolation_list yac_interp_stack_config_entry_get_type(
  union yac_interp_stack_config_entry const * interp_stack_entry) {

  return interp_stack_entry->general.type;
}

void yac_interp_stack_config_entry_get_average(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  enum yac_interp_avg_weight_type * reduction_type,
  int * partial_coverage) {

  YAC_ASSERT(
    interp_stack_entry->general.type == AVERAGE,
    "ERROR(yac_interp_stack_config_entry_get_average): "
    "wrong interpolation stack entry type");

  *reduction_type = interp_stack_entry->average.reduction_type;
  *partial_coverage = interp_stack_entry->average.partial_coverage;
}

void yac_interp_stack_config_entry_get_nnn(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  enum yac_interp_nnn_weight_type * type, size_t * n, double * scale) {

  YAC_ASSERT(
    (interp_stack_entry->general.type == N_NEAREST_NEIGHBOR) ||
    (interp_stack_entry->general.type == RADIAL_BASIS_FUNCTION),
    "ERROR(yac_interp_stack_config_entry_get_nnn): "
    "wrong interpolation stack entry type");

  *type = interp_stack_entry->n_nearest_neighbor.config.type;
  *n = interp_stack_entry->n_nearest_neighbor.config.n;
  *scale = interp_stack_entry->n_nearest_neighbor.config.data.rbf_scale;
}


void yac_interp_stack_config_entry_get_conservative(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  int * order, int * enforced_conserv, int * partial_coverage,
  enum yac_interp_method_conserv_normalisation * normalisation) {

  YAC_ASSERT(
    interp_stack_entry->general.type == CONSERVATIVE,
    "ERROR(yac_interp_stack_config_entry_get_conservative): "
    "wrong interpolation stack entry type");

  *order = interp_stack_entry->conservative.order;
  *enforced_conserv = interp_stack_entry->conservative.enforced_conserv;
  *partial_coverage = interp_stack_entry->conservative.partial_coverage;
  *normalisation = interp_stack_entry->conservative.normalisation;
}

void yac_interp_stack_config_entry_get_spmap(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  double * spread_distance, double * max_search_distance,
  enum yac_interp_spmap_weight_type * weight_type) {

  YAC_ASSERT(
    interp_stack_entry->general.type == SOURCE_TO_TARGET_MAP,
    "ERROR(yac_interp_stack_config_entry_get_spmap): "
    "wrong interpolation stack entry type");

  *spread_distance = interp_stack_entry->spmap.spread_distance;
  *max_search_distance = interp_stack_entry->spmap.max_search_distance;
  *weight_type = interp_stack_entry->spmap.weight_type;
}

void yac_interp_stack_config_entry_get_user_file(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** filename, char const ** src_grid_name,
  char const ** tgt_grid_name) {

  YAC_ASSERT(
    interp_stack_entry->general.type == USER_FILE,
    "ERROR(yac_interp_stack_config_entry_get_user_file): "
    "wrong interpolation stack entry type");

  *filename = interp_stack_entry->user_file.filename;
  *src_grid_name = interp_stack_entry->user_file.src_grid_name;
  *tgt_grid_name = interp_stack_entry->user_file.tgt_grid_name;
}

void yac_interp_stack_config_entry_get_fixed(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  double * value) {

  YAC_ASSERT(
    interp_stack_entry->general.type == FIXED_VALUE,
    "ERROR(yac_interp_stack_config_entry_get_fixed): "
    "wrong interpolation stack entry type");

  *value = interp_stack_entry->fixed.value;
}

void yac_interp_stack_config_entry_get_check(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** constructor_key, char const ** do_search_key) {

  YAC_ASSERT(
    interp_stack_entry->general.type == CHECK,
    "ERROR(yac_interp_stack_config_entry_get_check): "
    "wrong interpolation stack entry type");

  *constructor_key = interp_stack_entry->check.constructor_key;
  *do_search_key = interp_stack_entry->check.do_search_key;
}

void yac_interp_stack_config_entry_get_creep(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  int * creep_distance) {

  YAC_ASSERT(
    interp_stack_entry->general.type == CREEP,
    "ERROR(yac_interp_stack_config_entry_get_creep): "
    "wrong interpolation stack entry type");

  *creep_distance = interp_stack_entry->creep.creep_distance;
}

void yac_interp_stack_config_entry_get_user_callback(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** func_compute_weights_key) {

  YAC_ASSERT(
    interp_stack_entry->general.type == USER_CALLBACK,
    "ERROR(yac_interp_stack_config_entry_get_user_callback): "
    "wrong interpolation stack entry type");

  *func_compute_weights_key =
    interp_stack_entry->user_callback.func_compute_weights_key;
}
