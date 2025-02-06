// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <string.h>

#include "geometry.h"
#include "yac_mpi_internal.h"
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
    enum yac_interp_spmap_scale_type scale_type;
    double src_sphere_radius;
    double tgt_sphere_radius;
  } spmap;
  struct {
    enum yac_interpolation_list type;
  } hcsbb;
  struct {
    enum yac_interpolation_list type;
    char filename[YAC_MAX_FILE_NAME_LENGTH];
  } user_file;
  struct {
    enum yac_interpolation_list type;
    double value;
  } fixed;
  struct {
    enum yac_interpolation_list type;
    char constructor_key[YAC_MAX_ROUTINE_NAME_LENGTH];
    char do_search_key[YAC_MAX_ROUTINE_NAME_LENGTH];
  } check;
  struct {
    enum yac_interpolation_list type;
    int creep_distance;
  } creep;
  struct {
    enum yac_interpolation_list type;
    char func_compute_weights_key[YAC_MAX_ROUTINE_NAME_LENGTH];
  } user_callback;
  struct {
    enum yac_interpolation_list type;
    enum yac_interp_ncc_weight_type weight_type;
    int partial_coverage;
  } nearest_corner_cells;
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

static inline void check_interpolation_type(
  enum yac_interpolation_list type, char const * routine) {

  YAC_ASSERT_F(
    (type == YAC_AVERAGE) ||
    (type == YAC_RADIAL_BASIS_FUNCTION) ||
    (type == YAC_N_NEAREST_NEIGHBOR) ||
    (type == YAC_CONSERVATIVE) ||
    (type == YAC_SOURCE_TO_TARGET_MAP) ||
    (type == YAC_FIXED_VALUE) ||
    (type == YAC_BERNSTEIN_BEZIER) ||
    (type == YAC_USER_FILE) ||
    (type == YAC_CHECK) ||
    (type == YAC_CREEP) ||
    (type == YAC_USER_CALLBACK) ||
    (type == YAC_NEAREST_CORNER_CELLS) ||
    (type == YAC_UNDEFINED),
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
  YAC_ASSERT(
    (a->general.type != YAC_UNDEFINED) && (b->general.type != YAC_UNDEFINED),
    "ERROR(yac_interp_stack_config_entry_compare): "
    "interpolation type is undefined");

  if (a->general.type != b->general.type)
    return (a->general.type > b->general.type) -
           (a->general.type < b->general.type);

  switch(a->general.type) {
    default:
    case (YAC_BERNSTEIN_BEZIER):
      return 0;
    case (YAC_AVERAGE): {
      if (a->average.reduction_type != b->average.reduction_type)
        return (a->average.reduction_type > b->average.reduction_type) -
               (a->average.reduction_type < b->average.reduction_type);
      return (a->average.partial_coverage > b->average.partial_coverage) -
             (a->average.partial_coverage < b->average.partial_coverage);
    }
    case (YAC_RADIAL_BASIS_FUNCTION):
    case (YAC_N_NEAREST_NEIGHBOR): {
      if (a->n_nearest_neighbor.config.n !=
          b->n_nearest_neighbor.config.n)
        return (a->n_nearest_neighbor.config.n >
                b->n_nearest_neighbor.config.n) -
               (a->n_nearest_neighbor.config.n <
                b->n_nearest_neighbor.config.n);
      if (fabs(a->n_nearest_neighbor.config.max_search_distance -
               b->n_nearest_neighbor.config.max_search_distance) > yac_angle_tol)
        return (a->n_nearest_neighbor.config.max_search_distance >
                b->n_nearest_neighbor.config.max_search_distance) -
               (a->n_nearest_neighbor.config.max_search_distance <
                b->n_nearest_neighbor.config.max_search_distance);
      if (a->n_nearest_neighbor.config.type !=
          b->n_nearest_neighbor.config.type)
        return (a->n_nearest_neighbor.config.type >
                b->n_nearest_neighbor.config.type) -
               (a->n_nearest_neighbor.config.type <
                b->n_nearest_neighbor.config.type);
      switch (a->n_nearest_neighbor.config.type) {
        case (YAC_INTERP_NNN_GAUSS):
          return (a->n_nearest_neighbor.config.data.gauss_scale >
                  b->n_nearest_neighbor.config.data.gauss_scale) -
                 (a->n_nearest_neighbor.config.data.gauss_scale <
                  b->n_nearest_neighbor.config.data.gauss_scale);
        case (YAC_INTERP_NNN_RBF):
          return (a->n_nearest_neighbor.config.data.rbf_scale >
                  b->n_nearest_neighbor.config.data.rbf_scale) -
                 (a->n_nearest_neighbor.config.data.rbf_scale <
                  b->n_nearest_neighbor.config.data.rbf_scale);
        default:
          return 0;
      };
    }
    case (YAC_CONSERVATIVE): {
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
    case (YAC_SOURCE_TO_TARGET_MAP): {
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
      if (a->spmap.weight_type != b->spmap.weight_type)
        return (a->spmap.weight_type >
                b->spmap.weight_type) -
              (a->spmap.weight_type <
                b->spmap.weight_type);
      if (a->spmap.scale_type != b->spmap.scale_type)
        return (a->spmap.scale_type >
                b->spmap.scale_type) -
               (a->spmap.scale_type <
                b->spmap.scale_type);
      if (a->spmap.src_sphere_radius != b->spmap.src_sphere_radius)
        return (a->spmap.src_sphere_radius >
                b->spmap.src_sphere_radius) -
               (a->spmap.src_sphere_radius <
                b->spmap.src_sphere_radius);
      return (a->spmap.tgt_sphere_radius >
              b->spmap.tgt_sphere_radius) -
             (a->spmap.tgt_sphere_radius <
              b->spmap.tgt_sphere_radius);
    }
    case (YAC_FIXED_VALUE): {
      return (a->fixed.value >
              b->fixed.value) -
             (a->fixed.value <
              b->fixed.value);
    }
    case (YAC_USER_FILE): {
      return strcmp(a->user_file.filename, b->user_file.filename);
    }
    case (YAC_CHECK): {
      int ret;
      if ((ret = strncmp(a->check.constructor_key, b->check.constructor_key,
                         sizeof(a->check.constructor_key)))) return ret;
      return
        strncmp(
          a->check.do_search_key, b->check.do_search_key,
          sizeof(a->check.do_search_key));
    };
    case (YAC_CREEP): {
      return
        (a->creep.creep_distance > b->creep.creep_distance) -
        (a->creep.creep_distance < b->creep.creep_distance);
    }
    case (YAC_USER_CALLBACK): {
      return strcmp(a->user_callback.func_compute_weights_key,
                    b->user_callback.func_compute_weights_key);
    }
    case (YAC_NEAREST_CORNER_CELLS): {
      if (a->nearest_corner_cells.weight_type !=
          b->nearest_corner_cells.weight_type)
        return (a->nearest_corner_cells.weight_type >
                b->nearest_corner_cells.weight_type) -
               (a->nearest_corner_cells.weight_type <
                b->nearest_corner_cells.weight_type);
      return (a->nearest_corner_cells.partial_coverage >
              b->nearest_corner_cells.partial_coverage) -
             (a->nearest_corner_cells.partial_coverage <
              b->nearest_corner_cells.partial_coverage);
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
      interp_stack->config[i].general.type != YAC_UNDEFINED,
      "ERROR(yac_interp_stack_config_generate): "
      "unsupported interpolation method")
    switch((int)(interp_stack->config[i].general.type)) {
      default:
      case(YAC_AVERAGE): {
        enum yac_interp_avg_weight_type weight_type =
          interp_stack->config[i].average.reduction_type;
        int partial_coverage =
          (int)interp_stack->config[i].average.partial_coverage;
        method_stack[i] =
          yac_interp_method_avg_new(weight_type, partial_coverage);
        break;
      }
      case(YAC_CONSERVATIVE): {
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
      case(YAC_FIXED_VALUE): {
        double fixed_value =
          interp_stack->config[i].fixed.value;
        method_stack[i] = yac_interp_method_fixed_new(fixed_value);
        break;
      }
      case(YAC_USER_FILE): {
        char const * weight_file_name =
          interp_stack->config[i].user_file.filename;
        method_stack[i] =
          yac_interp_method_file_new(weight_file_name);
        break;
      }
      case(YAC_CHECK): {
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
      case (YAC_N_NEAREST_NEIGHBOR): {
        method_stack[i] =
          yac_interp_method_nnn_new(
            interp_stack->config[i].
              n_nearest_neighbor.config);
        break;
      }
      case (YAC_BERNSTEIN_BEZIER): {
        method_stack[i] =
          yac_interp_method_hcsbb_new();
        break;
      }
      case (YAC_RADIAL_BASIS_FUNCTION): {
        method_stack[i] =
          yac_interp_method_nnn_new(
            interp_stack->config[i].
              radial_basis_function.config);
        break;
      }
      case (YAC_SOURCE_TO_TARGET_MAP): {
        method_stack[i] =
          yac_interp_method_spmap_new(
            interp_stack->config[i].spmap.spread_distance,
            interp_stack->config[i].spmap.max_search_distance,
            interp_stack->config[i].spmap.weight_type,
            interp_stack->config[i].spmap.scale_type,
            interp_stack->config[i].spmap.src_sphere_radius,
            interp_stack->config[i].spmap.tgt_sphere_radius);
        break;
      }
      case (YAC_CREEP): {
        method_stack[i] =
          yac_interp_method_creep_new(
            interp_stack->config[i].creep.creep_distance);
        break;
      }
      case(YAC_USER_CALLBACK): {
        yac_func_compute_weights compute_weights_callback;
        void * user_data;
        yac_interp_method_callback_get_compute_weights_callback(
          interp_stack->config[i].user_callback.func_compute_weights_key,
          &compute_weights_callback, &user_data);

        method_stack[i] =
          yac_interp_method_callback_new(compute_weights_callback, user_data);
        break;
      }
      case(YAC_NEAREST_CORNER_CELLS): {
        method_stack[i] =
          yac_interp_method_ncc_new(
            interp_stack->config[i].nearest_corner_cells.weight_type,
            interp_stack->config[i].nearest_corner_cells.partial_coverage);
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
    entry->general.type != YAC_UNDEFINED,
    "ERROR(yac_interp_stack_config_get_entry_pack_size): "
    "invalid interpolation type")
  switch (entry->general.type) {
    default:
    case (YAC_AVERAGE):
      return (size_t)int_pack_size + // type
             (size_t)int_pack_size + // reduction_type
             (size_t)int_pack_size;  // partial_coverage
    case (YAC_RADIAL_BASIS_FUNCTION):
    case (YAC_N_NEAREST_NEIGHBOR):
      return (size_t)int_pack_size + // type
             (size_t)int_pack_size + // weight_type
             (size_t)int_pack_size + // n
             (size_t)dbl_pack_size + // max_search_distance
             (size_t)dbl_pack_size;  // scale
    case (YAC_CONSERVATIVE):
      return (size_t)int_pack_size + // type
             (size_t)int_pack_size + // order
             (size_t)int_pack_size + // enforced_conserv
             (size_t)int_pack_size + // partial_coverage
             (size_t)int_pack_size;  // normalisation
    case (YAC_SOURCE_TO_TARGET_MAP):
      return (size_t)int_pack_size + // type
             (size_t)dbl_pack_size + // spread_distance
             (size_t)dbl_pack_size + // max_search_distance
             (size_t)int_pack_size + // weight_type
             (size_t)int_pack_size + // scale_type
             (size_t)dbl_pack_size + // src_sphere_radius
             (size_t)dbl_pack_size;  // tgt_sphere_radius
    case (YAC_FIXED_VALUE):
      return (size_t)int_pack_size + // type
             (size_t)dbl_pack_size;  // value
    case (YAC_BERNSTEIN_BEZIER):
      return (size_t)int_pack_size;  // type
    case (YAC_USER_FILE):
      return (size_t)int_pack_size + // type
             yac_interp_stack_config_get_string_pack_size(
               entry->user_file.filename, comm);
    case (YAC_CHECK):
      return (size_t)int_pack_size + // type
             yac_interp_stack_config_get_string_pack_size(
               entry->check.constructor_key, comm) + // constructor_key
             yac_interp_stack_config_get_string_pack_size(
               entry->check.do_search_key, comm);    // do_search_key
    case (YAC_CREEP):
      return (size_t)int_pack_size + // type
             (size_t)int_pack_size;  // creep_distance
    case (YAC_USER_CALLBACK):
      return (size_t)int_pack_size + // type
             yac_interp_stack_config_get_string_pack_size(
               entry->user_callback.func_compute_weights_key, comm);
               // func_compute_weights_key
    case (YAC_NEAREST_CORNER_CELLS):
      return (size_t)int_pack_size + // type
             (size_t)int_pack_size + // weight_type
             (size_t)int_pack_size;  // partial_coverage
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
    entry->general.type != YAC_UNDEFINED,
    "ERROR(yac_interp_stack_config_pack_entry): "
    "invalid interpolation type")
  switch (entry->general.type) {
    default:
    case (YAC_AVERAGE): {
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
    case (YAC_RADIAL_BASIS_FUNCTION):
    case (YAC_N_NEAREST_NEIGHBOR): {
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
          &(entry->n_nearest_neighbor.config.max_search_distance), 1,
          MPI_DOUBLE, buffer, buffer_size, position, comm), comm);
      yac_mpi_call(
        MPI_Pack(
          &(entry->n_nearest_neighbor.config.data.rbf_scale), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
      break;
    }
    case (YAC_CONSERVATIVE): {
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
    case (YAC_SOURCE_TO_TARGET_MAP): {
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
      int scale_type = (int)(entry->spmap.scale_type);
      yac_mpi_call(
        MPI_Pack(
          &scale_type, 1, MPI_INT, buffer, buffer_size, position, comm),
        comm);
      yac_mpi_call(
        MPI_Pack(
          &(entry->spmap.src_sphere_radius), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
      yac_mpi_call(
        MPI_Pack(
          &(entry->spmap.tgt_sphere_radius), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
      break;
    }
    case (YAC_FIXED_VALUE): {
      yac_mpi_call(
        MPI_Pack(
          &(entry->fixed.value), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
      break;
    }
    case (YAC_BERNSTEIN_BEZIER):
      break;
    case (YAC_USER_FILE): {
      yac_interp_stack_config_pack_string(
        entry->user_file.filename, buffer, buffer_size, position, comm);
      break;
    }
    case (YAC_CHECK): {
      yac_interp_stack_config_pack_string(
        entry->check.constructor_key, buffer, buffer_size, position, comm);
      yac_interp_stack_config_pack_string(
        entry->check.do_search_key, buffer, buffer_size, position, comm);
      break;
    }
    case (YAC_CREEP): {
      yac_mpi_call(
        MPI_Pack(
          &(entry->creep.creep_distance), 1, MPI_INT,
          buffer, buffer_size, position, comm), comm);
      break;
    }
    case (YAC_USER_CALLBACK): {
      yac_interp_stack_config_pack_string(
        entry->user_callback.func_compute_weights_key,
        buffer, buffer_size, position, comm);
      break;
    }
    case (YAC_NEAREST_CORNER_CELLS): {
      int weight_type = (int)(entry->nearest_corner_cells.weight_type);
      yac_mpi_call(
        MPI_Pack(
          &weight_type, 1, MPI_INT, buffer, buffer_size, position, comm),
        comm);
      yac_mpi_call(
        MPI_Pack(
          &(entry->nearest_corner_cells.partial_coverage), 1, MPI_INT,
          buffer, buffer_size, position, comm), comm);
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
    entry->general.type != YAC_UNDEFINED,
    "ERROR(yac_interp_stack_config_unpack_entry): "
    "invalid interpolation type")
  switch (type) {
    default:
    case (YAC_AVERAGE): {
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
    case (YAC_RADIAL_BASIS_FUNCTION):
    case (YAC_N_NEAREST_NEIGHBOR): {
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
          &(entry->n_nearest_neighbor.config.max_search_distance),
          1, MPI_DOUBLE, comm), comm);
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->n_nearest_neighbor.config.data.rbf_scale),
          1, MPI_DOUBLE, comm), comm);
      break;
    }
    case (YAC_CONSERVATIVE): {
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
    case (YAC_SOURCE_TO_TARGET_MAP): {
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
      int scale_type;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &scale_type, 1, MPI_INT, comm),
        comm);
      entry->spmap.scale_type =
        (enum yac_interp_spmap_scale_type)scale_type;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->spmap.src_sphere_radius), 1, MPI_DOUBLE, comm), comm);
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->spmap.tgt_sphere_radius), 1, MPI_DOUBLE, comm), comm);
      break;
    }
    case (YAC_FIXED_VALUE): {
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->fixed.value), 1, MPI_DOUBLE, comm), comm);
      break;
    }
    case (YAC_BERNSTEIN_BEZIER):
      break;
    case (YAC_USER_FILE): {
      yac_interp_stack_config_unpack_n_string(
        buffer, buffer_size, position,
        entry->user_file.filename, YAC_MAX_FILE_NAME_LENGTH, comm);
      break;
    }
    case (YAC_CHECK): {
      yac_interp_stack_config_unpack_n_string(
        buffer, buffer_size, position,
        entry->check.constructor_key, YAC_MAX_ROUTINE_NAME_LENGTH, comm);
      yac_interp_stack_config_unpack_n_string(
        buffer, buffer_size, position,
        entry->check.do_search_key, YAC_MAX_ROUTINE_NAME_LENGTH, comm);
      break;
    }
    case (YAC_CREEP): {
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->creep.creep_distance), 1, MPI_INT, comm), comm);
      break;
    }
    case (YAC_USER_CALLBACK): {
      yac_interp_stack_config_unpack_n_string(
        buffer, buffer_size, position,
        entry->user_callback.func_compute_weights_key,
        YAC_MAX_ROUTINE_NAME_LENGTH, comm);
      break;
    }
    case (YAC_NEAREST_CORNER_CELLS): {
      int weight_type;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &weight_type, 1, MPI_INT, comm),
        comm);
      entry->nearest_corner_cells.weight_type =
        (enum yac_interp_ncc_weight_type)weight_type;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &(entry->nearest_corner_cells.partial_coverage),
          1, MPI_INT, comm), comm);
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

  entry->average.type = YAC_AVERAGE;
  entry->average.reduction_type = reduction_type;
  entry->average.partial_coverage = partial_coverage;
}

void yac_interp_stack_config_add_average_f2c(
  struct yac_interp_stack_config * interp_stack_config,
  int reduction_type, int partial_coverage) {

  YAC_ASSERT(
    (reduction_type == YAC_INTERP_AVG_ARITHMETIC) ||
    (reduction_type == YAC_INTERP_AVG_DIST) ||
    (reduction_type == YAC_INTERP_AVG_BARY),
    "ERROR(yac_interp_stack_config_add_average_f2c): "
    "reduction_type must be one of "
    "YAC_INTERP_AVG_ARITHMETIC/YAC_INTERP_AVG_DIST/YAC_INTERP_AVG_BARY");

  yac_interp_stack_config_add_average(
    interp_stack_config, (enum yac_interp_avg_weight_type)reduction_type,
    partial_coverage);
}

void yac_interp_stack_config_add_ncc(
  struct yac_interp_stack_config * interp_stack_config,
  enum yac_interp_ncc_weight_type weight_type, int partial_coverage) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->nearest_corner_cells.type = YAC_NEAREST_CORNER_CELLS;
  entry->nearest_corner_cells.weight_type = weight_type;
  entry->nearest_corner_cells.partial_coverage = partial_coverage;
}

void yac_interp_stack_config_add_ncc_f2c(
  struct yac_interp_stack_config * interp_stack_config,
  int weight_type, int partial_coverage) {

  YAC_ASSERT(
    (weight_type == YAC_INTERP_NCC_AVG) ||
    (weight_type == YAC_INTERP_NCC_DIST),
    "ERROR(yac_interp_stack_config_add_ncc_f2c): "
    "weight_type must be one of "
    "YAC_INTERP_NCC_AVG/YAC_INTERP_NCC_DIST");

  yac_interp_stack_config_add_ncc(
    interp_stack_config, (enum yac_interp_ncc_weight_type)weight_type,
    partial_coverage);
}

void yac_interp_stack_config_add_nnn(
  struct yac_interp_stack_config * interp_stack_config,
  enum yac_interp_nnn_weight_type type, size_t n,
  double max_search_distance, double scale) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  if (type == YAC_INTERP_NNN_RBF) {
    entry->radial_basis_function.type = YAC_RADIAL_BASIS_FUNCTION;
    entry->radial_basis_function.config =
      (struct yac_nnn_config){
        .type = YAC_INTERP_NNN_RBF, .n = n,
        .max_search_distance = max_search_distance,
        .data.rbf_scale = scale};
  } else {
    entry->n_nearest_neighbor.type = YAC_N_NEAREST_NEIGHBOR;
    entry->n_nearest_neighbor.config =
      (struct yac_nnn_config){
        .type = type, .n = n, .max_search_distance = max_search_distance,
        .data.rbf_scale = scale};
  }
}

void yac_interp_stack_config_add_nnn_f2c(
  struct yac_interp_stack_config * interp_stack_config,
  int type, size_t n, double max_search_distance, double scale) {

  YAC_ASSERT(
    (type == YAC_INTERP_NNN_AVG) ||
    (type == YAC_INTERP_NNN_DIST) ||
    (type == YAC_INTERP_NNN_GAUSS) ||
    (type == YAC_INTERP_NNN_RBF) ||
    (type == YAC_INTERP_NNN_ZERO),
    "ERROR(yac_interp_stack_config_add_nnn_f2c): "
    "type must be one of YAC_INTERP_NNN_AVG/YAC_INTERP_NNN_DIST/"
    "YAC_INTERP_NNN_GAUSS/YAC_INTERP_NNN_RBF/YAC_INTERP_NNN_ZERO.")

  yac_interp_stack_config_add_nnn(
    interp_stack_config, (enum yac_interp_nnn_weight_type)type, n,
    max_search_distance, scale);
}

void yac_interp_stack_config_add_conservative(
  struct yac_interp_stack_config * interp_stack_config,
  int order, int enforced_conserv, int partial_coverage,
  enum yac_interp_method_conserv_normalisation normalisation) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->conservative.type = YAC_CONSERVATIVE;
  entry->conservative.order = order;
  entry->conservative.enforced_conserv = enforced_conserv;
  entry->conservative.partial_coverage = partial_coverage;
  entry->conservative.normalisation = normalisation;
}

void yac_interp_stack_config_add_conservative_f2c(
  struct yac_interp_stack_config * interp_stack_config,
  int order, int enforced_conserv, int partial_coverage,
  int normalisation) {

  YAC_ASSERT(
    (normalisation == YAC_INTERP_CONSERV_DESTAREA) ||
    (normalisation == YAC_INTERP_CONSERV_FRACAREA),
    "ERROR(yac_interp_stack_config_add_conservative_f2c): "
    "type must be one of "
    "YAC_INTERP_CONSERV_DESTAREA/YAC_INTERP_CONSERV_FRACAREA.")

  yac_interp_stack_config_add_conservative(
    interp_stack_config, order, enforced_conserv, partial_coverage,
    (enum yac_interp_method_conserv_normalisation)normalisation);
}

void yac_interp_stack_config_add_spmap(
  struct yac_interp_stack_config * interp_stack_config,
  double spread_distance, double max_search_distance,
  enum yac_interp_spmap_weight_type weight_type,
  enum yac_interp_spmap_scale_type scale_type,
  double src_sphere_radius, double tgt_sphere_radius) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->spmap.type = YAC_SOURCE_TO_TARGET_MAP;
  entry->spmap.spread_distance = spread_distance;
  entry->spmap.max_search_distance = max_search_distance;
  entry->spmap.weight_type = weight_type;
  entry->spmap.scale_type = scale_type;
  entry->spmap.src_sphere_radius = src_sphere_radius;
  entry->spmap.tgt_sphere_radius = tgt_sphere_radius;
}

void yac_interp_stack_config_add_spmap_f2c(
  struct yac_interp_stack_config * interp_stack_config,
  double spread_distance, double max_search_distance,
  int weight_type, int scale_type,
  double src_sphere_radius, double tgt_sphere_radius) {

  YAC_ASSERT(
    (weight_type == YAC_INTERP_SPMAP_AVG) ||
    (weight_type == YAC_INTERP_SPMAP_DIST),
    "ERROR(yac_interp_stack_config_add_spmap_f2c): "
    "weight_type must be one of "
    "YAC_INTERP_SPMAP_AVG/YAC_INTERP_SPMAP_DIST.")

  YAC_ASSERT(
    (scale_type == YAC_INTERP_SPMAP_NONE) ||
    (scale_type == YAC_INTERP_SPMAP_SRCAREA) ||
    (scale_type == YAC_INTERP_SPMAP_INVTGTAREA) ||
    (scale_type == YAC_INTERP_SPMAP_FRACAREA),
    "ERROR(yac_interp_stack_config_add_spmap_f2c): "
    "scale_type must be one of "
    "YAC_INTERP_SPMAP_NONE/YAC_INTERP_SPMAP_SRCAREA/"
    "YAC_INTERP_SPMAP_INVTGTAREA/YAC_INTERP_SPMAP_FRACAREA.")

  yac_interp_stack_config_add_spmap(
    interp_stack_config, spread_distance, max_search_distance,
    (enum yac_interp_spmap_weight_type)weight_type,
    (enum yac_interp_spmap_scale_type)scale_type,
    src_sphere_radius, tgt_sphere_radius);
}

void yac_interp_stack_config_add_hcsbb(
  struct yac_interp_stack_config * interp_stack_config) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->hcsbb.type = YAC_BERNSTEIN_BEZIER;
}

static void check_string(
  char const * string, char const * file, int line, char const * routine,
  char const * variable) {

  YAC_ASSERT_F(
    string != NULL, "ERROR(%s:%d:%s): %s is NULL",
    file, line, routine, variable)
  YAC_ASSERT_F(
    strlen(string) < YAC_MAX_FILE_NAME_LENGTH,
    "ERROR(%s:%d:%s): %s is too long", file, line, routine, variable)
}

void yac_interp_stack_config_add_user_file(
  struct yac_interp_stack_config * interp_stack_config,
  char const * filename) {

  check_string(
    filename, __FILE__, __LINE__,
    "yac_interp_stack_config_add_user_file", "filename");

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->user_file.type = YAC_USER_FILE;
  strcpy(entry->user_file.filename, filename);
}

void yac_interp_stack_config_add_fixed(
  struct yac_interp_stack_config * interp_stack_config, double value) {

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->fixed.type = YAC_FIXED_VALUE;
  entry->fixed.value = value;
}

void yac_interp_stack_config_add_check(
  struct yac_interp_stack_config * interp_stack_config,
  char const * constructor_key, char const * do_search_key) {

  YAC_ASSERT_F(
    !constructor_key || strlen(constructor_key) < YAC_MAX_ROUTINE_NAME_LENGTH,
    "ERROR(yac_interp_stack_config_add_check): "
    "constructor_key name \"%s\" is too long "
    "(has to be smaller than %d)", constructor_key, YAC_MAX_ROUTINE_NAME_LENGTH);
  YAC_ASSERT_F(
    !do_search_key || strlen(do_search_key) < YAC_MAX_ROUTINE_NAME_LENGTH,
    "ERROR(yac_interp_stack_config_add_check): "
    "do_search_key name \"%s\" is too long "
    "(has to be smaller than %d)", do_search_key, YAC_MAX_ROUTINE_NAME_LENGTH);

  union yac_interp_stack_config_entry * entry =
    yac_interp_stack_config_add_entry(interp_stack_config);

  entry->check.type = YAC_CHECK;
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

  entry->creep.type = YAC_CREEP;
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

  entry->user_callback.type = YAC_USER_CALLBACK;
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
    interp_stack_entry->general.type == YAC_AVERAGE,
    "ERROR(yac_interp_stack_config_entry_get_average): "
    "wrong interpolation stack entry type");

  *reduction_type = interp_stack_entry->average.reduction_type;
  *partial_coverage = interp_stack_entry->average.partial_coverage;
}

void yac_interp_stack_config_entry_get_ncc(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  enum yac_interp_ncc_weight_type * type, int * partial_coverage) {

  YAC_ASSERT(
    (interp_stack_entry->general.type == YAC_NEAREST_CORNER_CELLS),
    "ERROR(yac_interp_stack_config_entry_get_ncc): "
    "wrong interpolation stack entry type");

  *type = interp_stack_entry->nearest_corner_cells.weight_type;
  *partial_coverage =
    interp_stack_entry->nearest_corner_cells.partial_coverage;
}

void yac_interp_stack_config_entry_get_nnn(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  enum yac_interp_nnn_weight_type * type, size_t * n,
  double * max_search_distance, double * scale) {

  YAC_ASSERT(
    (interp_stack_entry->general.type == YAC_N_NEAREST_NEIGHBOR) ||
    (interp_stack_entry->general.type == YAC_RADIAL_BASIS_FUNCTION),
    "ERROR(yac_interp_stack_config_entry_get_nnn): "
    "wrong interpolation stack entry type");

  *type = interp_stack_entry->n_nearest_neighbor.config.type;
  *n = interp_stack_entry->n_nearest_neighbor.config.n;
  *max_search_distance =
    interp_stack_entry->n_nearest_neighbor.config.max_search_distance;
  *scale = interp_stack_entry->n_nearest_neighbor.config.data.rbf_scale;
}

void yac_interp_stack_config_entry_get_conservative(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  int * order, int * enforced_conserv, int * partial_coverage,
  enum yac_interp_method_conserv_normalisation * normalisation) {

  YAC_ASSERT(
    interp_stack_entry->general.type == YAC_CONSERVATIVE,
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
  enum yac_interp_spmap_weight_type * weight_type,
  enum yac_interp_spmap_scale_type * scale_type,
  double * src_sphere_radius, double * tgt_sphere_radius) {

  YAC_ASSERT(
    interp_stack_entry->general.type == YAC_SOURCE_TO_TARGET_MAP,
    "ERROR(yac_interp_stack_config_entry_get_spmap): "
    "wrong interpolation stack entry type");

  *spread_distance = interp_stack_entry->spmap.spread_distance;
  *max_search_distance = interp_stack_entry->spmap.max_search_distance;
  *weight_type = interp_stack_entry->spmap.weight_type;
  *scale_type = interp_stack_entry->spmap.scale_type;
  *src_sphere_radius = interp_stack_entry->spmap.src_sphere_radius;
  *tgt_sphere_radius = interp_stack_entry->spmap.tgt_sphere_radius;
}

void yac_interp_stack_config_entry_get_user_file(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** filename) {

  YAC_ASSERT(
    interp_stack_entry->general.type == YAC_USER_FILE,
    "ERROR(yac_interp_stack_config_entry_get_user_file): "
    "wrong interpolation stack entry type");

  *filename = interp_stack_entry->user_file.filename;
}

void yac_interp_stack_config_entry_get_fixed(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  double * value) {

  YAC_ASSERT(
    interp_stack_entry->general.type == YAC_FIXED_VALUE,
    "ERROR(yac_interp_stack_config_entry_get_fixed): "
    "wrong interpolation stack entry type");

  *value = interp_stack_entry->fixed.value;
}

void yac_interp_stack_config_entry_get_check(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** constructor_key, char const ** do_search_key) {

  YAC_ASSERT(
    interp_stack_entry->general.type == YAC_CHECK,
    "ERROR(yac_interp_stack_config_entry_get_check): "
    "wrong interpolation stack entry type");

  *constructor_key = interp_stack_entry->check.constructor_key;
  *do_search_key = interp_stack_entry->check.do_search_key;
}

void yac_interp_stack_config_entry_get_creep(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  int * creep_distance) {

  YAC_ASSERT(
    interp_stack_entry->general.type == YAC_CREEP,
    "ERROR(yac_interp_stack_config_entry_get_creep): "
    "wrong interpolation stack entry type");

  *creep_distance = interp_stack_entry->creep.creep_distance;
}

void yac_interp_stack_config_entry_get_user_callback(
  union yac_interp_stack_config_entry const * interp_stack_entry,
  char const ** func_compute_weights_key) {

  YAC_ASSERT(
    interp_stack_entry->general.type == YAC_USER_CALLBACK,
    "ERROR(yac_interp_stack_config_entry_get_user_callback): "
    "wrong interpolation stack entry type");

  *func_compute_weights_key =
    interp_stack_entry->user_callback.func_compute_weights_key;
}
