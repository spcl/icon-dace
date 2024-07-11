/**
 * @file xt_config.c
 * @brief implementation of configuration object
 *
 * @copyright Copyright  (C)  2020 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 *
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#include <config.h>
#endif

#include <assert.h>
#include <errno.h>
#include <string.h>

#include <mpi.h>

#include <xt/xt_config.h>
#include <xt/xt_mpi.h>
#include "xt_config_internal.h"
#include "xt_exchanger_irecv_send.h"
#include "xt_exchanger_irecv_isend.h"
#include "xt_exchanger_mix_isend_irecv.h"
#include "xt_exchanger_irecv_isend_packed.h"
#include "xt_exchanger_irecv_isend_ddt_packed.h"
#include "xt_exchanger_neigh_alltoall.h"
#include "xt_idxlist_internal.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"

struct Xt_config_ xt_default_config = {
  .exchanger_new = xt_exchanger_mix_isend_irecv_new,
  .exchanger_team_share = NULL,
  .idxv_cnv_size = CHEAP_VECTOR_SIZE,
  .flags = 0,
};

Xt_config xt_config_new(void)
{
  Xt_config config = xmalloc(sizeof(*config));
  *config = xt_default_config;
  return config;
}

void xt_config_delete(Xt_config config)
{
  free(config);
}

static const struct {
  char name[32];
  Xt_exchanger_new f;
  int code;
} exchanger_table[] = {
  { "irecv_send",
    xt_exchanger_irecv_send_new, xt_exchanger_irecv_send },
  { "irecv_isend",
    xt_exchanger_irecv_isend_new, xt_exchanger_irecv_isend },
  { "irecv_isend_packed",
    xt_exchanger_irecv_isend_packed_new, xt_exchanger_irecv_isend_packed },
  { "irecv_isend_ddt_packed",
    xt_exchanger_irecv_isend_ddt_packed_new, xt_exchanger_irecv_isend_ddt_packed },
  { "mix_irecv_isend",
    xt_exchanger_mix_isend_irecv_new, xt_exchanger_mix_isend_irecv },
  { "neigh_alltoall",
#ifdef XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL
    xt_exchanger_neigh_alltoall_new,
#else
    (Xt_exchanger_new)0,
#endif
    xt_exchanger_neigh_alltoall },
};

enum {
  num_exchanger = sizeof (exchanger_table) / sizeof (exchanger_table[0]),
};

int
xt_exchanger_id_by_name(const char *name)
{
  for (size_t i = 0; i < num_exchanger; ++i)
    if (!strcmp(name, exchanger_table[i].name))
      return exchanger_table[i].code;
  return -1;
}

static inline size_t
exchanger_by_function(Xt_exchanger_new exchanger_new)
{
  for (size_t i = 0; i < num_exchanger; ++i)
    if (exchanger_table[i].f == exchanger_new)
      return i;
  return SIZE_MAX;
}


int xt_config_get_exchange_method(Xt_config config)
{
  Xt_exchanger_new exchanger_new = config->exchanger_new;
  size_t eentry = exchanger_by_function(exchanger_new);
  if (eentry != SIZE_MAX)
    return exchanger_table[eentry].code;
  static const char fmt[]
    = "error: unexpected exchanger function (%p)!";
  char buf[sizeof (fmt) + 3*sizeof(void *)];
  sprintf(buf, fmt, (void *)exchanger_new);
  Xt_abort(Xt_default_comm, buf, "xt_config.c", __LINE__);
}

Xt_exchanger_new
xt_config_get_exchange_new_by_comm(Xt_config config, MPI_Comm comm)
{
  Xt_exchanger_new exchanger_new = config->exchanger_new;
#ifdef XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL
  if (exchanger_new == xt_exchanger_neigh_alltoall_new) {
    int flag;
    xt_mpi_call(MPI_Comm_test_inter(comm, &flag), comm);
    if (flag)
      exchanger_new = xt_exchanger_mix_isend_irecv_new;
  }
#else
  (void)comm;
#endif
  return exchanger_new;
}

void xt_config_set_exchange_method(Xt_config config, int method)
{
  static const char fmt[]
    = "error: user-requested exchanger code (%d) does not exist!";
  char buf[sizeof (fmt) + 3*sizeof(int)];
  const char *msg = buf;
  for (size_t i = 0; i < num_exchanger; ++i)
    if (exchanger_table[i].code == method) {
      Xt_exchanger_new exchanger_new;
      if (exchanger_table[i].f) {
        exchanger_new = exchanger_table[i].f;
      } else {
        exchanger_new = xt_default_config.exchanger_new;
        size_t default_entry = exchanger_by_function(exchanger_new);
        if (default_entry == SIZE_MAX) {
          msg = "error: invalid default exchanger constructor!";
          goto abort;
        }
        fprintf(stderr, "warning: %s exchanger unavailable, using "
                "%s instead\n",
                exchanger_table[i].name, exchanger_table[default_entry].name);
      }
      config->exchanger_new = exchanger_new;
      return;
    }
  sprintf(buf, fmt, method);
abort:
  Xt_abort(Xt_default_comm, msg, "xt_config.c", __LINE__);
}

int xt_config_get_idxvec_autoconvert_size(Xt_config config)
{
  return config->idxv_cnv_size;
}

void
xt_config_set_idxvec_autoconvert_size(Xt_config config, int cnvsize)
{
  if (cnvsize > 3)
    config->idxv_cnv_size = cnvsize;
}

int
xt_config_get_redist_mthread_mode(Xt_config config)
{
  return (config->flags & (int32_t)xt_mthread_mode_mask) >> xt_mthread_mode_bit_ofs;
}

void
xt_config_set_redist_mthread_mode(Xt_config config, int mode)
{
  assert(mode >= XT_MT_NONE && mode <= XT_MT_OPENMP);
#ifndef _OPENMP
  if (mode == XT_MT_OPENMP)
    Xt_abort(Xt_default_comm, "error: automatic opening of OpenMP parallel regions requested,"
             " but OpenMP is not configured.\n", "xt_config.c", __LINE__);
#else
  if (mode == XT_MT_OPENMP) {
    int thread_support_provided = MPI_THREAD_SINGLE;
    xt_mpi_call(MPI_Query_thread(&thread_support_provided), Xt_default_comm);
    if (thread_support_provided != MPI_THREAD_MULTIPLE)
      Xt_abort(Xt_default_comm, "error: automatic opening of OpenMP parallel regions requested,\n"
             "        but MPI is not running in thread-safe mode.\n", "xt_config.c", __LINE__);
  }
#endif
  config->flags = (config->flags & ~(int32_t)xt_mthread_mode_mask)
    | ((int32_t)mode << xt_mthread_mode_bit_ofs);
}

void
xt_config_defaults_init(void)
{
  const char *config_env = getenv("XT_CONFIG_DEFAULT_EXCHANGE_METHOD");
  if (config_env) {
    int exchanger_id = xt_exchanger_id_by_name(config_env);
    if (exchanger_id != -1)
      xt_config_set_exchange_method(&xt_default_config, exchanger_id);
    else
      fprintf(stderr, "warning: Unexpected value "
              "for XT_CONFIG_DEFAULT_EXCHANGE_METHOD=%s\n", config_env);
  }
  config_env = getenv("XT_CONFIG_DEFAULT_IDXVEC_AUTOCONVERT_SIZE");
  if (config_env) {
    char *endptr;
    long v = strtol(config_env, &endptr, 0);
    if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN))
        || (errno != 0 && v == 0)) {
      perror("failed to parse value of "
             "XT_CONFIG_DEFAULT_IDXVEC_AUTOCONVERT_SIZE environment variable");
    } else if (endptr == config_env) {
      fputs("malformed value of XT_CONFIG_DEFAULT_IDXVEC_AUTOCONVERT_SIZE"
            " environment variable, no digits were found\n",
            stderr);
    } else if (v < 1 || v > INT_MAX) {
      fprintf(stderr, "value of XT_CONFIG_DEFAULT_IDXVEC_AUTOCONVERT_SIZE"
              " environment variable (%ld) out of range [1,%d]\n",
              v, INT_MAX);
    } else
      xt_config_set_idxvec_autoconvert_size(&xt_default_config, (int)v);
  }
  config_env = getenv("XT_CONFIG_DEFAULT_MULTI_THREAD_MODE");
  if (config_env) {
    char *endptr;
    long v = strtol(config_env, &endptr, 0);
    if (endptr != config_env) {
      if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN))
          || (errno != 0 && v == 0)) {
        perror("failed to parse value of "
               "XT_CONFIG_DEFAULT_MULTI_THREAD_MODE environment variable");
        goto dont_set_mt_mode;
      } else if (v < XT_MT_NONE || v > XT_MT_OPENMP) {
        fprintf(stderr, "numeric value of XT_CONFIG_DEFAULT_MULTI_THREAD_MODE"
                " environment variable (%ld) out of range [0,%d]\n",
                v, XT_MT_OPENMP);
        goto dont_set_mt_mode;
      } else if (*endptr) {
        fprintf(stderr, "trailing text '%s' found after numeric value (%*s) in "
                "XT_CONFIG_DEFAULT_MULTI_THREAD_MODE environment variable\n",
                endptr, (int)(endptr-config_env), config_env);
        goto dont_set_mt_mode;
      }
    } else {
      if (!strcasecmp(config_env, "XT_MT_OPENMP")) {
#ifndef _OPENMP
        fputs("multi-threaded operation requested via "
              "XT_CONFIG_DEFAULT_MULTI_THREAD_MODE, but OpenMP support is not"
              " compiled in!\n", stderr);
        goto dont_set_mt_mode;
#else
        v = XT_MT_OPENMP;
#endif
      } else if (!strcasecmp(config_env, "XT_MT_NONE")) {
        v = XT_MT_NONE;
      } else {
        fputs("unexpected value of XT_CONFIG_DEFAULT_MULTI_THREAD_MODE"
              " environment variable, unrecognized text or numeral\n",
              stderr);
        goto dont_set_mt_mode;
      }
    }
    xt_config_set_redist_mthread_mode(&xt_default_config, (int)v);
  }
dont_set_mt_mode:;
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
