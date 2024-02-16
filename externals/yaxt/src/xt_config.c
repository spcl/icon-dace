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
  char name[20];
  Xt_exchanger_new f;
  int code;
} exchanger_table[] = {
  { "irecv_send",
    xt_exchanger_irecv_send_new, xt_exchanger_irecv_send },
  { "irecv_isend",
    xt_exchanger_irecv_isend_new, xt_exchanger_irecv_isend },
  { "irecv_isend_packed",
    xt_exchanger_irecv_isend_packed_new, xt_exchanger_irecv_isend_packed },
  { "mix_irecv_isend",
    xt_exchanger_mix_isend_irecv_new, xt_exchanger_mix_isend_irecv },
  { "neigh_alltoall",
#if MPI_VERSION >= 3
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
#if MPI_VERSION >= 3
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
