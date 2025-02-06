/**
 * @file xt_config_internal.h
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

#include <stdint.h>

#include <mpi.h>

#include "core/ppm_visibility.h"
#include "xt/xt_redist.h"
#include "xt/xt_config.h"
#include "xt_exchanger.h"

enum xt_config_flags {
  exch_no_dt_dup = (1 << 0),
  xt_mthread_mode_bit_ofs = 1,
  xt_mthread_mode_mask = (1 << xt_mthread_mode_bit_ofs),
};


struct Xt_config_ {
  /**
   * constructor to use when creating the exchanger of a redist
   */
  Xt_exchanger_new exchanger_new;
  /**
   * pointer to exchanger team share data
   */
  void *exchanger_team_share;
  /**
   * automatically compress index lists of vector type at this size
   * into another representation to save on computation/memory overall
   */
  int idxv_cnv_size;
  /**
   * binary combination of xt_config_flags */
  int32_t flags;
};

extern struct Xt_config_ xt_default_config;

PPM_DSO_INTERNAL void
xt_config_defaults_init(void);

/**
 * Get appropriate exchanger constructor.
 *
 * @param config configuration object
 * @param comm communicator to use the constructor with
 * @returns configured exchanger constructor, or a fallback if the
 * configured constructor does not apply to \a comm.
 */
PPM_DSO_INTERNAL Xt_exchanger_new
xt_config_get_exchange_new_by_comm(Xt_config config, MPI_Comm comm);

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
