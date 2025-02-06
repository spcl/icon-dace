/**
 * @file ppm_c.c
 * @brief One-time initialization of PPM functionality
 *
 * @copyright Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
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

#include <stdbool.h>

#include "ppm.h"

#include "core/core.h"
#if USE_MPI
#include "core/ppm_extents_mp.h"
#include "core/ppm_std_type_kinds_mp.h"
#endif
#include "core/ppm_random.h"
#include "core/ppm_math_extensions.h"

static bool PPM_initialized_state = false,
  PPM_finalized_state = false;

void
PPM_initialize(MPI_Comm *default_comm, int *random_seed, unsigned *seed_output)
{
  int seed = random_seed ? *random_seed : 0;
  if (default_comm)
    PPM_set_default_comm(*default_comm);
  unsigned seed_out = PPM_ya_rand_init(PPM_default_comm, seed);
  if (seed_output)
    *seed_output = seed_out;
#ifdef USE_MPI
  PPM_create_extents_mp();
  PPM_initialize_std_type_kinds_mp();
  PPM_initialize_math_extensions_mp();
#endif
#pragma omp barrier
#pragma omp master
  PPM_initialized_state = true;
#pragma omp barrier
}

void
PPM_finalize(void)
{
#ifdef USE_MPI
  PPM_finalize_math_extensions_mp();
  PPM_finalize_std_type_kinds_mp();
  PPM_destroy_extents_mp();
#endif
  PPM_ya_rand_finish();
#pragma omp barrier
#pragma omp master
  PPM_finalized_state = true;
#pragma omp barrier
}

bool
PPM_initialized(void)
{
  bool state;
#pragma omp critical
  state = PPM_initialized_state;
  return state;
}


bool
PPM_finalized(void)
{
  bool state;
#pragma omp critical
  state = PPM_finalized_state;
  return state;
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
