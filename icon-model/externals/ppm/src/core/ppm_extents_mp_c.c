/**
 * @file ppm_extents_mp_c.c --- build MPI datatype for PPM_extent struct
 *
 * @copyright  (C)  2014  Thomas Jahns <jahns@dkrz.de>
 *
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
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
 *
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>

#include <mpi.h>

#include "core/ppm_visibility.h"
#include "core/ppm_extents_mp.h"
#include "core/ppm_xfuncs.h"

MPI_Datatype PPM_extent_mp = MPI_DATATYPE_NULL;

void
PPM_create_extents_mp(void)
{
#pragma omp single
  if (PPM_extent_mp == MPI_DATATYPE_NULL)
  {
    MPI_Datatype elemtype;
#if MPI_VERSION > 2 || ( MPI_VERSION == 2 && MPI_SUBVERSION > 1 )
    elemtype = MPI_INT32_T;
#else
    xmpi(MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof (int32_t), &elemtype));
#endif
    xmpi(MPI_Type_contiguous(2, elemtype, &PPM_extent_mp));
    xmpi(MPI_Type_commit(&PPM_extent_mp));
    MPI_Comm comm_self_clone;
    xmpi(MPI_Comm_dup(MPI_COMM_SELF, &comm_self_clone));
    enum {
      msg_count = 5,
    };
    struct PPM_extent a[msg_count] = { { 123456, 78901 } }, b[msg_count];
    for (size_t i = 1; i < msg_count; ++i)
      a[i] = (struct PPM_extent){ a[0].first + 333 * (int32_t)i,
                                  a[0].size + 555 * (int32_t)i };
    xmpi(MPI_Sendrecv(a, msg_count, PPM_extent_mp, 0, 1,
                      b, msg_count, PPM_extent_mp, 0, 1,
                      comm_self_clone, MPI_STATUS_IGNORE));
    xmpi(MPI_Comm_free(&comm_self_clone));
    bool transfer_worked = true;
    for (size_t i = 0; i < msg_count; ++i)
      transfer_worked &= (a[i].first == b[i].first) & (a[i].size == b[i].size);
    assert(transfer_worked);
  }
}

void
PPM_destroy_extents_mp(void)
{
#pragma omp single
  if (PPM_extent_mp != MPI_DATATYPE_NULL)
    MPI_Type_free(&PPM_extent_mp);
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
