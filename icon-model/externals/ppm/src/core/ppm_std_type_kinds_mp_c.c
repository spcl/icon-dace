/**
 * @file ppm_std_type_kinds_mp_c.c
 * @brief implementation constants and declarations for basic types to
 * be used in MPI programs, C implementation
 *
 * @copyright Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
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
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "core/ppm_std_type_kinds_mp.h"
#include "core/ppm_xfuncs.h"
#include "core/ppm_visibility.h"


#if !(MPI_VERSION > 2 || (MPI_VERSION == 2 && MPI_SUBVERSION > 1))
MPI_Datatype PPM_DT_C_DOUBLE_PAIR_MP PPM_DSO_API_EXPORT = MPI_DATATYPE_NULL;
#endif

void
PPM_initialize_std_type_kinds_mp(void)
{
#if !(MPI_VERSION > 2 || (MPI_VERSION == 2 && MPI_SUBVERSION > 1))
  xmpi(MPI_Type_contiguous(2, MPI_DOUBLE, &PPM_DT_C_DOUBLE_PAIR_MP));
  xmpi(MPI_Type_commit(&PPM_DT_C_DOUBLE_PAIR_MP));
#endif
}


void
PPM_finalize_std_type_kinds_mp(void)
{
#if !(MPI_VERSION > 2 || (MPI_VERSION == 2 && MPI_SUBVERSION > 1))
  xmpi(MPI_Type_free(&PPM_DT_C_DOUBLE_PAIR_MP));
#endif
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
