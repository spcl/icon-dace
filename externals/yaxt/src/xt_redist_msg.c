/**
 * @file xt_redist_msg.c
 *
 * @copyright Copyright  (C)  2021 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
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

#include "xt/xt_mpi.h"
#include "xt_redist_internal.h"

void
xt_redist_msgs_strided_copy(size_t n,
                            const struct Xt_redist_msg *restrict src,
                            size_t src_stride,
                            struct Xt_redist_msg *restrict dst,
                            size_t dst_stride,
                            MPI_Comm comm, bool dt_dup) {


  const unsigned char *restrict src_store = (const unsigned char *)src;
  unsigned char *restrict dst_store = (unsigned char *)dst;
  for (size_t i = 0; i < n; ++i) {
    const struct Xt_redist_msg *restrict src_msg
      = (const struct Xt_redist_msg *)(const void *)(src_store + i * src_stride);
    struct Xt_redist_msg *dst_msg
      = (struct Xt_redist_msg *)(void *)(dst_store + i * dst_stride);
    dst_msg->rank = src_msg->rank;
    if (dt_dup)
      xt_mpi_call(MPI_Type_dup(src_msg->datatype, &(dst_msg->datatype)), comm);
    else
      dst_msg->datatype = src_msg->datatype;
  }
}

void xt_redist_msgs_strided_destruct(size_t n, struct Xt_redist_msg *msgs,
                                     MPI_Comm comm, size_t ofs_stride) {

  unsigned char *restrict msgs_store = (unsigned char *)msgs;
  for (size_t i = 0; i < n; ++i) {
    MPI_Datatype *dt
      = &(((struct Xt_redist_msg *)(void *)(msgs_store + i * ofs_stride))->datatype);
    if (*dt != MPI_DATATYPE_NULL)
      xt_mpi_call(MPI_Type_free(dt), comm);
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
