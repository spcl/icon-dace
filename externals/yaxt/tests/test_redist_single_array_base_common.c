/**
 * @file test_redist_single_array_base_common.c
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
#  include <config.h>
#endif

#include "yaxt.h"
#include "tests.h"
#include "test_redist_common.h"

void
test_redist_single_array_base_(int nsend, const struct Xt_redist_msg *send_msgs,
                               int nrecv, const struct Xt_redist_msg *recv_msgs,
                               const void *src_data,
                               size_t num_dst,
                               void *dst_data,
                               prepare_dst dst_prep,
                               const void *dst_prep_info,
                               const void *ref_dst_data,
                               MPI_Datatype dst_data_dt,
                               MPI_Datatype ref_dst_data_dt,
                               MPI_Comm comm,
                               Xt_config config,
                               const char *file, int line)
{
  Xt_redist redist = xt_redist_single_array_base_custom_new(
    nsend, nrecv, send_msgs, recv_msgs, comm, config);
  // test number of send messages
  if (nsend != xt_redist_get_num_send_msg(redist))
    PUT_ERR("error in xt_redist_get_num_send_msg\n");

  // test number of recv messages
  if (nrecv != xt_redist_get_num_recv_msg(redist))
    PUT_ERR("error in xt_redist_get_num_recv_msg\n");

  // test communicator of redist
  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist), comm))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  check_redist_(redist, sync_mode_test_all, 1, &src_data, num_dst, &dst_data,
                dst_data, dst_prep, dst_prep_info, ref_dst_data,
                dst_data_dt, ref_dst_data_dt, file, line);
  {
    Xt_redist redist_copy = xt_redist_copy(redist);
    xt_redist_delete(redist);
    redist = redist_copy;
  }
  check_redist_(redist, sync_mode_test_all, 1, &src_data, num_dst, &dst_data,
                dst_data, dst_prep, dst_prep_info, ref_dst_data,
                dst_data_dt, ref_dst_data_dt, file, line);
  xt_redist_delete(redist);
}

/*
 * Local Variables:
 * coding: utf-8
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
 * license-default: "bsd"
 * End:
 */
