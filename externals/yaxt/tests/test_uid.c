/**
 * @file test_uid.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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

#include <mpi.h>
#include <yaxt.h>

#include "tests.h"

int main()
{
  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);
  xt_initialize(MPI_COMM_WORLD);
  Xt_int a_list[] = { 2, 3, 7, 10 };
  enum { a_list_size = sizeof (a_list) / sizeof (a_list[0]) };
  Xt_idxlist a = xt_idxvec_new(a_list, a_list_size);
  Xt_idxlist b = xt_idxlist_copy(a);
  Xt_uid a_id;
  if ((a_id = xt_idxlist_get_uid(a)) == xt_idxlist_get_uid(b))
    PUT_ERR("unexpected duplicate ids");
  enum { section_ndims = 3 };
  Xt_int g_size[section_ndims] = { 5, 5, 5 };
  int l_size[section_ndims] = { 2, 2, 2 };
  Xt_int l_start[section_ndims] = { 1, 2, 0 };
  Xt_idxlist c = xt_idxsection_new(0, section_ndims, g_size, l_size, l_start);
  if (xt_idxlist_get_uid(a) == xt_idxlist_get_uid(c)
      || xt_idxlist_get_uid(b) == xt_idxlist_get_uid(c))
    PUT_ERR("unexpected duplicate ids");
  if (xt_idxlist_get_uid(a) != a_id)
    PUT_ERR("unexpected id change");
  xt_idxlist_delete(c);
  xt_idxlist_delete(b);
  xt_idxlist_delete(a);
  xt_finalize();
  MPI_Finalize();
  return TEST_EXIT_CODE;
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
