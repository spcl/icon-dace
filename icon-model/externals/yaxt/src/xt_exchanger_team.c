/**
 * @file xt_exchanger_team.c
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

#include "xt_exchanger.h"



size_t
xt_exchanger_team_get_share_size(Xt_exchanger exchanger)
{
  return exchanger->vtable->team_share_size;
}

size_t
xt_exchanger_new_team_get_share_size(Xt_exchanger_new exchanger_new)
{
  return xt_exchanger_new_get_vtable(exchanger_new)->team_share_size;
}

void
xt_exchanger_team_share_default_init(Xt_exchanger exchanger,
                                     void *share)
{
  const struct xt_exchanger_vtable *vtab = exchanger->vtable;
  if (vtab && vtab->team_share_default_init)
    vtab->team_share_default_init(share);
}

void
xt_exchanger_new_team_share_default_init(Xt_exchanger_new exchanger_new,
                                         void *share)
{
  const struct xt_exchanger_vtable *vtab
    = xt_exchanger_new_get_vtable(exchanger_new);
  if (vtab && vtab->team_share_default_init)
    vtab->team_share_default_init(share);
}

void
xt_exchanger_team_share_destroy(Xt_exchanger exchanger, void *share)
{
  const struct xt_exchanger_vtable *vtab = exchanger->vtable;
  if (vtab && vtab->team_share_destroy)
    vtab->team_share_destroy(share);
}

void
xt_exchanger_new_team_share_destroy(Xt_exchanger_new exchanger_new,
                                    void *share)
{
  const struct xt_exchanger_vtable *vtab
    = xt_exchanger_new_get_vtable(exchanger_new);
  if (vtab && vtab->team_share_destroy)
    vtab->team_share_destroy(share);
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

