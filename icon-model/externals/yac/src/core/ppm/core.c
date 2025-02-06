// SPDX-FileCopyrightText: Thomas Jahns <jahns@dkrz.de>
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file core.c
 * @brief interface to user-adjustable core routines of scales ppm
 *
 * @copyright  (C)  2010,2011,2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords: ScalES PPM error handling
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
 * Commentary:
 *
 * The code in this file should be restricted to handle those parts
 * the user program should keep as much control about as possible,
 * like
 *
 * - error handling
 * - file handling
 *
 * Thus the facilities provided here should always come with hooks
 * for user-provided mechanisms.
 *
 * Code:
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>

#include "ppm/core.h"
#include "ppm/symprefix.h"

static MPI_Comm symprefix(default_comm) = MPI_COMM_WORLD;

void
symprefix(set_default_comm)(MPI_Comm comm)
{
  symprefix(default_comm) = comm;
}

// GCOVR_EXCL_START
static void
symprefix(abort_default)(MPI_Comm comm, const char *msg, const char *source, int line)
{
  fprintf(stderr, "Fatal error in %s, line %d: %s\n", source, line, msg);
#if defined (__xlC__) && defined (_AIX)
#pragma omp critical
#endif
  if (SymPrefix(mpi_calls_are_allowed)())
    MPI_Abort(comm, 1);
  abort();
}
// GCOVR_EXCL_STOP

static symprefix(abort_func) symprefix(abort_) =
  (symprefix(abort_func))symprefix(abort_default);

void
symprefix(abort)(MPI_Comm comm, const char *msg, const char *source, int line) {
  symprefix(abort_)(comm, msg, source, line);
}


void symprefix(abort_message)(const char *msg, const char *source, int line)
{
  symprefix(abort_)(symprefix(default_comm), msg, source, line);
}

void
symprefix(restore_default_abort_handler)(void)
{
  symprefix(abort_) = (symprefix(abort_func))symprefix(abort_default);
}

void
symprefix(set_abort_handler)(symprefix(abort_func) custom_abort) {
  symprefix(abort_) = custom_abort;
}

symprefix(abort_func)
symprefix(get_abort_handler)(void) {
  return symprefix(abort_);
}

symprefix(abort_func)
symprefix(get_default_abort_handler)(void) {
  return (symprefix(abort_func))symprefix(abort_default);
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
