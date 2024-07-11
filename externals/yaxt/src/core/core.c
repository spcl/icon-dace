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
 * URL: https://swprojects.dkrz.de/redmine/projects/scales-ppm
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
#  include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#ifndef WITHOUT_FORTRAN
#  if defined __clang__
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wreserved-id-macro"
#  endif
#  include <cfortran.h>
#  if defined __clang__
#    pragma GCC diagnostic pop
#  endif
#  if defined __GNUC__
#    undef COMMON_BLOCK_DEF
#    define COMMON_BLOCK_DEF(DEFINITION, NAME) extern DEFINITION NAME __attribute__((aligned(16),common))
#  endif
#endif
#ifdef USE_MPI
#  include <mpi.h>
#endif

#include "core/ppm_visibility.h"
#include "core/core.h"
#include "core/symprefix.h"

MPI_Comm SymPrefix(default_comm) = MPI_COMM_WORLD;

#ifndef WITHOUT_FORTRAN
#define F2C_Data COMMON_BLOCK(SYMPREFIX(F2C_DATA),symprefix(f2c_data))

typedef struct
{
  MPI_Fint symprefix(default_comm);
} SymPrefix(F2C_Def);

COMMON_BLOCK_DEF(SymPrefix(F2C_Def),F2C_Data);

SymPrefix(F2C_Def) F2C_Data
#  if defined __GNUC__
                   __attribute__((aligned(16),common))
#  endif
#endif
                   ;
void
SymPrefix(set_default_comm)(MPI_Comm comm)
{
#ifndef WITHOUT_FORTRAN
  MPI_Fint comm_f;
#  if defined(USE_MPI)
  comm_f = MPI_Comm_c2f(comm);
#  else
  comm_f = comm;
#  endif
  F2C_Data.symprefix(default_comm) = comm_f;
#endif
  SymPrefix(default_comm) = comm;
}

void
SymPrefix(abort_default)(MPI_Comm comm, const char *msg, const char *source, int line)
{
  fprintf(stderr, "Fatal error in %s, line %d: %s\n", source, line, msg);
#ifdef USE_MPI
#if defined (__xlC__) && defined (_AIX)
#pragma omp critical
#endif
  if (SymPrefix(mpi_calls_are_allowed)())
    MPI_Abort(comm, 1);
  else
    abort();
#else
  (void)comm;
#endif
  abort();
}

SymPrefix(abort_func) SymPrefix(abort) = SymPrefix(abort_default);

void
SymPrefix(restore_default_abort_handler)(void)
{
  SymPrefix(abort) = SymPrefix(abort_default);
}

/*
 * Local Variables:
 * license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
