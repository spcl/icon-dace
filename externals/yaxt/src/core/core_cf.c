/*
 * @file core_cf.c
 * @brief ScalES-PPM core library C/Fortran interface
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * Keywords:
 * @author Thomas Jahns <jahns@dkrz.de>
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
#  include <config.h>
#endif
#include <stdio.h>

#define FCALLSC_QUALIFIER PPM_DSO_INTERNAL

#if defined __clang__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wreserved-id-macro"
#endif
#define CF_USE_ALLOCA 1
#include <cfortran.h>
#if defined __clang__
#  pragma GCC diagnostic pop
#endif


#include "core/ppm_visibility.h"
#include "core.h"

static void
SymPrefix(set_default_comm_f)(MPI_Fint *comm_f)
{
#if defined(USE_MPI)
  int flag = 0;
  MPI_Comm comm_c;
#if defined (__xlC__) && defined (_AIX)
#pragma omp critical
#endif
  comm_c = (MPI_Initialized(&flag) == MPI_SUCCESS && flag)?
    MPI_Comm_f2c(*comm_f):SymPrefix(default_comm);
#else
  MPI_Comm comm_c = *comm_f;
#endif
  SymPrefix(default_comm) = comm_c;
}

FCALLSCSUB1(SymPrefix(set_default_comm_f), SYMPREFIX(SET_DEFAULT_COMM),
            symprefix(set_default_comm), PVOID)

static void
abort_f(MPI_Fint *comm_f, const char *msg,
                   const char *source, int line)
  __attribute__((noreturn));

static void
abort_f(MPI_Fint *comm_f, const char *msg,
                   const char *source, int line)
{
  MPI_Comm comm_c = MPI_COMM_NULL;
#if defined(USE_MPI)
  int flag = 0;
#if defined (__xlC__) && defined (_AIX)
#pragma omp critical
#endif
  if (MPI_Initialized(&flag) == MPI_SUCCESS && flag)
    comm_c = MPI_Comm_f2c(*comm_f);
#else
  comm_c = *comm_f;
#endif
  SymPrefix(abort)(comm_c, msg, source, line);
}

#undef CFattributes
#define CFattributes __attribute__((noreturn))
FCALLSCSUB4(abort_f, SYMPREFIX(ABORT), symprefix(abort),
            PVOID, STRING, STRING, INT)
#undef CFattributes
#define CFattributes

static void
SymPrefix(abort_handler_wrapper)(MPI_Comm comm, const char msg[],
                          const char source[], int line)
  __attribute__((noreturn));

static void
SymPrefix(set_abort_handler_f)(void (*abort_handler)(void));

#undef ROUTINE_1
#define ROUTINE_1 (void (*)(void))
FCALLSCSUB1(SymPrefix(set_abort_handler_f), SYMPREFIX(SET_ABORT_HANDLER),
            symprefix(set_abort_handler), ROUTINE)

static void
abort_default_f(MPI_Fint *comm_f, const char *msg, const char *source,
                int line)
  __attribute__((noreturn));

static void
abort_default_f(MPI_Fint *comm_f, const char *msg, const char *source,
                int line)
{
#if defined(USE_MPI)
  int flag = 0;
  MPI_Comm comm_c = (MPI_Initialized(&flag) == MPI_SUCCESS && flag)?
    MPI_Comm_f2c(*comm_f):SymPrefix(default_comm);
#else
  MPI_Comm comm_c = *comm_f;
#endif
  SymPrefix(abort_default)(comm_c, msg, source, line);
}

#undef FCALLSC_QUALIFIER
#define FCALLSC_QUALIFIER

FCALLSCSUB0(SymPrefix(restore_default_abort_handler),
            SYMPREFIX(RESTORE_DEFAULT_ABORT_HNDL),
            symprefix(restore_default_abort_hndl))

#undef CFattributes
#define CFattributes __attribute__((noreturn))
FCALLSCSUB4(abort_default_f,SYMPREFIX(ABORT_DEFAULT),
            symprefix(abort_default),
            PVOID,STRING,STRING,INT)
#undef CFattributes
#define CFattributes

/* this must be the last piece of code in the file because we
 * redefine a cfortran.h internal here, to allow calls to Fortran
 * function pointers */
#undef CFC_
#define CFC_(UN,LN) (UN)
#undef CFextern
#define CFextern typedef
__attribute__((noreturn))
PROTOCCALLSFSUB4(*SymPrefix(fortran_abort_func),,PVOID,STRING,STRING,INT)
#undef CFextern

static SymPrefix(fortran_abort_func) SymPrefix(fortran_abort_fp);

static void
SymPrefix(set_abort_handler_f)(void (*abort_handler)(void))
{
  SymPrefix(fortran_abort_fp)
    = (SymPrefix(fortran_abort_func))abort_handler;
  SymPrefix(abort) = SymPrefix(abort_handler_wrapper);
}

static void
SymPrefix(abort_handler_wrapper)(MPI_Comm comm, const char msg[],
                                 const char source[], int line)
{
#if defined(USE_MPI)
  int flag = 0;
  MPI_Fint comm_f;
#if defined (__xlC__) && defined (_AIX)
#pragma omp critical
#endif
  comm_f = (MPI_Initialized(&flag) == MPI_SUCCESS && flag)?
    MPI_Comm_c2f(comm):(MPI_Fint)0;
#else
  MPI_Fint comm_f = comm;
#endif
  /* cfortran.h does not understand const char * */
  char *msg_arg = (char *)msg, *source_arg = (char *)source;
#undef CPPPROTOCLSFSUB14
#define CPPPROTOCLSFSUB14(UN,LN,T1,T2,T3,T4,T5,T6,T7,T8,T9,TA,TB,TC,TD,TE)
  CCALLSFSUB4(*SymPrefix(fortran_abort_fp),,PVOID,STRING,STRING,INT,
              &comm_f, msg_arg, source_arg, line);
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
