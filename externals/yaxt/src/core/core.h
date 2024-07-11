/**
 * @file core.h --- interface to user-adjustable core routines of scales ppm
 *
 * @copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
 *
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
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
#ifndef PPM_CORE_H
#define PPM_CORE_H
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#else
/**
 * fall back to int in case no MPI implementation was found
 */
typedef int MPI_Comm;
/**
 * fall back to int in case no MPI implementation was found
 */
typedef int MPI_Fint;
/**
 * provide value to use in case no actual communicators exist
 */
enum {
  MPI_COMM_WORLD = 4711,
  MPI_COMM_NULL = 0,
};
#endif

#include "core/symprefix.h"

/* If we're not using GNU C, elide __attribute__ */
#ifndef __GNUC__
#  define  __attribute__(x)  /*NOTHING*/
#endif
#define XT_UNUSED(x) UNUSED_ ## x __attribute__((__unused__))


/**
 * functions used as error handler must conform to this interface
 */
typedef void (*SymPrefix(abort_func))(MPI_Comm comm, const char *msg,
                                      const char *source, int line)
  __attribute__((noreturn));
/**
 * Unless modified, this function pointer will reference PPM_abort_default.
 */
extern SymPrefix(abort_func) SymPrefix(abort);

/**
 * Restore default abort handler.
 */
void
SymPrefix(restore_default_abort_handler)(void);

/**
 * communicator object to use by default
 */
extern MPI_Comm SymPrefix(default_comm);

/**
 * This function
 * prints the message argument and file and line of the error
 * to standard error, and
 * calls either MPI_Abort or abort depending on whether
 * MPI is initialized.
 * @param comm MPI communcator object to use on call to MPI_Abort
 * @param msg message text to print
 * @param source string describing source file name
 * @param line line number of caller
 */
extern void
SymPrefix(abort_default)(MPI_Comm comm, const char *msg,
                         const char *source, int line)
  __attribute__((noreturn));

/**
 * change default communicator object
 */
extern void
SymPrefix(set_default_comm)(MPI_Comm comm);

#define die(msg) \
  SymPrefix(abort)(SymPrefix(default_comm), (msg), __FILE__, __LINE__)

#ifdef USE_MPI
static inline int
SymPrefix(mpi_calls_are_allowed)(void)
{
  int init_flag = 0, finished_flag = 0;
  return MPI_Initialized(&init_flag) == MPI_SUCCESS && init_flag
    && MPI_Finalized(&finished_flag) == MPI_SUCCESS && !finished_flag;
}
#endif

#endif
/*
 * Local Variables:
 * license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * license-markup: "doxygen"
 * End:
 */
