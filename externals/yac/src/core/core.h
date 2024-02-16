/**
 * @file core.h --- interface to user-adjustable core routines of scales ppm
 *
 * @copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
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

#include <mpi.h>
#include <yaxt.h>

#include "core/symprefix.h"

/* If we're not using GNU C, elide __attribute__ */
#ifndef __GNUC__
#  define  __attribute__(x)  /*NOTHING*/
#endif
// #define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#define UNUSED(x) (void)(x)

// idxtype from yaxt as global id type in yac
typedef Xt_int yac_int;
#define yac_int_dt Xt_int_dt

/**
 * functions used as error handler must conform to this interface
 */
typedef void (*symprefix(abort_func))(MPI_Comm comm, const char *msg,
                                      const char *source, int line)
  __attribute__((noreturn));

/**
 * Calls the currently set abort handler (\ref abort_default "symprefix(abort_default)" by default)
 * @param[in] comm   MPI communicator used to call MPI_Abort
 * @param[in] msg    message text to print
 * @param[in] source string describing source file name
 * @param[in] line   line number of caller
 */
void symprefix(abort)(MPI_Comm comm, const char *msg,
                      const char *source, int line)
  __attribute__((noreturn));

/**
 * Call the \ref abort "symprefix(abort)" function (providing the default communicator for
 * the comm argument).
 * @param msg    message text to print
 * @param source string describing source file name
 * @param line   line number of caller
 */
void symprefix(abort_message)(char const *msg, const char *source, int line);

/**
 * Restore default abort handler.
 */
void
symprefix(restore_default_abort_handler)(void);

/**
 * Set custom abort handler.
 * @param[in] custom_abort custom abort handler
 */
void
symprefix(set_abort_handler)(symprefix(abort_func) custom_abort);

/**
 * Get abort handler.
 * @return currently set abort handler
 */
symprefix(abort_func) symprefix(get_abort_handler)(void);

/**
 * Get default abort handler.
 * @return currently set abort handler
 */
symprefix(abort_func) symprefix(get_default_abort_handler)(void);

/**
 * change default communicator object
 */
extern void
symprefix(set_default_comm)(MPI_Comm comm);

#define die(msg) \
  symprefix(abort_message)((msg), __FILE__, __LINE__)

static inline int
SymPrefix(mpi_calls_are_allowed)(void)
{
  int init_flag = 0, finished_flag = 0;
  return MPI_Initialized(&init_flag) == MPI_SUCCESS && init_flag
    && MPI_Finalized(&finished_flag) == MPI_SUCCESS && !finished_flag;
}

#endif
/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * license-markup: "doxygen"
 * End:
 */
