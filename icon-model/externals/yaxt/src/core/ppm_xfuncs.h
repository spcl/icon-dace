/**
 * @file ppm_xfuncs.h
 * @brief add versions of standard API functions not returning on error
 *
 * @copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords: fail-safe wrappers
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
 */

#ifndef PPM_XFUNCS_INCLUDED
#define PPM_XFUNCS_INCLUDED

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#ifdef USE_MPI
#  include <mpi.h>
#endif

#include "core/symprefix.h"

void *
SymPrefix(xcalloc)(size_t nmemb, size_t size, const char *source, int line);

void *
SymPrefix(xmalloc)(size_t size, const char *source, int line);

void *
SymPrefix(xrealloc)(void *ptr, size_t size, const char *source, int line);

#define xcalloc(nmemb,size)                                     \
  SymPrefix(xcalloc)((nmemb), (size), __FILE__, __LINE__)
#define xmalloc(size) SymPrefix(xmalloc)((size), __FILE__, __LINE__)
#define xrealloc(ptr,size)                                      \
  SymPrefix(xrealloc)((ptr), (size), __FILE__, __LINE__)

FILE *
SymPrefix(xfopen)(const char *path, const char *mode,
                  const char *source, int line);

#define xfopen(path, mode) SymPrefix(xfopen)(path, mode, __FILE__, __LINE__)

void
SymPrefix(xfclose)(FILE *fp, const char *source, int line);

#define xfclose(fp) SymPrefix(xfclose)(fp, __FILE__, __LINE__)

int
SymPrefix(xfputc)(int c, FILE *stream, const char *source, int line);

#define xfputc(c,stream) SymPrefix(xfputc)((c),(stream), __FILE__, __LINE__)

#ifdef USE_MPI
void
SymPrefix(xmpi)(int errcode, const char *source, int line);

#define xmpi(errcode)                           \
  do {                                          \
    if (errcode == MPI_SUCCESS) ; else          \
      SymPrefix(xmpi)(errcode, __FILE__, __LINE__);     \
  } while (0)
#endif

#endif
/*
 * Local Variables:
 * license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
