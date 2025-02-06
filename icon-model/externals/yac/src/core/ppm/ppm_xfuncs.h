// SPDX-FileCopyrightText: Thomas Jahns <jahns@dkrz.de>
//
// SPDX-License-Identifier: BSD-3-Clause

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

#ifndef PPM_XFUNCS_INCLUDED
#define PPM_XFUNCS_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include "ppm/symprefix.h"

void *
symprefix(xcalloc)(size_t nmemb, size_t size, const char *source, int line);

void *
symprefix(xmalloc)(size_t size, const char *source, int line);

void *
symprefix(xrealloc)(void *ptr, size_t size, const char *source, int line);

#define xcalloc(nmemb,size)                                     \
  symprefix(xcalloc)((nmemb), (size), __FILE__, __LINE__)
#define xmalloc(size) symprefix(xmalloc)((size), __FILE__, __LINE__)
#define xrealloc(ptr,size)                                      \
  symprefix(xrealloc)((ptr), (size), __FILE__, __LINE__)

FILE *
symprefix(xfopen)(const char *path, const char *mode,
                  const char *source, int line);

#define xfopen(path, mode) symprefix(xfopen)(path, mode, __FILE__, __LINE__)

void
symprefix(xfclose)(FILE *fp, const char *source, int line);

#define xfclose(fp) symprefix(xfclose)(fp, __FILE__, __LINE__)

#endif
/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
