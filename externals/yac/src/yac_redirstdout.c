/**
 * @file yac_redirstdout.c
 *
 * @copyright Copyright  (C)  2014 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @author Rene Redler <rene.redler@mpimet.mpg.de>
 */

/*
 * Keywords:
 * Maintainer: Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
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

/* ---------------------------------------------------------------------
  Date          Programmer   Description
  ------------  ----------   -----------
  01. Dec 2003  R. Redler    created
  01. Dec 2005  H. Ritzdorf  revised
  03. Nov 2021  M. Hanke     refactored

   Copyright 2006-2010, NEC Europe Ltd., London, UK.
   All rights reserved. Use is subject to OASIS4 license terms.
   --------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "yac_interface.h"
#include "utils.h"

/**
 * \example test_redirstdout_c.c
 * This contains an example on how to use \ref yac_redirstdout.
 */

/**
 * \example test_redirstdout.F90
 * This contains an example on how to use \ref yac_redirstdout.
 */

void yac_redirstdout(const char *filestem, int parallel, int my_pe, int npes) {

  YAC_ASSERT(my_pe >= 0, "ERROR(yac_redirstdout): invalid rank")
  YAC_ASSERT(npes >= 1, "ERROR(yac_redirstdout): invalid number of ranks")
  YAC_ASSERT(filestem != NULL, "ERROR(yac_redirstdout): filestem is NULL")

  size_t filestem_size = strlen(filestem);
  YAC_ASSERT(filestem_size >= 1, "ERROR(yac_redirstdout): invalid filestem")

  /* get number of digits in rank */

  register int num = (npes > 0) ? npes : 1;
  register int n_digits = (int) (log10((double)(num) + (double) 0.5)) + 1;
  ASSERT (n_digits > 0)

  /* allocate file names */

  size_t len_alloc = filestem_size + 64;

  char * buffer = (char *) xmalloc (2 * len_alloc);
  char * sname = buffer;
  char * ename = buffer + len_alloc;

  char const * format_string[2][2] =
    {{"%s.log", "%s.err"}, {"%s.%0*d", "%s.err.%0*d"}};

  parallel = parallel == 1;
  snprintf(
    sname, len_alloc-1, format_string[parallel][0], filestem, n_digits, my_pe);
  snprintf(
    ename, len_alloc-1, format_string[parallel][1], filestem, n_digits, my_pe);

  /* redirect stdout/stderr */

  if (parallel || (my_pe == 0)) {
    YAC_ASSERT(
      freopen (sname, "w", stdout),
      "ERROR(yac_redirstdout): failed to redirect stdout")
    YAC_ASSERT(
      freopen (ename, "w", stderr),
      "ERROR(yac_redirstdout): failed to redirect stderr")
  }

  /* free file names */

  free (buffer);
}

