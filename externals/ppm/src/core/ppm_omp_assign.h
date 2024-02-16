/**
 * @file ppm_omp_assign.h
 * @brief multi-threaded array-filling core routine body include
 *
 * Since OpenMP 2.5 does not allow for unsigned loop counters and a
 * macro could not capture the OpenMP pragmas, this repetitive code
 * part is factored out into this file.
 *
 *
 * Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords: PRNG OpenMP
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

#ifndef RHS
#error "right hand side of assignment to a[i] must be defined as macro RHS"
#endif

/* the following code does assignment of RHS to a[i] for i in [0,n)
 * and typeof(n) == size_t */

#if defined(_OPENMP) && _OPENMP <= 200505
  ssize_t i;
  if (n > SSIZE_MAX)
  {
    ssize_t n2 = SSIZE_MAX;
#pragma omp for
    for (i = 0; i < n2; ++i)
      a[i] = RHS;
    n -= (size_t)SSIZE_MAX;
    a += SSIZE_MAX;
  }
#pragma omp for
  for (i = 0; i < (ssize_t)n; ++i)
    a[i] = RHS;
#else
  /* OpenMP 3.0 unsigned loop counters supported */
  size_t i;
#pragma omp for
  for (i = 0; i < n; ++i)
    a[i] = RHS;
#endif
/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
