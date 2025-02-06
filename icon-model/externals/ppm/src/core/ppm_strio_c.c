/**
 * @file ppm_strio_c.c
 * @brief flexible string I/O for Fortran
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
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
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _POSIX_C_SOURCE 200112L

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "cfortran.h"
#include "core/core.h"
#include "core/minmax.h"

struct PPM_fmt_elem
{
  int32_t argtype, flags;
  void *addr;
};

enum {
  arg_i4 = 1,
  arg_i8 = 2,
  arg_sp = 3,
  arg_dp = 4,
};

enum {
  no_error,
  parse_error = 1,
};


static ssize_t
PPM_sscana(const char *str, struct PPM_fmt_elem out[],
           size_t out_size, int *ierror)
{
  size_t i;
  char *next;
  assert(str);
  for (i = 0; i < out_size && *str; ++i)
  {
    switch(out[i].argtype) {
    case arg_i4:
      {
        long val;
        errno = 0;
        val = strtol(str, &next, 0);
        if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
            || (errno != 0 && val == 0) || (next == str)
            || val > INT32_MAX || val < INT32_MIN) {
          if (errno)
            perror("PPM_sscana");
          else
            fputs("PPM_sscana: parse error\n", stderr);
          *ierror = parse_error;
          return (ssize_t)MIN(i, (size_t)SSIZE_MAX);
        }
        *(int32_t *)(out[i].addr) = (int32_t)val;
        str = next;
      }
      break;
    case arg_i8:
      {
        long long val;
        errno = 0;
        val = strtoll(str, &next, 0);
        if ((errno == ERANGE && (val == LLONG_MAX || val == LLONG_MIN))
            || (errno != 0 && val == 0) || (next == str)
            || val > INT64_MAX || val < INT64_MIN) {
          if (errno)
            perror("PPM_sscana");
          else
            fputs("PPM_sscana: parse error\n", stderr);
          *ierror = parse_error;
          return (ssize_t) MIN(i, (size_t)SSIZE_MAX);
        }
        *(int64_t *)(out[i].addr) = (int64_t)val;
        str = next;
      }
      break;
    case arg_sp:
      {
        float val;
        errno = 0;
        val = strtof(str, &next);
        if ((errno == ERANGE && (val == HUGE_VALF || val == -HUGE_VALF))
            || (errno != 0 && val == 0.0f) || (next == str)) {
          if (errno)
            perror("PPM_sscana");
          else
            fputs("PPM_sscana: parse error\n", stderr);
          *ierror = parse_error;
          return (ssize_t) MIN(i, (size_t)SSIZE_MAX);
        }
        *(float *)(out[i].addr) = val;
        str = next;
      }
      break;
    case arg_dp:
      {
        double val;
        errno = 0;
        val = strtod(str, &next);
        if ((errno == ERANGE && (val == HUGE_VAL || val == -HUGE_VAL))
            || (errno != 0 && val == 0.0) || (next == str)) {
          if (errno)
            perror("PPM_sscana");
          else
            fputs("PPM_sscana: parse error\n", stderr);
          *ierror = parse_error;
          return (ssize_t) MIN(i, (size_t)SSIZE_MAX);
        }
        *(double *)(out[i].addr) = val;
        str = next;
      }
      break;
    default:
      PPM_abort(PPM_default_comm, "invalid parse type", __FILE__, __LINE__);
    }
  }
  *ierror = no_error;
  return (ssize_t)MIN(i, (size_t)SSIZE_MAX);
}


static int32_t
PPM_sscana_f(const char *str, struct PPM_fmt_elem out[],
             int32_t out_size, int32_t *ierror)
{
  ssize_t rv;
  int ierr;
  assert(out_size > 0);
  rv = PPM_sscana(str, out, (size_t)out_size, &ierr);
  *ierror = ierr;
  return (int32_t)rv;
}

FCALLSCFUN4(INT,PPM_sscana_f,PPM_SSCANA,ppm_sscana,STRING,PVOID,INT,PINT)

#if SIZEOF_INT_P == 8
static inline void
PPM_get_address(void *p, long long *a)
{
  *(int64_t *)a = (int64_t)p;
}

FCALLSCSUB2(PPM_get_address,PPM_GET_ADDRESS,ppm_get_address,PVOID,PLONGLONG)
#elif SIZEOF_INT_P == 4
static inline void
PPM_get_address(void *p, int32_t *a)
{
  *a = (int32_t)p;
}

FCALLSCSUB2(PPM_get_address,PPM_GET_ADDRESS,ppm_get_address,PVOID,PINT)
#else
#error "pointer size is unknown"
#endif
/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
