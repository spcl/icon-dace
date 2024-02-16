/**
 * @file ppm_extents_c.c
 * @brief extent implementation for C bindings
 *
 * Copyright  (C)  2010-2017  Thomas Jahns <jahns@dkrz.de>
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
 *
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <float.h>
#include <stdio.h>

#include "ppm_extents.h"

void
PPM_extents2iintervals(int ndims, struct PPM_iinterval dst[ndims],
                       const struct PPM_extent src[ndims])
{
  for (int i = 0; i < ndims; ++i)
    dst[i] = PPM_extent2iinterval(src[i]);
}

void
PPM_extents2iintervals64(int ndims, struct PPM_iinterval64 dst[ndims],
                         const struct PPM_extent64 src[ndims])
{
  for (int i = 0; i < ndims; ++i)
    dst[i] = PPM_extent2iinterval64(src[i]);
}



int
PPM_sprint_extent(char buf[], const struct PPM_extent *ext)
{
  int len;
  if (ext->size != 0)
    len = sprintf(buf, "[%"PRId32",%"PRId32"]", ext->first,
                  ext->first + ext->size - (ext->size<0?-1:1));
  else
  {
    buf[0] = '{';
    buf[1] = '}';
    buf[2] = '\0';
    len = 2;
  }
  return len;
}

int
PPM_sprint_extent64(char buf[], const struct PPM_extent64 *ext)
{
  int len;
  if (ext->size != 0)
    len = sprintf(buf, "[%"PRId64",%"PRId64"]", ext->first,
                  ext->first + ext->size - (ext->size<0?-1:1));
  else
  {
    buf[0] = '{';
    buf[1] = '}';
    buf[2] = '\0';
    len = 2;
  }
  return len;
}

int
PPM_sprint_iinterval(char buf[], const struct PPM_iinterval *iinterval)
{
  return sprintf(buf, "[%"PRId32",%"PRId32"]",
                 iinterval->first, iinterval->last);
}

int
PPM_sprint_iinterval64(char buf[], const struct PPM_iinterval64 *iinterval)
{
  return sprintf(buf, "[%"PRId64",%"PRId64"]",
                 iinterval->first, iinterval->last);
}

int
PPM_sprint_iinterval_sp(char buf[], const struct PPM_iinterval_sp *iinterval)
{
  return sprintf(buf, "[%.*g,%.*g]",
                 PPM_FLT_DECIMAL_DIG, iinterval->first,
                 PPM_FLT_DECIMAL_DIG, iinterval->last);
}

int
PPM_sprint_iinterval_dp(char buf[], const struct PPM_iinterval_dp *iinterval)
{
  return sprintf(buf, "[%.*g,%.*g]",
                 PPM_DBL_DECIMAL_DIG, iinterval->first,
                 PPM_DBL_DECIMAL_DIG, iinterval->last);
}


/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
