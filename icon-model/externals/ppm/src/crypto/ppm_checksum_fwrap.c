/**
 * @file ppm_checksum_fwrap.c
 * @brief message digest generation Fortran/C interface
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
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

#include <inttypes.h>

#include "cfortran.h"
#include "crypto/ppm_checksum.h"

static inline void
PPM_describe_digest_f(int digest, struct PPM_digest_description *hd)
{
  PPM_describe_digest((enum digest_type)digest, hd);
}

FCALLSCSUB2(PPM_describe_digest_f, PPM_DESCRIBE_DIGEST, ppm_describe_digest,
            INT,PVOID)

static inline void
PPM_hex_checksum_f(const void *buf, int buf_len, int elem_size,
                   struct PPM_digest_description *digest,
                   char *hexchecksum)
{
  unsigned char checksum[digest->size];
  PPM_checksum(buf, (size_t)buf_len * (size_t)elem_size, checksum, digest);
  PPM_md2hex(hexchecksum, checksum, (size_t)digest->size);
}

FCALLSCSUB5(PPM_hex_checksum_f,PPM_HEX_CHECKSUM_F,ppm_hex_checksum_f,
            PVOID,INT,INT,PVOID,PSTRING)

static inline void
PPM_hex_checksum_str_f(const char *buf, int buf_len,
                       struct PPM_digest_description *digest,
                       char *hexchecksum)
{
  unsigned char checksum[digest->size];
  PPM_checksum(buf, (size_t)buf_len, checksum, digest);
  PPM_md2hex(hexchecksum, checksum, (size_t)digest->size);
}

FCALLSCSUB4(PPM_hex_checksum_str_f,PPM_HEX_CHECKSUM_STR_F,
            ppm_hex_checksum_str_f,STRING,INT,PVOID,PSTRING)

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
