/**
 * @file ppm_checksum.h
 * @brief checksum computations interface description
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

enum digest_type {
  digest_type_md5 = 1,
  digest_type_sha1 = 2,
};

enum digest_family {
  digest_none = 0,
  digest_crypto = 1,
  digest_internal = 2,
};

struct PPM_digest_description
{
  void *helper;
  int32_t size, family;
};

/**
 * Fill in digest descriptor for desired message digest function.
 *
 * @param[in] digest_type desired function
 * @param[out] hd descriptor to fill
 */
void
PPM_describe_digest(enum digest_type digest, struct PPM_digest_description *hd);

/**
 * Compute checksum on data
 *
 * @param[in] buf data to compute checksum of
 * @param[in] buf_size number of octets in buf
 * @param[out] checksum binary output of checksum,
 * must provide for enough storage according to digest
 * @param[in] digest descriptor for desired digest
 */
void
PPM_checksum(const void *buf, size_t buf_size, unsigned char *checksum,
             struct PPM_digest_description *digest);

/**
 * Convert binary to hexadecimal string.
 *
 * @param[out] hex string to write hex transliteration of checksum to
 * @param[in] checksum binary string
 * @param[in] checksum_size number of octets in checksum
 */
void
PPM_md2hex(char *hex, const unsigned char *checksum, size_t checksum_size);

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
